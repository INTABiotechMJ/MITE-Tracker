# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import Popen, PIPE
import os
import time
import argparse
from math import ceil
import sys
import pandas as pd
from threading import Thread, Lock, active_count
import Queue
import logging

def makelog(stri):
    print(stri)
    logging.debug(stri)

def cur_time():
    global start_time
    elapsed_time = time.time() - start_time
    return "%f secs" % (elapsed_time,)

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genome", help="Genome file in fasta format", required=True)
parser.add_argument("-j","--jobname", help="Will create files under a folder called [jobname]", required=True)
parser.add_argument("-w","--workers", help="Max number of processes to use simultaneously", type=int, default=1)
parser.add_argument("--max_sep_len", help="IR max separation lenght", type=int, default=650)
parser.add_argument("--min_total_len", help="Min total lenght", type=int, default=50)
parser.add_argument("--align_min_len", help="TIR minimun aligmnent length", type=int, default=10)
args = parser.parse_args()#pylint: disable=invalid-name

#write results
if not os.path.isdir("results/" + args.jobname):
    os.mkdir("results/" + args.jobname)

filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "results/" + args.jobname + "/out.log")
logging.basicConfig(
    filename=filename,
    level=logging.DEBUG,
    handlers=[
    logging.FileHandler(filename),
    logging.StreamHandler()],
    format="%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")

max_sep_len = args.max_sep_len
min_total_len = args.min_total_len
align_min_len = args.align_min_len

def _findIR(q):
    global total_queue_count
    global intersecter
    global irs
    while True:
        try:
            seq, split_index, record_id,seq_len = q.get(timeout=5)
        except Queue.Empty:
            break
        splited_len = len(seq)
        seq_rc = str(Seq(seq).reverse_complement())
        record_q = SeqRecord(Seq(seq), id = record_id)
        record_s = SeqRecord(Seq(seq_rc), id = record_id + "_rc")
        query_filename = "tmp/query" + str(record_id + "_" + str(split_index))+".tmp"
        subject_filename = "tmp/subject" + str(record_id + "_" + str(split_index))+".tmp"
        SeqIO.write(record_q, query_filename, "fasta")
        SeqIO.write(record_s, subject_filename, "fasta")
        cmd_list = [
        'blastn',
        '-query',query_filename,
        '-subject',subject_filename,
        '-reward','2',
        '-max_target_seqs','1',
        '-penalty','-4',
        '-word_size','7',#'-ungapped',
        '-evalue','150','-strand',"plus",
        '-soft_masking','false' ,'-dust','no',
        '-outfmt',"'6 sstart send qstart qend score length mismatch gaps gapopen nident'"]
        cmd = ' '.join(cmd_list)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
        out,err = p.communicate()
        if err:
            q.task_done
            makelog("BLASTN error: %s" % (err, ) )
        os.remove(query_filename)
        os.remove(subject_filename)
        lines = out.splitlines()
        for row in lines:
            row = row.split()
            sstart = int(row[0])
            send = int(row[1])
            qstart = int(row[2])
            qend = int(row[3])
            score = int(row[4])
            length = int(row[5])
            mismatch = int(row[6])
            gaps = int(row[7])
            #filter valids IR
            if length < align_min_len:
                continue
            #subject transform cause it was reversed
            sstart = splited_len - sstart 
            send = splited_len - send 
            #obtain IR sequences
            seq_q = seq[qstart:qend]
            seq_q_prime = seq[send:sstart]
            #organice positions
            ir_start = min(qstart,qend,sstart,send)
            ir_end = max(qstart,qend,sstart,send)
            ir_len = ir_end - ir_start
            #length constraints
            if ir_len > max_sep_len:
                continue
            if ir_len < min_total_len:
                continue
            #move in genome, split index
            ir_seq = seq[ir_start:ir_end]
            ir_start += split_index
            ir_end += split_index
            #again validate complexity, a value of 1 means only two different nucleotides are present
            if lcc_simp(seq_q) <= 1.2:
                continue
            new_element = (ir_start, ir_end,ir_seq, record.id, ir_len, seq_q, seq_q_prime)
            with l_lock:
                irs.append(new_element)
        q.task_done()

start_time = time.time()
#Count sequences
makelog("Counting sequences: ")
fh = open(args.genome)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
seqs_count = n
makelog(n)
fh.close()

#initialize workers
q = Queue.Queue(maxsize=0)
for i in range(args.workers):
    worker = Thread(target=_findIR, args=(q,))
    worker.setDaemon(True)
    worker.start()
windows_size = 20000 #int(ceil(max_sep_len * 30))

#processes until certain amount of sequences
#stablish a balance between memory usage and processing
max_queue_size = 50
current_processing_size = 0
#initialize global variables
irs = []
l_lock = Lock()
#start adding sequences to process queue
count = 1
record_count = 0
fasta_seq = SeqIO.parse(args.genome, 'fasta')
for record in fasta_seq:
    processed = False
    total_queue_count = 1
    queue_count = 1
    record_count += 1
    porc_ant = 0
    split_index = 0
    clean_seq = ''.join(str(record.seq).splitlines())
    seq_len = len(clean_seq)
    params = (record.id, seq_len, record_count , seqs_count, (record_count * 100 / seqs_count), cur_time() )
    print ""
    makelog("Adding %s (len %i) %i/%i (%i%% of total sequences in %s)" % params)
    while split_index < seq_len - min_total_len:
        seq = clean_seq[split_index:split_index + windows_size]
        q.put((seq, split_index,record.id,seq_len,))
        queue_count += 1
        total_queue_count += 1
        split_index += windows_size - max_sep_len
        current_processing_size += windows_size
        if q.qsize() >= max_queue_size:
            current_processing_size = 0
            q.join()
            processed = True
        #in order to avoid overloading of memory, we add a join()
        #we do not direcly use the join() method so we can process several
        #small sequences at once

#In case of unprocessed sequences are left, let's wait
q.join()

makelog("Creating gff and fasta")
output_gff = open("results/" + args.jobname + "/IR.gff3","w") 
output_gff.write("##gff-version 3\n")

labels = ['start','end','seq','record','len','ir_1','ir_2']
df = pd.DataFrame.from_records(irs, columns=labels)

#filter out nested (keep larger)
l=[]
for idx, row in df.iterrows():
    res = df[(df.index != idx) & (df.start >= row.start) & (df.end <= row.end)]
    l.append(res)
res = pd.concat(l)
df.drop(res.index,inplace=True)

irs_seqs = []
df = df.sort_values(by=['start','end'])
for _, row in df.iterrows():
    count += 1
    name = 'IR_' + str(count)
    #append sequence record for biopython
    params = (row.record, row.start, row.end, row.len, row.ir_1, row.ir_2)
    description = "SEQ:%s START:%i END:%i IR_LEN:%i IR_1:%s IR_2:%s " % (params)
    ir_seq_rec = SeqRecord(Seq(row.seq), id=name, description = description)
    irs_seqs.append(ir_seq_rec)
    write_row = '\t'.join([ row.record, 'IRMatcher','inverted_repeat',str(row.start), str(row.end),'.','+','.','ID='+name ])
    output_gff.write(write_row + '\n')

makelog("Writing fasta")
SeqIO.write(irs_seqs, "results/" + args.jobname + "/IR.fasta" , "fasta")
print ""
makelog("Found %i inverted repeats in" % (count - 1,))
makelog(cur_time())