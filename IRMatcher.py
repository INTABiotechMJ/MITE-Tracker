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
    global last_porc
    global total_queue_count
    global irs
    while True:
        try:
            seq, split_index, record_id = q.get(timeout=10)
        except Queue.Empty:
            q.task_done()
            break
        splited_len = len(seq)
        if lcc_simp(seq) < 0.6: #Discard really low complexity IR
            q.task_done()
            continue
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
        '-evalue','1','-strand',"plus",
        #'-soft_masking','false' ,'-dust','no',
        '-outfmt',"'6 sstart send qstart qend score length mismatch gaps gapopen nident'"]
        cmd = ' '.join(cmd_list)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
        out,err = p.communicate()
        #blast_process = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        #out,err = blast_process.communicate()
        os.remove(query_filename)
        os.remove(subject_filename)
        lines = out.splitlines()
        found = 0
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
            seq_q_prime = seq[sstart:send]
            seq_s = Seq(seq_q_prime).reverse_complement()

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
            if lcc_simp(seq_q) <= 1.3:
                continue
            with l_lock:
                ir = {'ir': ir_seq, 'id':record.id, 
                        'start':ir_start, 'end':ir_end, 
                        'len':ir_len,
                        } 
                irs.append(ir)
                found += 1
        porc = (total_queue_count - q.unfinished_tasks) * 100 / total_queue_count
        if porc - last_porc >= 10:
            #print '%i%% ' % porc,
            last_porc = porc
            sys.stdout.flush()
        q.task_done()

start_time = time.time()
makelog("Counting sequences: ")
fh = open(args.genome)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
seqs_count = n
makelog(n)
fh.close()

fasta_seq = SeqIO.parse(args.genome, 'fasta')
irs = []
l_lock = Lock()
count = 1
record_count = 0

q = Queue.Queue(maxsize=0)
for i in range(args.workers):
    worker = Thread(target=_findIR, args=(q,))
    worker.setDaemon(True)
    worker.start()
windows_size = int(ceil(max_sep_len * 2))

#processes until certain amount of sequences
#stablish a balance between memory usage and processing
max_processing_size = windows_size * 30
current_processing_size = 0

for record in fasta_seq:
    total_queue_count = 1
    queue_count = 1
    last_porc = -5
    record_count += 1
    porc_ant = 0
    split_index = 0
    clean_seq = ''.join(str(record.seq).splitlines())
    seq_len = len(clean_seq)
    while split_index < seq_len - min_total_len:
        seq = clean_seq[split_index:split_index + windows_size]
        q.put((seq, split_index,record.id,))
        queue_count += 1
        total_queue_count += 1
        split_index += max_sep_len
        current_processing_size += max_sep_len
    #print ""
    params = (record.id, record_count , seqs_count, (record_count * 100 / seqs_count), cur_time() )
    makelog("Adding %s %i out of %i(%i%% of total in %s)" % params)
    #in order to avoid overloading of memory, we add a join()
    #we do not direcly use the join() method so we can process several
    #small sequences at once
    if current_processing_size >= max_processing_size:
        current_processing_size = 0
        q.join()
#In case of unprocessed sequences are left, let's wait
q.join()

#a centinell strategy, unused for now
#for i in range(args.workers):
#    q.put((-1, -1))

count = 1
ir_arr = []
gff_buff = []

#delete duplicated
id_1 = 1
unique_ir = []
for ir in irs:
    ir_2s = irs[id_1:]
    duplicated = False
    for ir_2 in ir_2s:
        if ir['id'] == ir_2['id'] and ir['start'] == ir_2['start'] and ir_2['end'] == ir['end']: 
            duplicated = True
    if not duplicated:
        unique_ir.append(ir)
    id_1 += 1
irs = unique_ir


#delete nested
id_1 = 0
for ir in irs:
    id_1 +=1
    nested = False
    id_2 = 0
    for ir_2 in irs:
        id_2 += 1
        if id_1 == id_2:
            continue
        if ir['id'] == ir_2['id'] and ir['start'] >= ir_2['start'] and ir_2['end'] >= ir['end']: 
            #ir is nested in ir_2
            nested = True
    if not nested:
        params = (ir['id'], ir['start'], ir['end'], ir['len'])
        description = "SEQ:%s START:%i END:%i ir_LEN:%i" % (params)
        ir_seq_rec = SeqRecord(Seq(ir['ir']), id='ir_' + str(count), description = description)
        ir_arr.append( ir_seq_rec)
        ir_new = {'from': ir['id'], 'id':'ir_' + str(count),'start':ir['start'], 'end':ir['end']}
        gff_buff.append(ir_new) 
        count += 1

#gff
output_gff = open("results/" + args.jobname + "/ir.gff3","w") 
output_gff.write("##gff-version 3\n")
for ir in gff_buff:
    write_row =  '\t'.join([ ir['from'], 'irParser','ir',str(ir['start']), str(ir['end']),'.','+','.','ID='+ir['id'] ]) 
    output_gff.write(write_row + '\n')

SeqIO.write(ir_arr, "results/" + args.jobname + "/ir.fasta" , "fasta")
print ""
makelog("Found %i ir in" % (count - 1,))
makelog(cur_time())