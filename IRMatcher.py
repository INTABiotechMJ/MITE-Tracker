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
parser.add_argument("--max_sep_len", help="MITE max lenght", type=int, default=650)
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


MIN_TSD_LEN = args.tsd_min_len
MAX_TSD_LEN = args.tsd_max_len
MITE_MAX_LEN = args.mite_max_len
MITE_MIN_LEN = args.mite_min_len
TIR_ALIGN_MIN_LEN = args.tir_align_min_len

def _parse(q):
    global last_porc
    global total_queue_count
    global mites
    while True:
        try:
            seq, split_index, record_id = q.get(timeout=10)
        except Queue.Empty:
            q.task_done()
            break
        splited_len = len(seq)
        if lcc_simp(seq) < 0.6: #Discard really low complexity MITE
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
        cmd = " ".join(cmd_list)
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
        out,err = p.communicate()
        #blast_process = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        #out,err = blast_process.communicate()
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
            if length < TIR_ALIGN_MIN_LEN:
                continue

            #subject transform cause it was reversed
            transform_send = splited_len - sstart 
            transform_sstart = splited_len - send 

            #obtain IR sequences
            seq_q = seq[qstart:qend]
            seq_q_prime = seq[transform_sstart:transform_send]
            seq_s = Seq(seq_q_prime).reverse_complement()

            #organice positions
            #subject further than query
            if (transform_send + transform_sstart) / 2 >  (qend + qstart) / 2:
                ir_1_start = qstart
                ir_1_end = qend
                ir_2_start = transform_sstart
                ir_2_end = transform_send
            else:
                ir_1_start = transform_sstart
                ir_1_end = transform_send
                ir_2_start = qstart
                ir_2_end = qend
            #validate distance between IR
            #print "-----"
            #print ir_1_start, ir_1_end
            #print ir_2_start, ir_2_end
            #print seq_q, seq_s
            if ir_2_end - ir_1_start > MITE_MAX_LEN:
             #   print "discarded by far"
                continue
            if ir_2_end - ir_1_start < MITE_MIN_LEN:
              #  print "discarded by close"
                continue

            #again validate complexity
            # a value of 1 means only two different nucleotides are present
            if lcc_simp(seq_q) <= 1.3:
               # print "discarded by low complexity"
                continue

            #validate TSD outside TIRs
            i = MAX_TSD_LEN
            valid_tsd = False
            while i >= MIN_TSD_LEN:
                tsd_one = seq[ir_1_start - i:ir_1_start]
                tsd_two = seq[ir_2_end:ir_2_end + i]
                if tsd_one.lower() == tsd_two.lower():
                    valid_tsd = True
                    mite_pos_one = ir_1_start - i
                    mite_pos_two = ir_2_end + i
                    tsd_in = 'no'
                #    print "valid"
                    break
                i -= 1
            #validate TSD inside TIRs
            #TSDs cannot be a large part of TIRs
            if not valid_tsd:
                i = MAX_TSD_LEN
                while i >= MIN_TSD_LEN:
                    tsd_one = seq[ir_1_start:ir_1_start+i]
                    tsd_two = seq[ir_2_end-i:ir_2_end]
                    if tsd_one.lower() == tsd_two.lower():
                        valid_tsd = True
                        mite_pos_one = ir_1_start
                        mite_pos_two = ir_2_end 
                        tsd_in = 'yes'
                 #       print "valid"
                        break
                    i -= 1
            if not valid_tsd:
                #print "no tsd"
                continue
            mite = seq[mite_pos_one:mite_pos_two]
            mite_len = mite_pos_two - mite_pos_one

            #calculate positions in full sequence
            mite_start_full = mite_pos_one + split_index
            mite_end_full = mite_pos_two + split_index 
            with l_lock:
                mite = {'mite': mite, 'id':record.id, 
                        'start':mite_start_full, 'end':mite_end_full, 
                        'len':mite_len, 'tsd':tsd_one, 'tsd_in':tsd_in
                        } 
                mites.append(mite)
        porc = (total_queue_count - q.unfinished_tasks) * 100 / total_queue_count
        if porc - last_porc >= 10:
            print '%i%% ' % porc,
            last_porc = porc
            sys.stdout.flush()
        q.task_done()
#-soft_masking false -dust no
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
mites = []
l_lock = Lock()
count = 1
record_count = 0

q = Queue.Queue(maxsize=0)
for i in range(args.workers):
    worker = Thread(target=_parse, args=(q,))
    worker.setDaemon(True)
    worker.start()
windows_size = int(ceil(MITE_MAX_LEN * 2))

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
    split_index = MAX_TSD_LEN
    clean_seq = ''.join(str(record.seq).splitlines())
    seq_len = len(clean_seq)

    while split_index < seq_len - MITE_MIN_LEN:
        seq = clean_seq[split_index:split_index + windows_size]
        q.put((seq, split_index,record.id,))
        queue_count += 1
        total_queue_count += 1
        split_index += MITE_MAX_LEN
        current_processing_size += MITE_MAX_LEN
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
mites_arr = []
gff_buff = []

#delete duplicated
id_1 = 1
unique_mites = []
for mite in mites:
    mites_2 = mites[id_1:]
    duplicated = False
    for mite_2 in mites_2:
        if mite['id'] == mite_2['id'] and mite['start'] == mite_2['start'] and mite_2['end'] == mite['end']: 
            duplicated = True
    if not duplicated:
        unique_mites.append(mite)
    id_1 += 1
mites = unique_mites

#delete nested
id_1 = 0
for mite in mites:
    id_1 +=1
    nested = False
    id_2 = 0
    for mite_2 in mites:
        id_2 += 1
        if id_1 == id_2:
            continue
        if mite['id'] == mite_2['id'] and mite['start'] >= mite_2['start'] and mite_2['end'] >= mite['end']: 
            #mite is nested in mite_2
            nested = True
    if not nested:
        params = (mite['id'], mite['start'], mite['end'], mite['len'], mite['tsd'], mite['tsd_in'], mite['tir_1'], mite['tir_2'] )
        description = "SEQ:%s START:%i END:%i MITE_LEN:%i TSD:%s TSD_IN:%s TIR_1:%s TIR_2:%s" % (params)
        mite_seq_rec = SeqRecord(Seq(mite['mite']), id='MITE_' + str(count), description = description)
        mites_arr.append( mite_seq_rec)
        mite_new = {'from': mite['id'], 'id':'MITE_' + str(count),'start':mite['start'], 'end':mite['end']}
        gff_buff.append(mite_new) 
        count += 1

#gff
output_gff = open("results/" + args.jobname + "/mites.gff3","w") 
output_gff.write("##gff-version 3\n")
for mite in gff_buff:
    write_row =  '\t'.join([ mite['from'], 'miteParser','MITE',str(mite['start']), str(mite['end']),'.','+','.','ID='+mite['id'] ]) 
    output_gff.write(write_row + '\n')

SeqIO.write(mites_arr, "results/" + args.jobname + "/mites.fasta" , "fasta")
print ""
makelog("Found %i MITES in" % (count - 1,))
makelog(cur_time())