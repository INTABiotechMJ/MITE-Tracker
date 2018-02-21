# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqRecord import SeqRecord
from itertools import chain
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

def makelog(stri, do_print=True):
    if do_print:
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
parser.add_argument("--max_sep_len", help="IR max separation lenght", type=int, default=800)
parser.add_argument("--min_total_len", help="Min total lenght", type=int, default=100)
parser.add_argument("--align_min_len", help="TIR minimun aligmnent length", type=int, default=10)
parser.add_argument("--tsd_min_len", help="TSD min lenght", type=int, default=2)
parser.add_argument("--tsd_max_len", help="TSD max lenght", type=int, default=10)
parser.add_argument("--FSL", help="Flanking seq length for comparison", type=int, default=50)
parser.add_argument("--min_copy_number", help="Minimum CN for families", type=int, default=3)

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
MIN_TSD_LEN = args.tsd_min_len
MAX_TSD_LEN = args.tsd_max_len
FSL = args.FSL
MIN_COPY_NUMBER = args.min_copy_number

def _findIR(q):
    global total_queue_count
    global intersecter
    global irs
    global flanking_seqs
    while True:
        try:
            seq, seq_fs, split_index, record_id,seq_len = q.get(timeout=5)
        except Queue.Empty:
            break
        splited_len = len(seq)
        seq_rc = str(Seq(seq).reverse_complement())
        complexity = lcc_simp(seq)
        if complexity < 1:
            q.task_done()
            continue
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
        #'-max_target_seqs','1',
        '-penalty','-4',
        '-word_size','7',
        #'-ungapped',
        '-evalue','145',
        '-strand',"plus",
        #'-soft_masking','false' ,'-dust','no',
        '-outfmt',"6 sstart send qstart qend score length mismatch gaps gapopen nident"]
        cmd = ' '.join(cmd_list)
        p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        out,err = p.communicate()
        if err:
            print(split_index, record_id,seq_len)
            makelog("BLASTN error: %s" % (err, ) )
        else:
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
            #ir_seq = seq[ir_start:ir_end]

            #again validate complexity, a value of 1 means only two different nucleotides are present
            if lcc_simp(seq_q) <= 1.2:
                continue

            #validate TSD outside TIRs
            i = MAX_TSD_LEN
            valid_tsd = False
            while i >= MIN_TSD_LEN:
                tsd_one = seq_fs[ir_start - i + FSL:ir_start + FSL]
                tsd_two = seq_fs[ir_end + FSL:ir_end + i + FSL]
                if tsd_one.lower() == tsd_two.lower():
                    valid_tsd = True
                    mite_pos_one = ir_start - i
                    mite_pos_two = ir_end + i
                    tsd_in = 'no'
                    break
                i -= 1
            #validate TSD inside TIRs
            #TSDs cannot be a large part of TIRs
            if not valid_tsd:
                i = MAX_TSD_LEN
                while i >= MIN_TSD_LEN:
                    tsd_one = seq_fs[ir_start + FSL:ir_start + i + FSL]
                    tsd_two = seq_fs[ir_end - i + FSL:ir_end + FSL]
                    if tsd_one.lower() == tsd_two.lower():
                        valid_tsd = True
                        mite_pos_one = ir_start
                        mite_pos_two = ir_end 
                        tsd_in = 'yes'
                        break
                    i -= 1
            #"no tsd"
            if not valid_tsd:
                continue

            ir_seq = seq[mite_pos_one:mite_pos_two]
            ir_len = mite_pos_two - mite_pos_one

            flanking_seq_left = seq_fs[mite_pos_one:mite_pos_one + FSL]
            flanking_seq_right = seq_fs[mite_pos_two+FSL:mite_pos_two + FSL + FSL]
           
            #calculate positions in full sequence
            mite_start_full = mite_pos_one + split_index
            mite_end_full = mite_pos_two + split_index 

            new_element = (mite_start_full, mite_end_full, ir_seq, record.id, ir_len, seq_q, seq_q_prime, tsd_one, tsd_in,flanking_seq_left,flanking_seq_right,length)
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
windows_size = 25000 #int(ceil(max_sep_len * 30))

#processes until certain amount of sequences
#stablish a balance between memory usage and processing
max_queue_size = 50
current_processing_size = 0
#initialize global variables
irs = []
flanking_seqs = []
l_lock = Lock()
#start adding sequences to process queue
record_count = 0
fasta_seq = SeqIO.parse(args.genome, 'fasta')
for record in fasta_seq:
    processed = False
    total_queue_count = 1
    queue_count = 1
    record_count += 1
    porc_ant = 0
    split_index = MAX_TSD_LEN + FSL
    clean_seq = ''.join(str(record.seq).splitlines())
    seq_len = len(clean_seq)
    params = (record.id, seq_len, record_count , seqs_count, (record_count * 100 / seqs_count), cur_time() )
    print ""
    makelog("Adding %s (len %i) %i/%i (%i%% of total sequences in %s)" % params)
    while split_index < seq_len - MIN_TSD_LEN - FSL:
        seq = clean_seq[split_index:split_index + windows_size]
        seq_fs = clean_seq[split_index-FSL :split_index + windows_size + FSL]
        q.put((seq, seq_fs, split_index,record.id,seq_len,))
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

makelog("Creating candidates fasta")

labels = ['start','end','seq','record','len','ir_1','ir_2','tsd','tsd_in','fs_left','fs_right', 'ir_length']
df = pd.DataFrame.from_records(irs, columns=labels)

makelog("Canidadates: " + str(len(df)))
#filter out nested (keep larger)
l=[]
for idx, row in df.iterrows():
    #filter all that are nested into this
    res = df[(df.record == row.record) & (df.start >= row.start + row.ir_length) & (df.end <= row.end - row.ir_length) & (df.index != idx) ]
    #res = df[(df.record == row.record) & (df.start >= row.start) & (df.end <= row.end ) & (df.index != idx) ]
    l.append(res)
res = pd.concat(l)
df.drop(res.index,inplace=True)
makelog("Canidadates (not nested): " + str(len(df)))
#remove all nested
#print len(df)
#df1 = df.reset_index()
#df1 = pd.merge(df1, df1, on='record')
#m = ((df1.start_x >= df1.start_y) & (df1.end_x <= df1.end_y) & (df1.index_x != df1.index_y))
#idx = df1.loc[m, 'index_x']
#df = df.drop(idx)
#print len(df)

fs_seqs = []
irs_seqs = []
df = df.sort_values(by=['record','start','end'])
count = 1
for _, row in df.iterrows():
    name = 'MITE_CAND_' + str(count)
    #append sequence record for biopython
    params = (row.record, row.start, row.end, row.tsd, row.tsd_in, row.len, row.ir_1, row.ir_2)
    description = "SEQ:%s START:%i END:%i TSD:%s TSD_IN:%s MITE_LEN:%i IR_1:%s IR_2:%s " % (params)
    ir_seq_rec = SeqRecord(Seq(row.seq), id=name, description = description)
    irs_seqs.append(ir_seq_rec)

    #flanking sequences
    fs_seq_rec = SeqRecord(Seq(row.fs_left), id=name + "L", description = "_") 
    fs_seqs.append(fs_seq_rec)
    #flanking sequences
    fs_seq_rec = SeqRecord(Seq(row.fs_right), id=name + "R", description = "_") 
    fs_seqs.append(fs_seq_rec)

    count += 1
makelog("Writing candidates sequences")
candidates_fasta = "results/" + args.jobname + "/mites.candidates.fasta"
SeqIO.write(irs_seqs, candidates_fasta , "fasta")

flanking_seqs_name = "results/" + args.jobname + "/flanking_seqs.candidates.fasta"
SeqIO.write(fs_seqs, flanking_seqs_name , "fasta")

makelog("Group elements")

cmd_list = [
'blastn',
'-query',candidates_fasta,
'-subject',candidates_fasta,
'-evalue','1e10',
'-outfmt',"6"]
p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
out,err = p.communicate()
if err:
    q.task_done
    makelog("BLASTN error: %s" % (err, ) )
lines = out.splitlines()

families = []
for row in lines:
    row = row.split()
    query = row[0]
    subject = row[1]
    identity = row[2]
    length = row[3]
    if query == subject:
        continue
    if float(identity < 80):
        continue
    if float(length < 80):
        continue
    res = []
    for family in families:
        if query in family or subject in family:
            res.append(family)
    if len(res) == 0:
        new_set = set([query, subject])
        families.append(new_set)
    elif len(res) == 1:
        res[0].add(query)
        res[0].add(subject)
    else:
        for r in res:
            try:
                families.remove(r)
            except ValueError:
                print ">>",families
                print "-->",res
                print "->",r
                exit()
        new_set = set(chain.from_iterable(res))
        families.append(new_set)

#group flanking sequences
makelog("Group flanking sequence")

cmd_list = [
'blastn',
'-query',flanking_seqs_name,
'-subject',flanking_seqs_name,
'-evalue','10',
'-outfmt',"6"]
p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
out,err = p.communicate()
if err:
    q.task_done
    makelog("BLASTN error: %s" % (err, ) )
lines = out.splitlines()
fs_families = []

for row in lines:
    row = row.split()
    query = row[0]
    subject = row[1]
    if query == subject:
        continue
    res = []
    for family in fs_families:
        if query in family or subject in family:
            res.append(family)
    if len(res) == 0:
        new_set = set([query[:-1], subject[:-1]])
        fs_families.append(new_set)
    elif len(res) == 1:
        res[0].add(query[:-1])
        res[0].add(subject[:-1])
    else:
        for r in res:
            try:
                fs_families.remove(r)
            except ValueError:
                print ">>",fs_families
                print "-->",res
                print "->",r
                exit()
        new_set = set(chain.from_iterable(res))
        fs_families.append(new_set)


makelog("Merging results")
#shared elements are in the same family and have the same flanking sequence
shared = set()
for family in families:
    for fs_family in fs_families:
        shared = shared | family.intersection(fs_family)

#remove families with low CN
for family in families[:]:
    family -= shared
    if len(family) < MIN_COPY_NUMBER:
        families.remove(family)

makelog("Saving results")
#save definitive elements
families = list(families)
irs_seqs = []
done_families = []
family_seqs = []
current_family = None
df = df.sort_values(by=['record','start','end'])
count = 1
count_real = 1
output_gff = open("results/" + args.jobname + "/mites.gff3","w") 
output_gff.write("##gff-version 3\n")
no_family, flank_seq_sim = 0, 0
for index, row in df.iterrows():
    if "MITE_CAND_" + str(count) in shared:
        count += 1
        flank_seq_sim += 1
        makelog("MITE_CAND_" + str(count) + " discarded because flanking sequence and inner sequence similarity", False)
        continue
    count += 1
    
    idx = 0
    family_number = 0
    for family in families:
        idx += 1
        if "MITE_CAND_" + str(count) in family:
            family_number = idx
            break

    #are not in any family
    if family_number == 0:
        no_family += 1
        makelog("MITE_CAND_" + str(count) + " discarded because is not in any family", False)
        continue
    df.loc[index,'family'] = int(family_number)
    df.loc[index,'id'] = "MITE_" + str(count_real)
    count_real += 1

for _, row in df.iterrows():
    if str(row.family) == 'nan':
        continue
    #append sequence record for biopython
    params = (row.record, row.start, row.end, int(row.family) , row.tsd, row.tsd_in, row.len, row.ir_1, row.ir_2)
    description = "SEQ:%s START:%i END:%i FAMILY:%s TSD:%s TSD_IN:%s MITE_LEN:%i IR_1:%s IR_2:%s " % (params)
    ir_seq_rec = SeqRecord(Seq(row.seq), id=row.id, description = description)
    irs_seqs.append(ir_seq_rec)

    if not row.family in done_families:
        done_families.append(row.family)
        params = (','.join(df[df.family == row.family].id))
        description += "ELEMENTS_IN_FAMILY: %s" % (params)
        ir_seq_rec = SeqRecord(Seq(row.seq), id=row.id, description = description)
        family_seqs.append(ir_seq_rec)

    write_row = '\t'.join([ row.record, 'miteParser','mite',str(row.start), str(row.end),'.','+','.','ID='+row.id+';FAMILY='+str(int(row.family))])
    output_gff.write(write_row + '\n')
    count_real += 1

makelog("Writing definitive sequences")
SeqIO.write(irs_seqs, "results/" + args.jobname + "/mites.all.fasta" , "fasta")
SeqIO.write(family_seqs, "results/" + args.jobname + "/mites.nr.fasta" , "fasta")

makelog("Found %i MITES and %i families MITEs in" % (count_real - 1,len(families),))
makelog("Discarded %i because they're in no family" % (no_family,) )
makelog("Discarded %i because flanking sequence similarity" % (flank_seq_sim,) )
makelog(cur_time())