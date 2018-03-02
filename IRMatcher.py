# -*- coding: utf-8 -*-
from __future__ import division
from Bio import SeqIO
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import itertools
from Bio.Seq import Seq
from subprocess import Popen, PIPE
import os
import time
import argparse
from math import ceil
import clusterutils
from collections import OrderedDict
import sys
import pandas as pd
from threading import Thread, Lock, active_count
import Queue
import logging
import findir

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genome", help="Genome file in fasta format", required=True)
parser.add_argument("-j","--jobname", help="Will create files under a folder called [jobname]", required=True)
parser.add_argument("-w","--workers", help="Max number of processes to use simultaneously", type=int, default=1)
parser.add_argument("--mite_max_len", help="MITE max lenght", type=int, default=800)
parser.add_argument("--mite_min_len", help="Min total lenght", type=int, default=100)
parser.add_argument("--align_min_len", help="TIR minimun aligmnent length", type=int, default=10)
parser.add_argument("--tsd_min_len", help="TSD min lenght", type=int, default=2)
parser.add_argument("--tsd_max_len", help="TSD max lenght", type=int, default=10)
parser.add_argument("--FSL", help="Flanking seq length for comparison", type=int, default=50)
parser.add_argument("--min_copy_number", help="Minimum CN for families", type=int, default=3)

parser.add_argument("--only_cluster", help="Cluster candidates step",action='store_true')

args = parser.parse_args()#pylint: disable=invalid-name

#write results
if not os.path.isdir("results/" + args.jobname):
    os.mkdir("results/" + args.jobname)
if not os.path.isdir("results/" + args.jobname + "/temp/"):
    os.mkdir("results/" + args.jobname + "/temp/")

file_candidates_fasta = "results/" + args.jobname + "/candidates.fasta"
file_candidates_cluster = "results/" + args.jobname + "/candidates.fasta.cluster"
file_candidates_csv = "results/" + args.jobname + "/candidates.csv"
file_candidates_partial_prefix = "results/" + args.jobname + "/temp/candidates.partial."
candidates_partial_cluster = "results/" + args.jobname + "/candidates.cluster.fasta"
candidates_partial_repr_prefix =  "results/" + args.jobname + "/temp/candidates.representative.partial."
all_file = "results/" + args.jobname + "/all.fasta"
families_file = "results/" + args.jobname + "/families.fasta"

filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "results/" + args.jobname + "/out.log")
logging.basicConfig(
    filename=filename,
    level=logging.DEBUG,
    handlers=[
    logging.FileHandler(filename),
    logging.StreamHandler()],
    format="%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")

MITE_MAX_LEN = args.mite_max_len
align_min_len = args.align_min_len
MIN_TSD_LEN = args.tsd_min_len
MAX_TSD_LEN = args.tsd_max_len
MITE_MIN_LEN = args.mite_min_len

start_time = time.time()

def makelog(stri, do_print=True):
    if do_print:
        print(stri)
    logging.debug(stri)

def cur_time():
    global start_time
    elapsed_time = time.time() - start_time
    return "%f secs" % (elapsed_time,)

if not args.only_cluster:
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
    l_lock = Lock()
    perc_seq = {}
    last_perc_seq = {}
    irs = {}

    for i in range(args.workers):
        worker = Thread(target=findir.findIR, args=(q,args,l_lock, irs, perc_seq, last_perc_seq))
        worker.setDaemon(True)
        worker.start()
    windows_size = MITE_MAX_LEN * 2
    #windows_size = 5000 #int(ceil(MITE_MAX_LEN * 30))

    #processes until certain amount of sequences
    #stablish a balance between memory usage and processing
    max_queue_size = 50
    current_processing_size = 0
    #initialize global variables
    #start adding sequences to process queue
    record_count = 0
    fasta_seq = SeqIO.parse(args.genome, 'fasta')
    for record in fasta_seq:
        processed = False
        total_queue_count = 1
        queue_count = 1
        record_count += 1
        porc_ant = 0
        split_index = MAX_TSD_LEN + args.FSL
        clean_seq = ''.join(str(record.seq).splitlines())
        seq_len = len(clean_seq)
        params = (record.id, seq_len, record_count , seqs_count, (record_count * 100 / seqs_count), cur_time() )
        makelog("Adding %s (len %i) %i/%i (%i%% of total sequences in %s)" % params)
        while split_index < seq_len - MAX_TSD_LEN - args.FSL:
            seq = clean_seq[split_index:split_index + windows_size]
            seq_fs = clean_seq[split_index - args.FSL :split_index + windows_size + args.FSL]
            q.put((seq, seq_fs, split_index,record.id,seq_len,))
            queue_count += 1
            total_queue_count += 1
            split_index += MITE_MAX_LEN
            current_processing_size += MITE_MAX_LEN
            if q.qsize() >= max_queue_size:
                current_processing_size = 0
                q.join()
                processed = True
            #in order to avoid overloading of memory, we add a join()
            #we do not direcly use the join() method so we can process several
            #small sequences at once

    #In case of unprocessed sequences are left, let's wait
    q.join()

    makelog("Search for nested elements")

    labels = ['start','end','seq','record','len','ir_1','ir_2','tsd','tsd_in','fs_left','fs_right', 'ir_length','candidate_id','status','cluster']
    df = pd.DataFrame.from_records(irs.values(), columns=labels)
    makelog("Initial candidates: " + str(len(df)))
    makelog(cur_time())
    #filter out nested (keep larger)
    l=[]
    for idx, row in df.iterrows():
        #filter all that are nested into this
        res = df[(df.record == row.record) & (df.start + row.ir_length >= row.start) & (df.end  - row.ir_length <= row.end) & (df.index != idx) ]
        df.drop(res.index,inplace=True)
        #res = df[(df.record == row.record) & (df.start >= row.start) & (df.end <= row.end ) & (df.index != idx) ]
        l.append(res)
    res = pd.concat(l)
    makelog("Valid candidates (not nested): " + str(len(df)))
    makelog(cur_time())
    count = 1
  
    """
    fs_seqs = []
    irs_seqs = []
    df = df.sort_values(by=['record','start','end'])
    count = 1
    #positions = {}
    for index, row in df.iterrows():
        name = 'MITE_CAND_' + str(count)
        #positions[name] = (row.start, row.end)
        df.loc[index, 'candidate_id'] = name
        #append sequence record for biopython
        params = (row.record, row.start, row.end, row.tsd, row.tsd_in, row.len, row.ir_1, row.ir_2)
        description = "SEQ:%s START:%i END:%i TSD:%s TSD_IN:%s MITE_LEN:%i IR_1:%s IR_2:%s " % (params)
        ir_seq_rec = SeqRecord(Seq(row.seq), id=name, description = description)
        irs_seqs.append(ir_seq_rec)
        count += 1
    makelog("Writing candidates sequences")
    #SeqIO.write(irs_seqs, file_candidates_fasta , "fasta")
    df.to_csv(file_candidates_csv, index=False)
    makelog(cur_time())"""
#    df.to_csv(file_candidates_csv, index=False)
    count = 0
    positions = {}
    df = df.sort_values('start')
    for index, row in df.iterrows():
        name = 'MITE_CAND_' + str(count)
        df.loc[index, 'candidate_id'] = name
        positions[name] = (row.start, row.end)
        count += 1

    df.to_csv(file_candidates_csv, index=False)

if args.only_cluster:
    df = pd.read_csv(file_candidates_csv)
    count = 0
    positions = {}
    for index, row in df.sort_values('start').iterrows():
        name = 'MITE_CAND_' + str(count)
        df.loc[index, 'candidate_id'] = name
        positions[name] = (row.start, row.end)
        count += 1

max_len = int(df[['len']].max())
min_len = int(df[['len']].min())
num_files = 10
sep_size = ((max_len - min_len) / num_files) 
margin = sep_size * 0.25
last = min_len
count = 0
filtered_clusters = {}
irs_seqs_total = []
print min_len, max_len, sep_size, margin
valid_seqs = []
named_seqs = {}
for i in range(1, num_files + 1):
    makelog("Creating file for clustering " + str(i))
    curr = (i * sep_size) + min_len # separations
    print last, curr
    current_seqs = df[(df.len >= last) & (df.len <= curr)]
    irs_seqs = []
    for index, row in current_seqs.iterrows():
        params = (row.record, row.start, row.end, row.tsd, row.tsd_in, row.len)
        description = "SEQ:%s START:%i END:%i TSD:%s TSD_IN:%s MITE_LEN:%i" % (params)
        named_seqs[row.candidate_id] = (description, row.seq)
        ir_seq_rec = SeqRecord(Seq(row.seq), id=row.candidate_id, description = description)
        irs_seqs.append(ir_seq_rec)
        #irs_seqs_total.append(ir_seq_rec)
    file_candidates_partial = file_candidates_partial_prefix + str(i) + ".fasta"
    SeqIO.write(irs_seqs, file_candidates_partial , "fasta")
    makelog("Start clustering file  " + str(i))
    candidates_partial_repr = candidates_partial_repr_prefix + str(i) + ".fasta"
    cmd_list = [
    './cdhit/cd-hit-est',
    '-i',file_candidates_partial,
    '-o',candidates_partial_repr,
    '-c', '0.80','-n','7','-d','0','-T','0','-aL','0.8','-s','0.8','-M','0']
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    for c in iter(lambda: p.stdout.read(), ''):
        makelog(c)
    makelog("Done clustering file  " + str(i))
    clusters_dic = clusterutils.loadcluster(candidates_partial_repr + ".clstr")
    new_clusters = clusterutils.filtercluster(clusters_dic, args.min_copy_number,positions)
    #import ipdb; ipdb.set_trace()
    valid_seqs += [item for sublist in new_clusters.values() for item in sublist]
#    filtered_clusters = dict(filtered_clusters.items() + new_clusters.items())
    makelog("Done processing file  " + str(i))

    last = (i * sep_size) - margin + min_len

#write only valid sequences
valid_seq_records = []
for k,v in named_seqs.items():
    if k in valid_seqs:
        description, seq = v
        ir_seq_rec = SeqRecord(Seq(seq), id=k, description=description)
        valid_seq_records.append(ir_seq_rec)

SeqIO.write(valid_seq_records, file_candidates_fasta , "fasta")

import cdhitcluster
cdhitcluster.cluster(file_candidates_fasta, file_candidates_cluster, positions, min_copy_number, df)

makelog("Clustering valid sequences")
cmd_list = [
'./cdhit/cd-hit-est',
'-i',file_candidates_fasta,
'-o',file_candidates_cluster,
'-c', '0.80','-n','7','-d','0','-T','0','-aL','0.8','-s','0.8','-M','0']
p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
for c in iter(lambda: p.stdout.read(), ''):
    makelog(c)
#out,err = p.communicate()
makelog("Clustering done")

clusters_dic = clusterutils.loadcluster(file_candidates_cluster + ".clstr")
filtered_clusters = clusterutils.filtercluster(clusters_dic, args.min_copy_number, positions)
unique_clusters = set(filtered_clusters.keys())
num_clusters = len(unique_clusters)

#loop through clusters
for current_cluster in unique_clusters:
    #search candidates for that cluster
    #all possible 2-combinations of candidates
    candidates = filtered_clusters[current_cluster]
    combinations = [(x,y) for x,y in itertools.combinations(candidates, 2)]
    dist_fs = {}
    for seq_id in combinations:
        x,y = seq_id

        if not x in candidates or not y in candidates:
            continue
        cand_x = df[(df.candidate_id == x)]
        cand_y = df[(df.candidate_id == y)]

        #if they're partially overlapped, ignore flanking sequence comparison
        if cand_x.iloc[0].end >= cand_y.iloc[0].start and cand_y.iloc[0].end >= cand_x.iloc[0].start:
            continue

        fs_right_1 = cand_x.iloc[0].fs_right
        fs_left_1 = cand_x.iloc[0].fs_left

        fs_right_2 = cand_y.iloc[0].fs_right
        fs_left_2 = cand_y.iloc[0].fs_left
        
        score_r1_r2 = pairwise2.align.localms(fs_right_1, fs_right_2, 1, -1, -1, -1,score_only=True)
        score_l1_l2 = pairwise2.align.localms(fs_left_1, fs_left_2, 1, -1, -1, -1,score_only=True)
        score_r1_l2 = pairwise2.align.localms(fs_right_1, fs_left_2, 1, -1, -1, -1,score_only=True)
        score_r2_l1 = pairwise2.align.localms(fs_right_2, fs_left_1, 1, -1, -1, -1,score_only=True)
        max_score = max(score_r1_r2,score_l1_l2,score_r1_l2,score_r2_l1)

        if max_score == []:
            max_score = 0
        max_score /= args.FSL
        #todo validate scoring
        if max_score < 0.5:
            dist_fs[x] = 1
            dist_fs[y] = 1
    if len(dist_fs) < args.min_copy_number:
        df.loc[df['candidate_id'].isin(filtered_clusters[current_cluster]), 'status'] =  'low_cn_flank_seq'
        del filtered_clusters[current_cluster]

#again to remove < MIN_COPY_NUMBER elements
#filtered_clusters = clusterutils.filtercluster(filtered_clusters, args.min_copy_number, positions, df, 'low_copy_number_2')
ordered_cluster = OrderedDict(sorted(filtered_clusters.items(), key=lambda t: t[1]))

makelog("Clusters: " + str(len(filtered_clusters)) + " writing sequences")

fasta_seq = SeqIO.parse(file_candidates_fasta, 'fasta')
buffer_rec = []
for record in fasta_seq:
    for clus, seqs in filtered_clusters.items():
        if record.id in seqs:
            df.loc[df['candidate_id'] == record.id, 'cluster'] =  clus
            df.loc[df['candidate_id'] == record.id, 'status'] =  'valid'
            record.description = "%s CLUSTER:%s" % (record.description, clus)
            buffer_rec.append(record)
            continue

SeqIO.write(buffer_rec, all_file , "fasta")

df.to_csv(file_candidates_csv, index=False)
clusterutils.cluster2seq(ordered_cluster, file_candidates_fasta, families_file )
makelog(cur_time())