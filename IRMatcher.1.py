# -*- coding: utf-8 -*-
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
import cdhitutils
from collections import OrderedDict
import sys
import pandas as pd
from threading import Thread, Lock, active_count
import Queue
import logging
import findir

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
parser.add_argument("--mite_max_len", help="MITE max lenght", type=int, default=800)
parser.add_argument("--mite_min_len", help="Min total lenght", type=int, default=100)
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

MITE_MAX_LEN = args.mite_max_len
align_min_len = args.align_min_len
MIN_TSD_LEN = args.tsd_min_len
MAX_TSD_LEN = args.tsd_max_len
MITE_MIN_LEN = args.mite_min_len
FSL = args.FSL
MIN_COPY_NUMBER = args.min_copy_number


cluster_candidates_file = "results/" + args.jobname + "/candidates.representative.fasta"

clusters_dic = cdhitutils.loadcluster(cluster_candidates_file + ".clstr")
filtered_clusters = cdhitutils.filtercluster(clusters_dic, args.min_copy_number,positions)
unique_clusters = set(filtered_clusters.keys())
num_clusters = len(unique_clusters)
#loop through clusters
for current_cluster in unique_clusters:
    #search candidates for that cluster
    #all possible 2-combinations of candidates
    candidates = filtered_clusters[current_cluster]
    combinations = [(x,y) for x,y in itertools.combinations(candidates, 2)]
    dist_fs = 0
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
        max_score /= FSL
        #todo validate scoring
        if max_score < 0.5:
            dist_fs += 1
    if dist_fs < args.min_copy_number:
        makelog(str(max_score) + x + y)
        print filtered_clusters[current_cluster]
        print current_cluster
        del filtered_clusters[current_cluster]#.remove(x)
        #filtered_clusters[current_cluster].remove(y)

#again to remove < MIN_COPY_NUMBER elements
filtered_clusters = cdhitutils.filtercluster(filtered_clusters, args.min_copy_number, positions)
ordered_cluster = OrderedDict(sorted(filtered_clusters.items(), key=lambda t: t[1]))

makelog("Clusters: " + str(len(filtered_clusters)))

fasta_seq = SeqIO.parse(candidates_fasta, 'fasta')
buffer_rec = []
for record in fasta_seq:
    for seqs in filtered_clusters.values():
        if record.id in seqs:
            buffer_rec.append(record)
            continue
all_file = "results/" + args.jobname + "/all.fasta"
SeqIO.write(buffer_rec, all_file , "fasta")

cdhitutils.cluster2seq(ordered_cluster, candidates_fasta, "results/" + args.jobname + "/families.fasta" )
makelog(cur_time())