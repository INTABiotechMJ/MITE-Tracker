# -*- coding: utf-8 -*-
from __future__ import division
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
from Bio.Seq import Seq
from subprocess import Popen, PIPE
import os, shutil
import time
import argparse
from math import ceil
from collections import OrderedDict
import sys
import pandas as pd
from threading import Thread, Lock, active_count
import queue
import logging
import findir

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genome", help="Genome file in fasta format", required=True)
parser.add_argument("-j","--jobname", help="Will create files under a folder called [jobname]", required=True)
parser.add_argument("-w","--workers", help="Max number of processes to use simultaneously", type=int, default=1)
parser.add_argument("--mite_max_len", help="MITE max lenght", type=int, default=800)
parser.add_argument("--mite_min_len", help="Min total lenght", type=int, default=50)
parser.add_argument("--align_min_len", help="TIR minimun aligmnent length", type=int, default=10)
parser.add_argument("--tsd_min_len", help="TSD min lenght", type=int, default=2)
parser.add_argument("--tsd_max_len", help="TSD max lenght", type=int, default=10)
parser.add_argument("--FSL", help="Flanking seq length for comparison", type=int, default=50)
parser.add_argument("--min_copy_number", help="Minimum CN for families", type=int, default=3)
parser.add_argument("--task", help="Task: all|candidates|cluster (default=all)", default='all')
args = parser.parse_args()#pylint: disable=invalid-name


#file names
file_names = {}
file_names['file_representative'] = "results/" + args.jobname + "/families_nr.fasta"
file_names['file_candidates_fasta'] = "results/" + args.jobname + "/candidates.fasta"
file_names['file_gff'] = "results/" + args.jobname + "/all.gff3"
file_names['file_candidates_cluster'] = "results/" + args.jobname + "/candidates.fasta.cluster"
file_names['file_candidates_csv'] = "results/" + args.jobname + "/candidates.csv"
file_names['file_temp_cluster_dir'] = "results/" + args.jobname + "/temp/"
file_names['file_temp_cluster'] = file_names['file_temp_cluster_dir'] + "clust"
file_names['file_candidates_partial_prefix'] = "results/" + args.jobname + "/temp/candidates.partial."
file_names['candidates_partial_cluster'] = "results/" + args.jobname + "/candidates.cluster.fasta"
file_names['candidates_partial_repr_prefix'] =  "results/" + args.jobname + "/temp/candidates.representative.partial."
file_names['all_file'] = "results/" + args.jobname + "/all.fasta"
file_names['families_file'] = "results/" + args.jobname + "/families.fasta"
file_names['file_candidates_dir'] = "results/" + args.jobname + "/temp/"

#write dir for results
if not os.path.isdir("results/" + args.jobname):
    os.mkdir("results/" + args.jobname)

if not os.path.isdir("results/" + args.jobname + "/tmp"):
    os.mkdir("results/" + args.jobname + "/tmp")

if os.path.isdir(file_names['file_temp_cluster_dir']):
    shutil.rmtree(file_names['file_temp_cluster_dir'])
os.mkdir(file_names['file_temp_cluster_dir'])


filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "results/" + args.jobname + "/out.log")
logging.basicConfig(filename=filename,level=logging.DEBUG,format='%(asctime)s %(message)s')

MITE_MAX_LEN = args.mite_max_len
MIN_TSD_LEN = args.tsd_min_len
MAX_TSD_LEN = args.tsd_max_len
MITE_MIN_LEN = args.mite_min_len
start_time = time.time()

def makelog(stri, do_print=True):
    if do_print:
        print(stri)
    logging.debug(stri)

if not args.task in ('all', 'cluster', 'candidates'):
    makelog('task parameter not valid')
    exit()

def cur_time():
    global start_time
    elapsed_time = time.time() - start_time
    return "%f secs" % (elapsed_time,)

if args.task == 'all' or args.task == 'candidates':
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
    q = queue.Queue(maxsize=0)
    l_lock = Lock()
    perc_seq = {}
    last_perc_seq = {}
    candidates = {}
    workers = []
    for i in range(args.workers):
        worker = Thread(target=findir.findIR, args=(q,args,l_lock, candidates, perc_seq, last_perc_seq))
        worker.setDaemon(True)
        workers.append(worker)
        #worker.start()
    windows_size = MITE_MAX_LEN * 2
    #windows_size = 5000 #int(ceil(MITE_MAX_LEN * 30))

    #processes until certain amount of sequences
    #stablish a balance between memory usage and processing
    max_queue_size = 50
    #initialize global variables
    #start adding sequences to process queue
    record_count = 0
    fasta_seq = SeqIO.parse(args.genome, 'fasta')
    worker_started = False
    for record in fasta_seq:
        #it takes too long to reach here and workers may timeout if we started them early
        if not worker_started:
            worker_started = True
            for worker in workers:
                worker.start()
        processed = False
        queue_count = 1
        record_count += 1
        porc_ant = 0
        split_index = MAX_TSD_LEN + args.FSL
        clean_seq = ''.join(str(record.seq).splitlines())
        seq_len = len(clean_seq)
        params = (record.id, record_count , seqs_count, (record_count * 100 / seqs_count), cur_time() )
        makelog("Adding %s %i/%i (%i%% of total sequences in %s)" % params)
        while split_index < seq_len - MAX_TSD_LEN - args.FSL:
            seq = clean_seq[split_index:split_index + windows_size]
            seq_fs = clean_seq[split_index - args.FSL :split_index + windows_size + args.FSL]
            q.put((seq, seq_fs, split_index,record.id,seq_len,queue_count,))
            queue_count += 1
            split_index += MITE_MAX_LEN
            if q.qsize() >= max_queue_size:
                q.join()
                processed = True
    #In case of unprocessed sequences are left, let's wait
    q.join()
    total_candidates = {}
    count = 0
    irs_seqs = []
    for part in candidates.values():
        for candidate in part:
            #organize and name
            name = 'MITE_CAND_%i|%s|%s|%s|%s|%s' % (count,candidate['record'], candidate['start'], candidate['end'], candidate['tsd'],candidate['tir_len'],)
            count += 1
            candidate['candidate_id'] = name
            total_candidates[name] = candidate
            #record
            params = (candidate['tsd_in'], candidate['mite_len'],candidate['tir_len'])
            description = "TSD_IN:%s MITE_LEN:%i TIR_LEN:%i" % (params)
            candidate['description'] = description
            ir_seq_rec = SeqRecord(Seq(candidate['seq']), id=candidate['candidate_id'], description=description)
            irs_seqs.append(ir_seq_rec)
    makelog("Candidates: " + str(count))
    #write candidate fasta
    SeqIO.write(irs_seqs, file_names['file_candidates_fasta'] , "fasta")
    #write candidates csv
    to_csv = []
    for k,v in total_candidates.items():
        to_csv.append(v)
    df = pd.DataFrame(to_csv)
    df.to_csv(file_names['file_candidates_csv'], index=False)

if args.task == 'cluster':
    df = pd.read_csv(file_names['file_candidates_csv'])
    total_candidates = {}
    for index, row in df.iterrows():
        total_candidates[row.candidate_id] = row.to_dict()
        
if args.task == 'all' or args.task == 'cluster':
    import vsearchcluster
    vsearchcluster.cluster(file_names, total_candidates, args.min_copy_number, args.FSL, args.workers)
makelog(cur_time())