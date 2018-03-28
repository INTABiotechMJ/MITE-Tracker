#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
#from matplotlib import pyplot as plt
#from matplotlib_venn import venn3, venn3_circles
import argparse
from collections import OrderedDict
import pandas as pd
import os

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-t", "--tracker", help="All elements from program to analyze", required=True)
parser.add_argument("-n", "--hunter", help="All elements from program to analyze", required=True)
parser.add_argument("-d", "--detect", help="All elements from program to analyze", required=True)
args = parser.parse_args()#pylint: disable=invalid-name


#extract elements ids from program
fasta_index = SeqIO.index(args.tracker, 'fasta')
fasta_index = sorted(fasta_index.items())
elem_tracker = []
for record_1 in fasta_index:
    elem_tracker.append(record_1)
df_tracker = pd.DataFrame({'qseqid': elem_tracker})

#extract elements ids from program
fasta_index = SeqIO.index(args.hunter, 'fasta')
fasta_index = sorted(fasta_index.items())
elem_hunter = []
for record_1 in fasta_index:
    elem_hunter.append(record_1)
df_hunter = pd.DataFrame({'qseqid': elem_hunter})

#extract elements ids from program
fasta_index = SeqIO.index(args.detect, 'fasta')
fasta_index = sorted(fasta_index.items())
elem_detect = []
for record_1 in fasta_index:
    elem_detect.append(record_1)
df_detect = pd.DataFrame({'qseqid': elem_detect})

#extract elements from program matching some element in database
cmd = 'blastn -query %s -max_target_seqs 1 -subject %s -outfmt "6 qseqid sseqid"  > t_d'
cmd = cmd % (args.tracker, args.detect)
os.system(cmd)
df_t_d = pd.read_csv('t_d', delimiter="\t", header=None)
df_t_d.columns = ['qseqid','sseqid']

#extract elements from program matching some element in database
cmd = 'blastn -query %s -max_target_seqs 1 -subject %s -outfmt "6 qseqid sseqid"  > t_h'
cmd = cmd % (args.tracker, args.hunter)
os.system(cmd)
df_t_h = pd.read_csv('t_h', delimiter="\t", header=None)
df_t_h.columns = ['qseqid','sseqid']

#extract elements from program matching some element in database
cmd = 'blastn -query %s -max_target_seqs 1 -subject %s -outfmt "6 qseqid sseqid"  > d_h'
cmd = cmd % (args.detect, args.hunter)
os.system(cmd)
df_d_h = pd.read_csv('d_h', delimiter="\t", header=None)
df_d_h.columns = ['qseqid','sseqid']

df_join_t_d_h = pd.merge(df_t_d, df_t_h, on='qseqid', how='outer')
df_join_t_d_h = pd.merge(df_join_t_d_h, df_d_h, on='qseqid', how='outer')
#df_join_t_d_h.to_csv("res.csv")

df_join_t_d_h = df_join_t_d_h[(~df_join_t_d_h.sseqid_x.isnull()) & (~df_join_t_d_h.sseqid_y.isnull())]
df_join_t_d_h.to_csv("res.csv")

print 'tracker:',len((elem_tracker))
print 'hunter:',len((elem_hunter))
print 'detect',len((elem_detect))
print 'TxD:',len(pd.unique(df_t_d.sseqid))
print 'TxH:',len(pd.unique(df_t_h.sseqid))
print 'DxH:',len(pd.unique(df_d_h.sseqid))
print 'TxDxH:',len(pd.unique(df_join_t_d_h.sseqid_y))

#print pd.unique(df_join_t_d_h.qseqid)
#set1 = set(df_db_mites.sseqid)
#set2 = set(df_program.qseqid)

#venn2([set1, set2], ('Repbase', args.program))
#plt.savefig(args.program + ".png", dpi=600)
