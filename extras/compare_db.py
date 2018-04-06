#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import argparse
from collections import OrderedDict
import pandas as pd
import os

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-m", "--mites", help="All elements from program to analyze", required=True)
parser.add_argument("-d", "--db", help="MITE elements from database", required=True)
parser.add_argument("-p", "--program", help="Program name", required=True)
parser.add_argument("-n", "--dbname", help="Database name", required=True)
parser.add_argument("-l", "--label", help="label name", required=True)
args = parser.parse_args()#pylint: disable=invalid-name


#extract elements ids from program
fasta_index = SeqIO.index(args.mites, 'fasta')
fasta_index = OrderedDict(sorted(fasta_index.items()))
elem = []
for record_1 in fasta_index:
    elem.append(record_1)

df_all_mites = pd.DataFrame({'qseqid': elem})

#extract elements ids from database
fasta_index = SeqIO.index(args.db, 'fasta')
fasta_index = OrderedDict(sorted(fasta_index.items()))
elem_repbase = []
for record_1 in fasta_index:
    elem_repbase.append(record_1)
df_db_mites = pd.DataFrame({'sseqid': elem_repbase})

#extract elements from program matching some element in database
cmd = 'blastn -query %s  -subject %s -outfmt 6  > %s'
cmd = cmd % (args.mites, args.db, args.program + '_db')
print cmd
os.system(cmd)
df_blast_p_db = pd.read_csv(args.program + '_db', delimiter="\t", header=None)
df_blast_p_db.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]
#query = program
#subject = database

#db covered elements
values_covered = pd.unique(df_blast_p_db.sseqid)
df_covered = pd.DataFrame({'sseqid': values_covered})
#all results from database elements

#db uncovered elements
df_uncovered = df_all_mites[~df_all_mites.qseqid.isin(df_blast_p_db.qseqid)]
values = pd.unique(df_uncovered.qseqid)
df_uncovered = pd.DataFrame({'qseqid': values})
#all elements from program that are not in results

#db uncovered elements sseqname
df_uncovered_db = df_db_mites[~df_db_mites.sseqid.isin(df_blast_p_db.sseqid)]
values = pd.unique(df_uncovered_db.sseqid)
df_uncovered_db = pd.DataFrame({'qseqid': values})
#all elements from database that are not in results

elems = (df_uncovered.qseqid.tolist() + df_covered.sseqid.tolist())
df_program = pd.DataFrame({'qseqid': elems})

set1 = set(df_db_mites.sseqid)
set2 = set(df_program.qseqid)

v = venn2([set1, set2], (args.dbname, args.program))

for text in v.set_labels:
    text.set_fontsize(17)
for text in v.subset_labels:
    text.set_fontsize(15)


db_elements_len = len(values_covered)
pm_elements_len = len(pd.unique(df_blast_p_db.qseqid))
perc_coverage = db_elements_len * 100 / len(elem_repbase)
v.get_label_by_id('11').set_text('%i|%i\n(%i%%)' % (db_elements_len,pm_elements_len,perc_coverage ))
plt.title('%s)' % (args.label,))
outname = args.dbname + "_" + args.program + ".png"
print outname
plt.savefig(outname, dpi=600)
