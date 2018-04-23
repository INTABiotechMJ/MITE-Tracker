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
parser.add_argument("-e1", "--elements_1", help="All elements from program to analyze")
parser.add_argument("-e2", "--elements_2", help="MITE elements from database")
parser.add_argument("-e1n", "--elements_1_name", help="Program name")
parser.add_argument("-e2n", "--elements_2_name", help="Database name")
parser.add_argument("-l", "--label", help="display label name")
parser.add_argument("-o", "--output", help="New elements fasta")
args = parser.parse_args()#pylint: disable=invalid-name

def compare(elements_1, elements_2, elements_1_name, elements_2_name, label,ax,savefigname = False):
    #extract elements ids from program
    fasta_index = SeqIO.index(elements_1, 'fasta')
    fasta_index = OrderedDict(sorted(fasta_index.items()))
    elem = []
    for record_1 in fasta_index:
        elem.append(record_1)

    df_program_mites = pd.DataFrame({'qseqid': elem})

    #extract elements ids from database
    fasta_index = SeqIO.index(elements_2, 'fasta')
    fasta_index = OrderedDict(sorted(fasta_index.items()))
    elem_repbase = []
    for record_1 in fasta_index:
        elem_repbase.append(record_1)
    df_db_mites = pd.DataFrame({'sseqid': elem_repbase})

    #extract elements from program matching some element in database
    #cmd = 'blastn -word_size 12 -qcov_hsp_perc 80  -query %s  -subject %s -outfmt 6  > %s'
    cmd = 'blastn -task blastn -evalue 10e-3 -qcov_hsp_perc 80 -query %s  -subject %s -outfmt 6  > %s'
    #cmd = 'blastn -query %s  -subject %s -outfmt 6  > %s'
    dbname = elements_2_name.replace(' ','_') + '_' + elements_1_name.replace(' ','_') + '_db'
    cmd = cmd % (elements_1.replace(' ','_'), elements_2.replace(' ','_'),  dbname)
    print(cmd)
    os.system(cmd)
    df_blast_p_db = pd.read_csv(dbname, delimiter="\t", header=None)
    #os.remove(dbname)
    df_blast_p_db.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    #again backwards
    cmd = 'blastn -task blastn -evalue 10e-3 -qcov_hsp_perc 80 -subject %s  -query %s -outfmt 6  > %s'
    #cmd = 'blastn -query %s  -subject %s -outfmt 6  > %s'
    dbname = elements_2_name.replace(' ','_') + '_' + elements_1_name.replace(' ','_') + '_db'
    cmd = cmd % (elements_1.replace(' ','_'), elements_2.replace(' ','_'),  dbname)
    print(cmd)
    os.system(cmd)
    df_blast_p_db_back = pd.read_csv(dbname, delimiter="\t", header=None)
    #os.remove(dbname)
    df_blast_p_db_back.columns = ['sseqid','qseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]
    frames = [df_blast_p_db, df_blast_p_db_back]
    df_blast_p_db = pd.concat(frames)

    #query = program
    #subject = database

    #db covered elements
    values_covered = pd.unique(df_blast_p_db.sseqid)
    df_covered = pd.DataFrame({'sseqid': values_covered})
    #all results from database elements

    #db uncovered elements
    df_uncovered = df_program_mites[~df_program_mites.qseqid.isin(df_blast_p_db.qseqid)]
    values = pd.unique(df_uncovered.qseqid)
    df_uncovered = pd.DataFrame({'qseqid': values})
    #all elements from program that are not in results

    #db uncovered elements sseqname
    df_uncovered_db = df_db_mites[~df_db_mites.sseqid.isin(df_blast_p_db.sseqid)]
    values_uncovered_db = pd.unique(df_uncovered_db.sseqid)
    df_uncovered_db = pd.DataFrame({'qseqid': values_uncovered_db})
    #all elements from database that are not in results

    elems = (df_uncovered.qseqid.tolist() + df_covered.sseqid.tolist())
    df_program = pd.DataFrame({'qseqid': elems})

    set1 = set(df_db_mites.sseqid)
    set2 = set(df_program.qseqid)

    #new elements found
    df_new_elements = df_program_mites[~df_program_mites.qseqid.isin(df_blast_p_db.qseqid)]

    v = venn2([set1, set2], (elements_2_name, elements_1_name),ax=ax)

    for text in v.set_labels:
        text.set_fontsize(10)
    for text in v.subset_labels:
        text.set_fontsize(8)


    db_elements_len = len(values_covered)
    pm_elements_len = len(pd.unique(df_blast_p_db.qseqid))

    new_elements_list = pd.unique(df_new_elements.qseqid)
    new_elements_len = len(new_elements_list)
    perc_coverage = db_elements_len * 100 / len(elem_repbase)
    v.get_label_by_id('11').set_text('%i|%i\n(%i%%)' % (db_elements_len,pm_elements_len,perc_coverage ))
    if label:
        title = '%s)' % (label,)
        ax.set_title(title)

    outname = elements_2_name + "_" + elements_1_name + ".png"
    print(outname)
    print("Program matches elements: %i" % (pm_elements_len,))
    print("DB matches elements: %i" % (db_elements_len,))
    print("New elements: %i" % (new_elements_len,))
    print("Uncovered db elements: ")
    print(values_uncovered_db)

    if savefigname:
        plt.savefig(savefigname, dpi=800,bbox_inches='tight')


    fasta_file = SeqIO.parse(elements_1, 'fasta')

    new_records = []
    for record in fasta_file:
        if record.id in new_elements_list:
            new_records.append(record)
            #print record
    SeqIO.write(new_records, "newseqs" + elements_1_name + " " + elements_2_name, 'fasta')

if __name__ == "__main__":
    #compare(args.elements_1, args.elements_2, args.elements_1_name, args.elements_2_name, args.label, args.output)
    figure, axes = plt.subplots(2, 1)

    compare("data/tracker/families_nr.fasta","data/mitehunter/all.fa","MITE Tracker" ,"MITE Hunter", "A",axes[0])
    compare("data/tracker/families_nr.fasta","data/detectmite/rice.mite.fasta","MITE Tracker" ,"detectMITE", "B",axes[1],"a.png")
    
    plt.clf()
    figure, axes = plt.subplots(3, 1)

    compare("data/tracker/families_nr.fasta","data/repbase/repbase_nonautonomous.fasta","MITE Tracker" ,"Repbase", "A",axes[0])
    compare("data/detectmite/rice.mite.fasta","data/repbase/repbase_nonautonomous.fasta","detectMITE" ,"Repbase", "B",axes[1])
    compare("data/mitehunter/all.fa","data/repbase/repbase_nonautonomous.fasta","MITE Hunter" ,"Repbase", "C",axes[2],"b.png")


    plt.clf()
    figure, axes = plt.subplots(1, 1)

    compare("data/tracker/families_nr.fasta","data/all/all.fasta","MITE Tracker" ,"All", "A", axes)
    