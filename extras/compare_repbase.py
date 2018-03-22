import pandas as pd
from Bio import SeqIO
import os


eles = []
eles.append( ('MITE Finder', '../../data/mite-tracker/finder_repbase.csv','../../data/mite-tracker/finder_repbase_filtered.csv') )
eles.append( ('detectMITE', '../../data/mite-tracker/detect_repbase.csv','../../data/mite-tracker/detect_repbase_filtered.csv') )
eles.append( ('MITE Tracker', '../../data/mite-tracker/tracker_repbase.csv','../../data/mite-tracker/tracker_repbase_filtered.csv') )

for ele in eles:
    name, seq_total, seq_filtered = ele
    df_repbase_total = pd.read_csv(seq_total, delimiter="\t")
    df_repbase_total.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    df_repbase_filtered = pd.read_csv(seq_filtered, delimiter="\t")
    df_repbase_filtered.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    df_res = df_repbase_total[~df_repbase_total.sseqid.isin(df_repbase_filtered.sseqid)]

    total = len(pd.unique(df_repbase_total.qseqid))
    mites = len(pd.unique(df_repbase_filtered.qseqid))
    diff = len(pd.unique(df_res.qseqid))

    print '*' * 10
    print name
    print 'Total: %s' % (total, )
    print 'MITEs: %s' % (mites, )
    print 'Total Repbase NOT MITEs: %s' % (diff, )
    print  (diff * 100 / total)