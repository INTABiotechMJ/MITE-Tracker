import pandas as pd
from Bio import SeqIO
import os


eles = []
eles.append( ('detectMITE', 'data/detect_repbase.csv','data/detect_repbase_filtered.csv') )
eles.append( ('MITE Tracker', 'data/tracker_repbase.csv','data/tracker_repbase_filtered.csv') )

for ele in eles:
    name, seq_total, seq_filtered = ele
    df_repbase_total = pd.read_csv(seq_total, delimiter="\t")
    df_repbase_total.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    df_repbase_filtered = pd.read_csv(seq_filtered, delimiter="\t")
    df_repbase_filtered.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    df_res = df_repbase_total[~df_repbase_total.qseqid.isin(df_repbase_filtered.qseqid)]

    total_query = len(pd.unique(df_repbase_total.qseqid))
    total_subject = len(pd.unique(df_repbase_total.sseqid))
    #mites = len(pd.unique(df_repbase_filtered.qseqid))
    diff = len(pd.unique(df_res.qseqid))

    print '*' * 10
    print name
    print 'Total query: %s' % (total_query, )
    print 'Total subject: %s' % (total_subject, )
    #print 'MITEs: %s' % (mites, )
    print 'Total Repbase NOT MITEs: %s' % (diff, )
    print 'Wrong mites', pd.unique(df_res.qseqid)
    print 'Wrong mites', pd.unique(df_res.sseqid)
    print  (diff * 100 / total_query)