import pandas as pd
from Bio import SeqIO
import os


file_repbase_total = "rebpase_detect_total.csv"
file_repbase_filtered = "repbase_detect_filtered.csv"

df_repbase_total = pd.read_csv(file_repbase_total, delimiter="\t")
df_repbase_total.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

df_repbase_filtered = pd.read_csv(file_repbase_filtered, delimiter="\t")
df_repbase_filtered.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

df_res = df_repbase_total[~df_repbase_total.sseqid.isin(df_repbase_filtered.sseqid)]


total = len(pd.unique(df_repbase_total.qseqid))
mites = len(pd.unique(df_repbase_filtered.qseqid))
diff = len(pd.unique(df_res.qseqid))
print 'Total: %s' % (total, )
print 'MITEs: %s' % (mites, )
print 'Diff: %s' % (diff, )
print  (diff * 100 / total)
#print df_res.qseqid


file_repbase_total = "rebpase_total.csv"
file_repbase_filtered = "repbase_filtered.csv"

df_repbase_total = pd.read_csv(file_repbase_total, delimiter="\t")
df_repbase_total.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

df_repbase_filtered = pd.read_csv(file_repbase_filtered, delimiter="\t")
df_repbase_filtered.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

df_res = df_repbase_total[~df_repbase_total.sseqid.isin(df_repbase_filtered.sseqid)]

total = len(pd.unique(df_repbase_total.qseqid))
mites = len(pd.unique(df_repbase_filtered.qseqid))
diff = len(pd.unique(df_res.qseqid))
print 'Total: %s' % (total, )
print 'MITEs: %s' % (mites, )
print 'Diff: %s' % (diff, )
print  (diff * 100 / total)