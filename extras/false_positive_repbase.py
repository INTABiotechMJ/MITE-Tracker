import pandas as pd
from Bio import SeqIO
import os


eles = []
eles.append( ('detectMITE', 'data/detectmite/rice.mite.fasta') )
eles.append( ('MITE_Tracker', 'data/tracker/families_nr.fasta') )
eles.append( ('MITE_Hunter', 'data/mitehunter/all.fa') )
repbase = 'data/repbase/repbase_oryzativa.fasta'
repbase_filtered = 'data/repbase/repbase_nonautonomous.fasta'

for ele in eles:
    name, families = ele
    name_repbase = name + '_repbase'
    name_repbase_filtered = name + '_repbase_filtered'
    cmd = 'blastn -task blastn -evalue 10e-3 -qcov_hsp_perc 80  -query %s  -subject %s -outfmt 6  > %s'
    cmd = cmd % (families, repbase, name_repbase)
    print cmd
    os.system(cmd)


    cmd = 'blastn -task blastn -evalue 10e-3 -qcov_hsp_perc 80  -query %s  -subject %s -outfmt 6  > %s'
    cmd = cmd % (families, repbase_filtered, name_repbase_filtered)
    print cmd
    os.system(cmd)

    df_repbase_total = pd.read_csv(name_repbase, delimiter="\t")
    df_repbase_total.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    df_repbase_filtered = pd.read_csv(name_repbase_filtered, delimiter="\t")
    df_repbase_filtered.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    df_res = df_repbase_total[~df_repbase_total.qseqid.isin(df_repbase_filtered.qseqid)]

    total_query = len(pd.unique(df_repbase_total.qseqid))
    total_subject = len(pd.unique(df_repbase_total.sseqid))
    mites_covered = len(pd.unique(df_repbase_filtered.sseqid))
    diff = len(pd.unique(df_res.qseqid))

    print '*' * 10
    print name
    print 'Total query: %s' % (total_query, )
    print 'Total subject: %s' % (total_subject, )
    print 'MITEs: %s' % (mites_covered, )
    print 'Total Repbase NOT MITEs: %s' % (diff, )
    print 'Wrong mites', pd.unique(df_res.sseqid)
    print  str((diff * 100 / total_query)) + "%"