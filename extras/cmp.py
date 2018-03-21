import pandas as pd
from Bio import SeqIO
import os



file_repbase = "../data/repbase_nonautonomous_rice.fa"
file_repbase_filtered = "../data/repbase_nonautonomous_rice.50_800.fa"
file_blast_repbase_target = "blast_repbase_target.fa"

fasta_seq = SeqIO.parse(file_repbase, 'fasta')
seqs = []
seq_repbase = []
for record in fasta_seq:
    clean_seq = ''.join(str(record.seq).splitlines())
    seq_len = len(clean_seq)    
    #filter repbase 50 >= x >= 800
    if seq_len < 50 or seq_len > 800:
        continue
    seq_repbase.append(record.id)
    seqs.append(record)
SeqIO.write(seqs, file_repbase_filtered , 'fasta')

targets = []
MITETracker = ('MITE Tracker',"results/rice_wz_6/all.fasta")
#MITETracker = ('detectMITE',"../rice_new/all.fasta")
targets.append(MITETracker)

print 'Repbase total seqs %s' % (str(len(seqs)))
for target in targets:
    name, file_target = target
    #get all elements in target
    seq_target = []
    fasta_seq = SeqIO.parse(file_target, 'fasta')
    for record in fasta_seq:
        seq_target.append(record.id)

    #blast repbase and target
    os.system('blastn -query %s -subject %s -outfmt 6 > %s' % (file_target,file_repbase_filtered,file_blast_repbase_target))

    df_target = pd.DataFrame.from_records({'id':seq_target},columns=['id'])
    df_repbase = pd.DataFrame.from_records({'id':seq_repbase},columns=['id'])

    df_blast = pd.read_csv(file_blast_repbase_target, delimiter="\t")
    df_blast.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

    print '%s total seqs %s' % (name, str(len(seq_target)))
    print 'repbase total seqs %s' % (str(len(df_repbase)))
    print '%s has %s hits into repbase' % (name, str(len(df_blast.qseqid.unique())))
    print 'repbase has %s hits from %s' % (str(len(df_blast.sseqid.unique())), name)
    import ipdb; ipdb.set_trace()
    #print '%s has %i sequences that are not in repbase' % (name, str(len(df_blast)))
    #print 'Repbase has %i sequences that are not in %s' % (str(len(seq_target)), seq_target)