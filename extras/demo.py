import pandas as pd

#TREP Comparison

for i in ("trep_detect.csv", "trep_tracker.csv"):
    print "-----File: ",i
    #blastn -query ../rice_wz_6_mcn_3/families_nr.fasta   -max_target_seqs 1     -subject ../../../data/otrep.fa  -outfmt 6 > trep_tracker.csv
    #blastn -query ../detectmite/rice.mite.fasta     -max_target_seqs 1     -subject ../../../data/otrep.fa  -outfmt 6 > trep_detect.csv          
    df = pd.read_csv(i, delimiter="\t")
    df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]
    rm, fm = 0, {}
    for k,row in df.iterrows():
        mite = row.sseqid
        #print mite
        if mite.startswith('XXX'):
            continue
        if mite.startswith('DT'):
            rm += 1
        else:
            print "False positive: ",row.qseqid, '-> '  ,mite
            fm[row.qseqid] = 1
    tm = rm + len(fm)
    print "Total elements: ", tm
    print "False positives: ", len(fm)
    print "% false positive:",100 - (rm * 100.0 / tm)