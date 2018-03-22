import pandas as pd

#TREP Comparison

for i in ("../../data/mite-tracker/tracker_trep.csv", "../../data/mite-tracker/detect_trep.csv", "../../data/mite-tracker/finder_trep.csv"):
    print "-----File: ",i
    #blastn -query ../rice_wz_6_mcn_3/families_nr.fasta   -max_target_seqs 1     -subject ../../../data/otrep.fa  -outfmt 6 > trep_tracker.csv
    #blastn -query ../detectmite/rice.mite.fasta     -max_target_seqs 1     -subject ../../../data/otrep.fa  -outfmt 6 > trep_detect.csv          
    df = pd.read_csv(i, delimiter="\t")
    df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]
    rm, fm = 0, {}
    for id_ in pd.unique(df.qseqid):
        mite = df[df.qseqid==id_].iloc[0].sseqid
        #print mite
        if mite.startswith('XXX'):
            continue
        if mite.startswith('DT'):
            rm += 1
        else:
            print "False positive: ",id_, '-> '  ,mite
            fm[id_] = 1
    tm = rm + len(fm)
    print "Total elements: ", tm
    print "False positives: ", len(fm)
    print "% false positive:",100 - (rm * 100.0 / tm)