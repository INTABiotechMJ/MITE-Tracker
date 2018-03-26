import pandas as pd

#TREP Comparison

for i in ("data/tracker_trep.csv", "data/detect_trep.csv"):
    print "-----File: ",i
    df = pd.read_csv(i, delimiter="\t")
    df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]
    rm, fm = 0, {}
    for id_ in pd.unique(df.qseqid):
        mite = df[df.qseqid==id_].iloc[0].sseqid
        #print mite
        if mite.startswith('XXX'):
            print id_
        if mite.startswith('DT'):
            rm += 1
        else:
            print "False positive: ",id_, '-> '  ,mite
            fm[id_] = 1
    tm = rm + len(fm)
    print "Total elements: ", tm
    print "False positives: ", len(fm)
    print "% false positive:",100 - (rm * 100.0 / tm)
    print "% covered", len(pd.unique(df.sseqid))