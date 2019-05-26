import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="BLASTn outfmt 6 file", required=True)
parser.add_argument("-n", "--name", help="Program name (required for second column)", required=True)
parser.add_argument("-f", "--feature", help="Feature name (required for third column)", required=True)
parser.add_argument("-o", "--output", help="gff3 file name output", required=True)
args = parser.parse_args()

#read blast output
df = pd.read_csv(args.input,sep="\t", header=None)
df.columns = ['qseqid','sseqid','qstart','qend','sstart','send','score','length','mismatch','gaps','gapopen','nident','pident','evalue','qlen','slen','qcovs']

#add some dummy columns
df['source'] = args.name
df['frame'] = '.'
df['score'] = '.'
df['strand'] = '+'
df['feature'] = args.feature
df = df[['sseqid','source','feature','sstart','send','score','strand','frame','qseqid']]

#gff3 -> ['seqname' , 'source' , 'feature' , 'start' , 'end' , 'score' , 'strand' , 'frame' , 'attribute']
#make it happend
df.to_csv(args.output, sep="\t", header=None, index=None)
