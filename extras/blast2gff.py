import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="BLASTN outfmt 6 file", required=True)
parser.add_argument("-n", "--name", help="Program name (required for second column)", required=True)
parser.add_argument("-o", "--output", help="gff3 file name output", required=True)
args = parser.parse_args()

#read blast output
df = pd.read_csv(args.input,sep="\t", header=None)
df.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

#add some dummy columns
df['source'] = args.name
df['frame'] = '.'
df['score'] = '.'
df['strand'] = '+'
df['feature'] = 'TE'

output_gff = open(args.output,"w") 
output_gff.write("##gff-version 3\n")
count = 1
for k,v in df.iterrows():
    start = str( min(v.sstart, v.send) )
    end = str( max(v.sstart, v.send) )
    strand = '+' if start < end else '-'
    description = 'ID=%s;occurrence=%i' % ( v.qseqid, count, )
    params = [v.sseqid, args.name,'MITE',start, end, '.', strand, '.', description]
    write_row =  '\t'.join(params) 
    output_gff.write(write_row + '\n')
    count += 1
print "Epa! wrote %i elements" % ( (count-1) , )