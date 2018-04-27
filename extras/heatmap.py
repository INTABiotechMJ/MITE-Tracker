"""
Grabs a blastn -outfmt 6 output and plots a heatmap by chromosome 
"""
import argparse
import pandas as pd

parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genome", help="Genome description gff", required=True)
parser.add_argument("-b", "--blast", help="allhits file outfmt 6", required=True)
parser.add_argument("-s", "--step", help="Divide chromosomes into steps", default=5000, type=int)
parser.add_argument("-o", "--out", help="Out file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

#loads genome description
df_genome = pd.read_csv(args.genome, index_col=False, sep='\t', header=None)
df_genome.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

df_blast = pd.read_csv(args.blast, index_col=False, sep='\t', header=None)
df_blast.columns = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore',]

result = []
max_hits_count = 0
max_chr_start = ""
max_chr_end = ""
max_chr = ""
print(len(df_genome))
for k,  chromosome in df_genome.iterrows():
    df_chromosome = df_blast[(df_blast.sseqid == chromosome.seqname)]
    print chromosome.seqname
    print(len(df_chromosome))
    for i in range(0, chromosome.end, args.step):
        start = i
        end = i + args.step
        if end > chromosome.end:
            continue
#            start = chromosome.end - args.step 
#            end = chromosome.end + 1
        #print start, end
        hits_count = len(df_chromosome[(df_chromosome.sstart >= start) & (df_chromosome.sstart <= end) & (df_chromosome.pident >= 80)])
        if hits_count > max_hits_count:
            max_chr_start = start
            max_chr_end = end
            max_chr = chromosome.seqname
            max_hits_count = hits_count
            print '->',chromosome.seqname,max_hits_count
        #print(chromosome.seqname, start, end, hits_count)
        position = start / 1000000.0
        result.append( [chromosome.seqname,position, hits_count] )
        #print result
df = pd.DataFrame(result)
df.columns = ['chromosome','position','hits']
result = pd.pivot_table(data=df,
                    index='position',
                    values='hits',
                    columns='chromosome')
import seaborn as sns
import matplotlib.pyplot as plt
from numpy import arange
x = arange(25).reshape(5, 5)
cmap = sns.dark_palette("#2ecc71", as_cmap=True)

#fmt="g",
cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma']),
         ('Sequential', [
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
         ('Sequential (2)', [
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper']),
         ('Diverging', [
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
         ('Qualitative', [
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']),
         ('Miscellaneous', [
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]
print('st',max_chr_start)
print('end',max_chr_end)
print('chr',max_chr)
print('hit',max_hits_count)
import numpy as np
for cm in cmaps:
    for cmap in cm[1]:

        #num_ticks = 10
        # the index of the position of yticks
        #yticks = np.linspace(0, len(result) - 1, num_ticks, dtype=np.int)
        # the content of labels of these yticks
        #yticklabels = [result[idx] for idx in yticks]
        yticklabels = []
        for i in range(45):
            yticklabels.append(str(i) + "Mb")
            #yticks.append(i * 1.0)
        ax = sns.heatmap(result, cmap=cmap,yticklabels=7)
        #ax = sns.heatmap(result, cmap=cmap)
        #ax.set_yticks(np.arange(0, 45, step=5),tuple(yticklabels))

        #plt.show()
        plt.savefig('heatmap/heatmap_'+str(args.step)+'_'+cmap, dpi=800,bbox_inches='tight')
        plt.clf()
        exit()
