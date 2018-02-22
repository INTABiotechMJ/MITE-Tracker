import cdhitutils
from collections import OrderedDict


filtered_clusters = cdhitutils.filtercluster("results/rice_chr01_t2/mites.cand.cluster.clstr",3)
num_clusters = len(set(filtered_clusters.values()))
ordered_cluster = OrderedDict(sorted(filtered_clusters.items(), key=lambda t: t[1]))

for seq_id, cluster_id in ordered_cluster.items():
    df[(df.record == row.record)]

seqs = cdhitutils.cluster2seq("sda.fa","results/rice_chr01_t2/mites.candidates.fasta",ordered_cluster)