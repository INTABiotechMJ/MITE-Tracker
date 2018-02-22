# -*- coding: utf-8 -*-
import itertools
import os

def loadcluster(cluster_file):
    cluster_file, cluster_dic = open(cluster_file), {}
    # parse through the cluster file and store the cluster name + sequences in the dictionary
    cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
    for cluster in cluster_groups:
        name = cluster.next().strip()
        seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.next()]
        cluster_dic[name[1:]] = seqs
    return cluster_dic

def filtercluster(cluster_dic, minimum):
    filtered_dic = {}
    for cluster in set(cluster_dic.keys()):
        if len(cluster_dic[cluster]) >= minimum:
            filtered_dic[cluster] = cluster_dic[cluster]
    return filtered_dic

def cluster2seq(cluster_dic, fasta, outfile):
    seq_file, seq_list, = open(fasta), {}
    # parse through the fasta file and obtain the sequence
    seq_groups = (x[1] for x in itertools.groupby(seq_file, key=lambda line: line[0] == '>'))
    for header in seq_groups:
        description = header.next().strip()
        header = description.split(" ")[0]
        sequence = ''.join(seq_line.strip() for seq_line in seq_groups.next())
        seq_list[header[1:]] = (description,sequence)
    # close the sequence file
    seq_file.close()
    # open the filter_output file
    filter_file = open(outfile, 'w')
    # loop through the sequence list
    last_cluster = None
    for cluster_id, cluster_seqs in cluster_dic.items():
        for seq_id in cluster_seqs:
            description, sequence = seq_list[seq_id]
            if not cluster_id == last_cluster:
                filter_file.write('{0}{1}\n'.format('-'*10,cluster_id ))
                last_cluster = cluster_id
            filter_file.write('{0}\n{1}\n'.format(description, '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])))
    # close the filtered results file
    filter_file.close()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()#pylint: disable=invalid-name
    parser.add_argument("-c", "--cluster", help="cdhit cluster file (.clstr)", required=True)
    parser.add_argument('-m', '--minimum', help='The minimum cluster size.', required=True, type=int)

    args = parser.parse_args()#pylint: disable=invalid-name
    filtered_dic = filtercluster(args.cluster, args.minimum)
    for k,i in filtered_dic.items():
        print k,i
