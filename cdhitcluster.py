# -*- coding: utf-8 -*-
import itertools
import os
import logging

def makelog(stri, do_print=True):
    if do_print:
        print(stri)
    logging.debug(stri)

def loadcluster(cluster_file):
    cluster_file, cluster_dic = open(cluster_file), {}
    # parse through the cluster file and store the cluster name + sequences in the dictionary
    cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
    for cluster in cluster_groups:
        name = cluster.next().strip()
        seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.next()]
        cluster_dic[name[1:]] = seqs
    return cluster_dic

def filtercluster(cluster_dic, minimum,positions):
    filtered_dic = {}
    for cluster in set(cluster_dic.keys()):
        #let's get first only valid groups to avoid unnecesary processing
        if len(cluster_dic[cluster]) >= minimum:
            #group into overlapped (overlapped groups count as one individual)
            #ie should have more than minimum elements non overlapped
            cluster_positions = [v for k,v in positions.items() if k in cluster_dic[cluster]]
            #res = df.loc[df['candidate_id'].isin(cluster_dic[cluster])]
            #res = res[['start','end']]
            #cluster_positions = [tuple(x) for x in res.values]
            merged_overlaped = merge_overlap(cluster_positions)
            if len(merged_overlaped) >= minimum:
                filtered_dic[cluster] = cluster_dic[cluster]
    return filtered_dic

def merge_overlap(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged

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
    family_number = 1
    for cluster_id, cluster_seqs in cluster_dic.items():
        for seq_id in cluster_seqs:
            description, sequence = seq_list[seq_id]
            description += " family:" + str(family_number)
            if not cluster_id == last_cluster:
                filter_file.write('{0}\n'.format('-'*10))
                last_cluster = cluster_id
            filter_file.write('{0}\n{1}\n'.format(description, '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])))
            family_number += 1
    # close the filtered results file
    filter_file.close()

def cluster(file_names, positions, min_copy_number, df, FSL):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio import pairwise2
    from subprocess import Popen, PIPE
    from collections import OrderedDict
    max_len = int(df[['len']].max())
    min_len = int(df[['len']].min())
    num_files = 10
    sep_size = ((max_len - min_len) / num_files) 
    margin = sep_size * 0.25
    last = min_len
    count = 0
    filtered_clusters = {}
    irs_seqs_total = []
    print min_len, max_len, sep_size, margin
    valid_seqs = []
    named_seqs = {}
    print '*' * 20
    for i in range(1, num_files + 1):
        makelog("Clustering file " + str(i))
        curr = (i * sep_size) + min_len # separations
        print last, curr
        current_seqs = df[(df.len >= last) & (df.len <= curr)]
        irs_seqs = []
        for index, row in current_seqs.iterrows():
            params = (row.record, row.start, row.end, row.tsd, row.tsd_in, row.len)
            description = "SEQ:%s START:%i END:%i TSD:%s TSD_IN:%s MITE_LEN:%i" % (params)
            named_seqs[row.candidate_id] = (description, row.seq)
            ir_seq_rec = SeqRecord(Seq(row.seq), id=row.candidate_id, description = description)
            irs_seqs.append(ir_seq_rec)
            #irs_seqs_total.append(ir_seq_rec)
        file_candidates_partial = file_names['file_candidates_partial_prefix'] + str(i) + ".fasta"
        SeqIO.write(irs_seqs, file_candidates_partial , "fasta")
        candidates_partial_repr = file_names['candidates_partial_repr_prefix'] + str(i) + ".fasta"
        cmd_list = [
        './cdhit/cd-hit-est',
        '-i', file_candidates_partial,
        '-o', candidates_partial_repr,
        '-c', '0.80','-n','7','-d','0','-T','0','-aL','0.8','-s','0.8','-M','0']
        p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        for c in iter(lambda: p.stdout.read(), ''):
            continue
            makelog(c)
        clusters_dic = loadcluster(candidates_partial_repr + ".clstr")
        new_clusters = filtercluster(clusters_dic, min_copy_number,positions)
        valid_seqs += [item for sublist in new_clusters.values() for item in sublist]
        last = (i * sep_size) - margin + min_len

    #write only valid sequences
    valid_seq_records = []
    for k,v in named_seqs.items():
        if k in valid_seqs:
            description, seq = v
            ir_seq_rec = SeqRecord(Seq(seq), id=k, description=description)
            valid_seq_records.append(ir_seq_rec)

    SeqIO.write(valid_seq_records, file_names['file_candidates_fasta'] , "fasta")

    makelog("Clustering valid sequences using cdhit")
    cmd_list = [
    './cdhit/cd-hit-est',
    '-i',file_names['file_candidates_fasta'],
    '-o',file_names['file_candidates_cluster'],
    '-c', '0.80','-n','7','-d','0','-T','0','-aL','0.8','-s','0.8','-M','0']
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    for c in iter(lambda: p.stdout.read(), ''):
        continue
        makelog(c)
    #out,err = p.communicate()
    makelog("Clustering done")

    clusters_dic = loadcluster(file_names['file_candidates_cluster'] + ".clstr")
    filtered_clusters = filtercluster(clusters_dic, min_copy_number, positions)
    unique_clusters = set(filtered_clusters.keys())
    num_clusters = len(unique_clusters)

    #loop through clusters
    for current_cluster in unique_clusters:
        #search candidates for that cluster
        #all possible 2-combinations of candidates
        candidates = filtered_clusters[current_cluster]
        combinations = [(x,y) for x,y in itertools.combinations(candidates, 2)]
        dist_fs = {}
        for seq_id in combinations:
            x,y = seq_id

            if not x in candidates or not y in candidates:
                continue
            cand_x = df[(df.candidate_id == x)]
            cand_y = df[(df.candidate_id == y)]

            #if they're partially overlapped, ignore flanking sequence comparison
            if cand_x.iloc[0].end >= cand_y.iloc[0].start and cand_y.iloc[0].end >= cand_x.iloc[0].start:
                continue

            fs_right_1 = cand_x.iloc[0].fs_right
            fs_left_1 = cand_x.iloc[0].fs_left

            fs_right_2 = cand_y.iloc[0].fs_right
            fs_left_2 = cand_y.iloc[0].fs_left
            
            score_r1_r2 = pairwise2.align.localms(fs_right_1, fs_right_2, 1, -1, -1, -1,score_only=True)
            score_l1_l2 = pairwise2.align.localms(fs_left_1, fs_left_2, 1, -1, -1, -1,score_only=True)
            score_r1_l2 = pairwise2.align.localms(fs_right_1, fs_left_2, 1, -1, -1, -1,score_only=True)
            score_r2_l1 = pairwise2.align.localms(fs_right_2, fs_left_1, 1, -1, -1, -1,score_only=True)
            max_score = max(score_r1_r2,score_l1_l2,score_r1_l2,score_r2_l1)

            if max_score == []:
                max_score = 0
            max_score /= FSL
            #todo validate scoring
            if max_score < 0.5:
                dist_fs[x] = 1
                dist_fs[y] = 1
        if len(dist_fs) < min_copy_number:
            df.loc[df['candidate_id'].isin(filtered_clusters[current_cluster]), 'status'] =  'low_cn_flank_seq'
            del filtered_clusters[current_cluster]

    #again to remove < MIN_COPY_NUMBER elements
    #filtered_clusters = clusterutils.filtercluster(filtered_clusters, args.min_copy_number, positions, df, 'low_copy_number_2')
    ordered_cluster = OrderedDict(sorted(filtered_clusters.items(), key=lambda t: t[1]))

    makelog("Clusters: " + str(len(filtered_clusters)) + " writing sequences")

    fasta_seq = SeqIO.parse(file_names['file_candidates_fasta'], 'fasta')
    buffer_rec = []
    for record in fasta_seq:
        for clus, seqs in filtered_clusters.items():
            if record.id in seqs:
                df.loc[df['candidate_id'] == record.id, 'cluster'] =  clus
                df.loc[df['candidate_id'] == record.id, 'status'] =  'valid'
                record.description = "%s CLUSTER:%s" % (record.description, clus)
                buffer_rec.append(record)
                continue

    SeqIO.write(buffer_rec, file_names['all_file'] , "fasta")
    df.to_csv(file_names['file_candidates_csv'], index=False)
    cluster2seq(ordered_cluster, file_names['file_candidates_fasta'], file_names['families_file'] )
    return True