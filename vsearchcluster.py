# -*- coding: utf-8 -*-
import itertools
import logging
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqUtils import GC
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

def filtercluster(cluster_dic, minimum, candidates):
    filtered_dic = {}
    for cluster in set(cluster_dic.keys()):
        #let's get first only valid groups to avoid unnecesary processing
        if len(cluster_dic[cluster]) >= minimum:
            #group into overlapped (overlapped groups count as one individual)
            #ie should have more than minimum elements non overlapped
            cluster_positions = {}
            for k in cluster_dic[cluster]:
                if isinstance(candidates[k]['start'], int) and isinstance(candidates[k]['end'], int):
                    if not candidates[k]['record'] in cluster_positions:
                        cluster_positions[candidates[k]['record']] = []
                    cluster_positions[candidates[k]['record']].append( (candidates[k]['start'], candidates[k]['end']) )
            #cluster_candidates = [v for k,v in candidates.items() if k in cluster_dic[cluster]]
            #for candidate in cluster_candidates:
            #    cluster_positions.append( (candidate['start'], candidate['end']) )
            total_len = 0
            for k,v in cluster_positions.items():
                overlapped_cluster = merge_overlap(v)
                total_len += len(overlapped_cluster)
            if total_len >= minimum:
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

def cluster2seq(cluster_dic, candidates, outfile):
    # open the filter_output file
    filter_file = open(outfile, 'w')
    # loop through the sequence list
    last_cluster = None
    #family_number = 1
    for cluster_id, cluster_seqs in cluster_dic.items():
        for seq_id in cluster_seqs:
            cand_id, description, sequence = candidates[seq_id]['id'], candidates[seq_id]['description'],candidates[seq_id]['seq']
            #header = ">%s %s FAMILY:%s" % (cand_id, description, str(family_number)) 
            header = ">%s %s" % (cand_id, description) 
            if not cluster_id == last_cluster:
                filter_file.write('{0}\n'.format('-'*10))
                last_cluster = cluster_id
            filter_file.write('{0}\n{1}\n'.format(header, '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])))
        #family_number += 1
    # close the filtered results file
    filter_file.close()


def complex_enough(seq):
    complexity = lcc_simp(seq.upper())
    if complexity < 1.25:
        return False
    gc = GC(seq.upper())
    if gc < 15 or gc > 95:
        return False
    return True

def cluster(file_names, candidates, min_copy_number, FSL, workers):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio import pairwise2
    from subprocess import Popen, PIPE
    from collections import OrderedDict
    import os, shutil
    import math

    makelog("Clustering")
    cmd_list = [
    './vsearch-2.7.1/bin/vsearch',
    '--cluster_fast',file_names['file_candidates_fasta'],
    #'--consout',file_names['file_representative'],
    '--threads',str(workers),
    '--strand','both',
    '--clusters',file_names['file_temp_cluster'],
    '--iddef','1',
    '-id', '0.8']
    p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
    out,err = p.communicate()
    #for stdout_line in iter(popen.stdout.readline, ""):
    #    yield stdout_line
    #popen.stdout.close()
    #return_code = popen.wait()
    ##if return_code:
    #    raise subprocess.CalledProcessError(return_code, cmd)

    #for c in iter(lambda: p.stdout.read(), ''):
        #makelog(c)
    makelog("Clustering done")
    makelog("Filtering clusters")
    #count for minimum file length
    clusters_dic = {}
    list_dir = os.listdir(file_names['file_temp_cluster_dir'])
    makelog("Initial clusters: %i" % (len(list_dir),))
    for fn in list_dir:
        if os.path.isfile(file_names['file_temp_cluster_dir'] + fn):
            fh = open(file_names['file_temp_cluster_dir'] + fn)
            n = 0
            for line in fh:
                if line.startswith(">"):
                    n += 1
                    id_seq = line[1:line.find(" ")]
                    if fn in clusters_dic:
                        clusters_dic[fn].append(id_seq)
                    else:
                        clusters_dic[fn] = [id_seq]
            fh.close()
    shutil.rmtree(file_names['file_temp_cluster_dir'])
    
    #        os.unlink(file_names['file_temp_cluster_dir'] + fn)
    #        if n < args.min_copy_number:
    #            df.loc[df['candidate_id'] == 'id_seq', 'status'] =  'low_cn'
    #            os.remove(fn)
    #            continue
    #clusters_dic = loadcluster(cluster_candidates_file + ".clstr")
    filtered_clusters = filtercluster(clusters_dic, min_copy_number, candidates)
    unique_clusters = set(filtered_clusters.keys())
    num_clusters = len(unique_clusters)
    #loop through clusters
    for current_cluster in unique_clusters:
        #search candidates for that cluster
        #all possible 2-combinations of candidates
        candidates_in_cluster = filtered_clusters[current_cluster]
        #porc_of_clusters = int(math.ceil(len(candidates_in_cluster) * 0.4))
        #new_min_copy_number = max(min_copy_number,porc_of_clusters)
        new_min_copy_number = min_copy_number
        sum_diff_fs_cluster = 0
        for x in candidates_in_cluster:
            totally_different_fs = True
            cand_x = candidates[x]
            fs_right_1 = cand_x['fs_right'].upper()
            fs_left_1 = cand_x['fs_left'].upper()
            if fs_left_1 == '' or fs_right_1 == '' or not isinstance(fs_left_1,str) or not isinstance(fs_right_1,str):
                totally_different_fs = False
                continue
            if not complex_enough(fs_right_1) or not complex_enough(fs_left_1):
                totally_different_fs = False
                continue

            at_least_one = False
            for y in candidates_in_cluster:
                cand_y = candidates[y]
                if cand_x['candidate_id'] == cand_y['candidate_id']:
                    continue
                # R1 x R2
                # L1 x L2
                # L1RC x R2
                # R1RC x L2
                fs_right_2 = cand_y['fs_right'].upper()
                fs_left_2 = cand_y['fs_left'].upper()
                #some MITE could be at the end or begining of the sequence and this not having flanking seqs
                if fs_right_2 == '' or fs_left_2 == '':
                    continue
                #empty strings in some versions of pandas are returned as nan, so we make sure the flanking seqs are strings
                if not isinstance(fs_right_2,str) or not isinstance(fs_left_2,str):
                    continue
                if not complex_enough(fs_right_2) or not complex_enough(fs_left_2):
                    continue


                fs_left_1_rc = Seq(fs_left_1).reverse_complement()
                fs_right_1_rc = Seq(fs_right_1).reverse_complement()
                
                #calculate scores
                score_r1_r2 = pairwise2.align.localms(fs_right_1, fs_right_2, 1, -1, -1, -1,score_only=True)
                score_l1_l2 = pairwise2.align.localms(fs_left_1, fs_left_2, 1, -1, -1, -1,score_only=True)
                score_l1rc_r2 = pairwise2.align.localms(fs_left_1_rc, fs_right_2, 1, -1, -1, -1,score_only=True)
                score_r1rc_l2 = pairwise2.align.localms(fs_right_1_rc, fs_left_2, 1, -1, -1, -1,score_only=True)
                #get max score
                #max_score = max(score_r1_r2,score_l1_l2)
                max_score = max(score_r1_r2,score_l1_l2,score_l1rc_r2,score_r1rc_l2)
                if max_score == []:
                    max_score = 0
                max_score /= FSL
                at_least_one = True
                if max_score > 0.5:
                    totally_different_fs = False
                    break
               
            if totally_different_fs and at_least_one:
                sum_diff_fs_cluster += 1
            if sum_diff_fs_cluster >= new_min_copy_number:
                break

        if sum_diff_fs_cluster < new_min_copy_number:
            #makelog(' '.join(filtered_clusters[current_cluster]) + " filtered by flanking sequence")
            del filtered_clusters[current_cluster]
        #else:
        #    makelog(' '.join(filtered_clusters[current_cluster]) + " not filtered by flanking sequence")

    #again to remove < MIN_COPY_NUMBER elements
    #filtered_clusters = filtercluster(filtered_clusters, args.min_copy_number, positions, df, 'low_copy_number_2')
    ordered_cluster = OrderedDict(sorted(filtered_clusters.items(), key=lambda t: t[0]))

    makelog("Clusters: " + str(len(filtered_clusters)))
    buffer_rec = []
    #import ipdb; ipdb.set_trace()
    #for candidate in candidates.values():
    count = 1
    family_number = 1#ordered_cluster.keys().index(clus)
    buffer_nr = []
    output_gff = open(file_names['file_gff'],"w") 
    output_gff.write("##gff-version 3\n")

    for clus, seqs in ordered_cluster.items():
        one_per_family = False
        tsd_family = []
        for seq in seqs:
            candidate = candidates[seq]
            candidate['id'] = "MITE_T_%s|%s|%s|%s|%s|%s|F%s" % (str(count),candidate['record'],candidate['start'],candidate['end'],candidate['tsd'],candidate['tir_len'],family_number)
            candidate['description'] = "%s CANDIDATE_ID:%s" % (candidate['description'], candidate['candidate_id'].split('|')[0])
            record = SeqRecord(Seq(candidate['seq']), id=candidate['id'], description=candidate['description'])
            buffer_rec.append(record)
            
            write_row =  '\t'.join([candidate['record'], 'MITE_Tracker','MITE',str(candidate['start']), str(candidate['end']),'.','+','.','ID='+candidate['id'] ]) 
            output_gff.write(write_row + '\n')

            tsd_family.append(candidate['tsd'])
            if not one_per_family:
                one_per_family = True
                record_family = record
                buffer_nr.append(record_family)
            count += 1
        from statistics import mode,StatisticsError
        try:
            tsd_consensus = mode(tsd_family)
            record_family.description = '%s COMMON_TSD:%s' % (record_family.description, tsd_consensus)
        except StatisticsError:
            pass
        family_number += 1


    SeqIO.write(buffer_nr, file_names['file_representative'] , "fasta")
    SeqIO.write(buffer_rec, file_names['all_file'] , "fasta")
    cluster2seq(ordered_cluster, candidates, file_names['families_file'])