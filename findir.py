from Bio.Seq import Seq
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from subprocess import Popen, PIPE
import os
from Bio.SeqUtils import GC
import queue
import logging

def complex_enough(seq):
    complexity = lcc_simp(seq.upper())
    if complexity < 1.25:
        return False
    gc = GC(seq.upper())
    if gc < 15 or gc > 95:
        return False
    return True

def makelog(stri, do_print=True):
    if do_print:
        print(stri)
    logging.debug(stri)

def findIR(q, args,l_lock, candidates, perc_seq, last_perc_seq):
    while True:
        try:
            seq, seq_fs, split_index, record_id, seq_len, count = q.get(timeout=15)
        except queue.Empty:
            break
        splited_len = len(seq)
        seq_rc = str(Seq(seq).reverse_complement())
        complexity = lcc_simp(seq.upper())
        if not complex_enough(seq):
            q.task_done()
            continue
        record_q = SeqRecord(Seq(seq), id = record_id)
        record_s = SeqRecord(Seq(seq_rc), id = record_id + "_rc")
        query_filename = "results/" + args.jobname + "/tmp/query" + str(record_id + "_" + str(split_index))+".tmp"
        subject_filename = "results/" + args.jobname + "/tmp/subject" + str(record_id + "_" + str(split_index))+".tmp"
        SeqIO.write(record_q, query_filename, "fasta")
        SeqIO.write(record_s, subject_filename, "fasta")
        cmd_list = [
        'blastn',
        '-query',query_filename,
        '-subject',subject_filename,
        '-reward','2',
        #'-max_target_seqs','1',
        '-penalty','-4',
        '-word_size','4',
        #'-ungapped',
        #'-evalue','140',
        '-strand','plus',
        #'-soft_masking','false',
        #'-dust','no',
        '-outfmt',"6 sstart send qstart qend score length mismatch gaps gapopen nident"]
        p = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        out,err = p.communicate()
        if err:
            print(split_index, record_id,seq_len)
            makelog("BLASTN error: %s" % (err, ) )
        os.remove(query_filename)
        os.remove(subject_filename)
        lines = out.splitlines()
        for row in lines:
            row = row.split()
            sstart = int(row[0])
            send = int(row[1])
            qstart = int(row[2])
            qend = int(row[3])
            score = int(row[4])
            length = int(row[5])
            mismatch = int(row[6])
            gaps = int(row[7])

            #filter valids IR
            if length < args.align_min_len:
                continue

            #subject transform cause it was reversed
            sstart = splited_len - sstart 
            send = splited_len - send

            #obtain IR sequences
            seq_q = seq[qstart:qend]
            seq_q_prime = seq[send:sstart]

            qstart, qend = min(qstart, qend),max(qstart, qend)
            sstart, send = min(sstart, send),max(sstart, send)

            #organice positions
            ir_start = min(qstart,qend,sstart,send)
            ir_end = max(qstart,qend,sstart,send)
            
            #calculate length
            ir_len = ir_end - ir_start

            #length constraints
            if ir_len > args.mite_max_len:
                continue
            if ir_len < args.mite_min_len:
                continue
            #move in genome, split index
            #ir_seq = seq[ir_start:ir_end]

            #again validate complexity, a value of 1 means only two different nucleotides are present
            if not complex_enough(seq_q):
                continue

            if not complex_enough(seq_q_prime):
                continue

            #validate TSD outside TIRs
            i = args.tsd_max_len
            valid_tsd = False
            while i >= args.tsd_min_len:
                tsd_one = seq_fs[ir_start - i + args.FSL:ir_start + args.FSL]
                tsd_two = seq_fs[ir_end + args.FSL:ir_end + i + args.FSL]
                if tsd_one.lower() == tsd_two.lower():
                    valid_tsd = True
                    mite_pos_one = ir_start + args.FSL
                    mite_pos_two = ir_end + args.FSL
                    tsd_in = 'no'
                    break
                i -= 1

            #validate TSD inside TIRs
            #TSDs cannot be a large part of TIRs
            if not valid_tsd:
                i = args.tsd_max_len
                while i >= args.tsd_min_len:
                    tsd_one = seq_fs[ir_start + args.FSL:ir_start + i + args.FSL]
                    tsd_two = seq_fs[ir_end - i + args.FSL:ir_end + args.FSL]
                    if tsd_one.lower() == tsd_two.lower():
                        valid_tsd = True
                        mite_pos_one = ir_start + args.FSL + i
                        mite_pos_two = ir_end  + args.FSL - i
                        tsd_in = 'yes'
                        break
                    i -= 1
            #"no tsd"
            if not valid_tsd:
                continue

            ir_seq = seq_fs[mite_pos_one:mite_pos_two]
            if not complex_enough(ir_seq):
                continue
            #ir_seq = seq_fs[mite_pos_one - args.FSL:mite_pos_two + args.FSL]
            ir_len = mite_pos_two - mite_pos_one

            fs_start = max(0,mite_pos_one - args.FSL)
            flanking_seq_left = seq_fs[fs_start:mite_pos_one]
            flanking_seq_right = seq_fs[mite_pos_two:mite_pos_two + args.FSL]

            if len(flanking_seq_right) < args.FSL or len(flanking_seq_left) < args.FSL:
                continue

            #calculate positions in full sequence
            mite_start_full = mite_pos_one + split_index - args.FSL 
            mite_end_full = mite_pos_two + split_index - args.FSL
            
            #new_element = (mite_start_full, mite_end_full, ir_seq, record_id, ir_len, seq_q, seq_q_prime, tsd_one, tsd_in,flanking_seq_left,flanking_seq_right,length,'','unfiltered','')
            new_element = {
                'start': mite_start_full,
                'end': mite_end_full, 
                'end': mite_end_full, 
                'seq': ir_seq, 
                'record': record_id, 
                'mite_len': ir_len, 
                'tir1_start': mite_start_full,
                'tir1_end': mite_start_full + length,
                'tir2_start': mite_end_full - length,
                'tir2_end': mite_end_full,
                'tir1_seq': seq_q, 
                'tir2_seq': seq_q_prime,
                'tsd': tsd_one,
                'tsd_in': tsd_in,
                'fs_left': flanking_seq_left,
                'fs_right': flanking_seq_right,
                'tir_len': length,
            }
            with l_lock:
                index = "%s_%i" % (record_id, (count)) 
                #we don't want overlapped TIRs, save the broader
                is_nested = False
                has_nested = False
                for curr_count in range(count - 1, count + 2):
                    curr_index = "%s_%i" % (record_id, (curr_count)) 
                    if curr_index in candidates:
                        for candidate in candidates[curr_index]:
                            #if new element TIR is nested in other TIR
                            if new_element['start'] >= candidate['tir1_start'] and \
                                new_element['start'] <= candidate['tir1_end'] and \
                                new_element['end'] >= candidate['tir2_start'] and \
                                new_element['end'] <= candidate['tir2_end'] and \
                                new_element['record'] == candidate['record']:
                                is_nested = True
                            #if other element is nested inside new element TIR
                            if candidate['start'] >= new_element['tir1_start'] and \
                                candidate['start'] <= new_element['tir1_end'] and \
                                candidate['end'] >= new_element['tir2_start'] and \
                                candidate['end'] <= new_element['tir2_end'] and \
                                candidate['record'] == candidate['record']:
                                has_nested = True
                                candidates[curr_index].remove(candidate)
                if has_nested or not is_nested:
                    if not index in candidates:
                        candidates[index] = []
                    candidates[index].append(new_element)
                        
            curr_perc = int(split_index * 100 / seq_len)

            if not record_id in perc_seq or not record_id in last_perc_seq:
                perc_seq[record_id] = curr_perc
                last_perc_seq[record_id] = curr_perc
                makelog(record_id + " " + str(curr_perc) + "%")

            if perc_seq[record_id] - last_perc_seq[record_id] >= 10:
                makelog(record_id + " " + str(curr_perc) + "%")
                last_perc_seq[record_id] = curr_perc

            perc_seq[record_id] = curr_perc

        q.task_done()
