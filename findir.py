from Bio.Seq import Seq
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from subprocess import Popen, PIPE
import os
import Queue
import logging

def makelog(stri, do_print=True):
    if do_print:
        print(stri)
    logging.debug(stri)

def findIR(q, args,l_lock, irs, perc_seq, last_perc_seq):
    
    while True:
        try:
            seq, seq_fs, split_index, record_id,seq_len = q.get(timeout=5)
        except Queue.Empty:
            break
        splited_len = len(seq)
        seq_rc = str(Seq(seq).reverse_complement())
        complexity = lcc_simp(seq)
        if complexity < 1:
            q.task_done()
            continue
        record_q = SeqRecord(Seq(seq), id = record_id)
        record_s = SeqRecord(Seq(seq_rc), id = record_id + "_rc")
        query_filename = "tmp/query" + str(record_id + "_" + str(split_index))+".tmp"
        subject_filename = "tmp/subject" + str(record_id + "_" + str(split_index))+".tmp"
        SeqIO.write(record_q, query_filename, "fasta")
        SeqIO.write(record_s, subject_filename, "fasta")
        cmd_list = [
        'blastn',
        '-query',query_filename,
        '-subject',subject_filename,
        '-reward','2',
        #'-max_target_seqs','1',
        '-penalty','-4',
        '-word_size','7',
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
        else:
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
            if lcc_simp(seq_q) <= 1.2:
                continue

            if lcc_simp(seq_q_prime) <= 1.2:
                continue

            #validate TSD outside TIRs
            i = args.tsd_max_len
            valid_tsd = False
            while i >= args.tsd_min_len:
                tsd_one = seq_fs[ir_start - i + args.FSL:ir_start + args.FSL]
                tsd_two = seq_fs[ir_end + args.FSL:ir_end + i + args.FSL]
                if tsd_one.lower() == tsd_two.lower():
                    valid_tsd = True
                    mite_pos_one = ir_start - i
                    mite_pos_two = ir_end + i
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
                        mite_pos_one = ir_start
                        mite_pos_two = ir_end 
                        tsd_in = 'yes'
                        break
                    i -= 1
            #"no tsd"
            if not valid_tsd:
                continue

            ir_seq = seq_fs[mite_pos_one + args.FSL :mite_pos_two  + args.FSL ]
            ir_len = mite_pos_two - mite_pos_one

            flanking_seq_left = seq_fs[mite_pos_one:mite_pos_one + args.FSL]
            flanking_seq_right = seq_fs[mite_pos_two+args.FSL:mite_pos_two + args.FSL + args.FSL]
           
            #calculate positions in full sequence
            mite_start_full = mite_pos_one + split_index
            mite_end_full = mite_pos_two + split_index 
            new_element = (mite_start_full, mite_end_full, ir_seq, record_id, ir_len, seq_q, seq_q_prime, tsd_one, tsd_in,flanking_seq_left,flanking_seq_right,length,'')

            with l_lock:
                irs[record_id + "_" + str(mite_start_full) + "_" + str(mite_end_full)] = new_element

            curr_perc = split_index * 100 / seq_len

            if not record_id in perc_seq or not record_id in last_perc_seq:
                perc_seq[record_id] = curr_perc
                last_perc_seq[record_id] = curr_perc
                makelog(record_id + " " + str(curr_perc) + "%")

            if perc_seq[record_id] - last_perc_seq[record_id] >= 5:
                makelog(record_id + " " + str(curr_perc) + "%")
                last_perc_seq[record_id] = curr_perc

            perc_seq[record_id] = curr_perc

        q.task_done()
