from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_simp
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import Popen, PIPE
import os 
from Bio.Blast.Applications import NcbiblastnCommandline
import argparse
import time
import sys

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-g", "--genome", help="Genome file in fasta format", required=True)
parser.add_argument("-o","--outfile", help="Output fasta file", required=True)
#parser.add_argument("-c","--cores", help="Max number of processes to use simultaneously", type=int, default=1)
#parser.add_argument("--tir_perc_len", help="TIR lenght proportion relative to MITE length (100% means half the MITE minus the TSD)", type=int, default=30)
#parser.add_argument("--tir_min_len", help="TIR min lenght (will overwrite tir_perc_len)", type=int, default=0) #TODO remove min word
#parser.add_argument("--tsd_min_len", help="TSD min lenght", type=int, default=2)
#parser.add_argument("--tsd_max_len", help="TSD max lenght", type=int, default=10)
#parser.add_argument("--mite_min_len", help="MITE min lenght", type=int, default=50)
#parser.add_argument("--mite_max_len", help="MITE max lenght", type=int, default=800)
#parser.add_argument("--flank_seq_len", help="Flanking sequence length", type=int, default=30)
#parser.add_argument("-m","--method", help="0:pairwise(default) - 1:mismatches", type=int, default=0)
#parser.add_argument("-d","--database", help="Database for classification", default=False)
args = parser.parse_args()#pylint: disable=invalid-name

MIN_TSD_LEN = 2
MAX_TSD_LEN = 10
MITE_MAX_LEN = 600
MITE_MIN_LEN = 30

#-soft_masking false -dust no
start_time = time.time()
windows_size = MITE_MAX_LEN / 2
fasta_seq = SeqIO.parse(args.genome, 'fasta')
mites = [] 
count = 1
record_count = 1
for record in fasta_seq:
    print record.id, record_count
    record_count += 1
    porc_ant = 0
    split_index = MAX_TSD_LEN
    clean_seq = ''.join(str(record.seq).splitlines())
    seq_len = len(clean_seq)
    while split_index < seq_len - windows_size:
        seq = clean_seq[split_index:split_index + windows_size]
        splited_len = len(seq)
        if lcc_simp(seq) < 0.5: #Discard really low complexity MITE
            split_index += windows_size
            continue
        seq_rc = str(Seq(seq).reverse_complement())
        record_q = SeqRecord(Seq(seq), id = record.id)
        record_s = SeqRecord(Seq(seq_rc), id = record.id + "_rc")
        SeqIO.write(record_q, "query", "fasta")
        SeqIO.write(record_s, "subject", "fasta")
        cmd_list = ['blastn','-query',"query",'-subject',"subject",
        '-reward','2',
        '-max_target_seqs','1',
        '-penalty','-4',
        '-word_size','7','-ungapped',
        '-evalue','10','-strand',"plus",
        '-soft_masking','false' ,'-dust','no',
        '-outfmt',"6 sstart send qstart qend score length mismatch gaps gapopen nident"]
        #]
        blast_process = Popen(cmd_list, stdout=PIPE, stderr=PIPE)
        out,err = blast_process.communicate()

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
            if length < 10:
                continue

            #subject transform cause it was reversed
            transform_send = splited_len - sstart 
            transform_sstart = splited_len - send 

            #obtain IR sequences
            seq_q = seq[qstart-1:qend]
            seq_q_prime = seq[transform_sstart:transform_send+1]
            seq_s = Seq(seq_q_prime).reverse_complement()

            #organice positions
            #subject larger than query
            if (transform_send + transform_sstart) / 2 >  (qend + qstart) / 2:
                ir_1_start = qstart
                ir_1_end = qend
                ir_2_start = transform_sstart
                ir_2_end = transform_send
            else:
                ir_1_start = transform_sstart
                ir_1_end = transform_send
                ir_2_start = qstart
                ir_2_end = qend
            #validate distance between IR
            if ir_2_end - ir_1_start > MITE_MAX_LEN:
                continue
            if ir_2_end - ir_1_start < MITE_MIN_LEN:
                continue
            #again validate complexity
            if lcc_simp(seq_q) < 0.95 or lcc_simp(seq_s) < 0.95:
                continue

            #validate TSD outside TIRs
            i = MAX_TSD_LEN
            valid_tsd = False
            while i >= MIN_TSD_LEN:
                tsd_one = seq[ir_1_start - i:ir_1_start]
                tsd_two = seq[ir_2_end:ir_2_end + i]
                if tsd_one.lower() == tsd_two.lower():
                    valid_tsd = True
                    mite_pos_one = ir_1_start - i
                    mite_pos_two = ir_2_end + i
                    tsd_in = 'no'
                    break
                i -= 1
            #validate TSD inside TIRs
            if not valid_tsd:
                i = MAX_TSD_LEN
                while i >= MIN_TSD_LEN:
                    tsd_one = seq[ir_1_start:ir_1_start+i]
                    tsd_two = seq[ir_2_end-i:ir_2_end]
                    if tsd_one.lower() == tsd_two.lower():
                        valid_tsd = True
                        mite_pos_one = ir_1_start
                        mite_pos_two = ir_2_end 
                        tsd_in = 'yes'
                        break
                    i -= 1
            if not valid_tsd:
                continue
            mite = seq[mite_pos_one:mite_pos_two]
            mite_len = mite_pos_two - mite_pos_one

            #calculate positions in full sequence
            mite_start_full = mite_pos_one + split_index
            mite_end_full = mite_pos_two + split_index 
            mites.append( (mite, record.id, mite_start_full, mite_end_full, mite_len, tsd_one, tsd_in, ) )
            porc = split_index * 100 / seq_len
            if porc - porc_ant > 5:
                print "%i%%" % (porc, ), 
                sys.stdout.flush()
                porc_ant = porc
            count += 1
        split_index += windows_size
    #endfor record fasta

#discover nested mites
count = 1
mites_arr = []

for mite in mites:
    nested = False
    for mite_2 in mites:
        #print mite
        #print mite_2
        #print
        if mite[1] == mite_2[1] and mite[2] > mite_2[2] and mite_2[3] > mite[3]: 
            #mite_2 is nested in mite_1
            nested = True
    if not nested:
        params = (mite[1], mite[2], mite[3], mite[4], mite[5], mite[6] )
        description = "SEQ:%s START:%i END:%i MITE_LEN:%i TSD:%s TSD_IN:%s" % (params)
        mite_seq_rec = SeqRecord(Seq(mite[0]), id='MITE_' + str(count), description = description)
        mites_arr.append(mite_seq_rec)
        count += 1
print "Found %i MITES" % (count - 1, )

SeqIO.write(mites_arr, args.outfile, "fasta")
elapsed_time = time.time() - start_time
print elapsed_time