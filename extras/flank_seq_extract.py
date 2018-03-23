# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-f", "--fasta", help="Fasta file", required=True)
parser.add_argument("-g", "--genome", help="Genome file", required=True)
parser.add_argument("-o", "--out", help="Out file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

buffer_rec = []
fasta_seq = SeqIO.parse(args.fasta, 'fasta')
for record in fasta_seq:
    list_id = record.id.split('|')
    start = int(list_id[2])
    end = int(list_id[3])
    curr_id = list_id[1]
    start -= 50
    end += 50
    list_id[2] = str(start)
    list_id[3] = str(end)
    list_id = '|'.join(list_id)
    genome = SeqIO.parse(args.genome, 'fasta')
    for record_g in genome:
        if record_g.id == curr_id:
            clean_seq = ''.join(str(record_g.seq).splitlines())
            seq = clean_seq[start:end]
            record = SeqRecord(Seq(seq), id=list_id, description=record.description)
            buffer_rec.append(record)
SeqIO.write(buffer_rec, args.out , "fasta")