# -*- coding: utf-8 -*-
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-e", "--embl", help="EMBL file", required=True)
parser.add_argument("-f", "--fasta", help="FASTA file", required=True)
parser.add_argument("-o", "--out", help="Out file", required=True)

args = parser.parse_args()#pylint: disable=invalid-name

content_file = open(args.embl, 'r')
content = content_file.read()
content_split = content.split('//\n')

print len("total",content_split)
nonaut = []
for con in content_split:    
    file_id = con.strip().split('\n')[0].replace('ID','',1).strip().split(' ')[0]
    if file_id == '':
        continue
    if ('nonautonomous' in con or 'Nonautonomous' in con or 'non-autonomous' in con or 'Non-autonomous' in con or 'MITE' in con) and \
        not 'Retrotransposon' in con:
        nonaut.append(file_id)
print len("nonaut",nonaut)

nonaut_mite = []
nonaut_mite_rec = []
fasta_seq = SeqIO.parse(args.fasta, 'fasta')
for record in fasta_seq:
    if record.id in nonaut:
        clean_seq = ''.join(str(record.seq).splitlines())
        if len(clean_seq) < 800:
            nonaut_mite.append(record.id)
            nonaut_mite_rec.append(record)

print len("nonaut < 800",nonaut_mite)
SeqIO.write(nonaut_mite_rec, args.out , "fasta")