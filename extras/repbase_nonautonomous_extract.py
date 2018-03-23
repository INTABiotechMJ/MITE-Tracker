# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import argparse

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-o", "--out", help="Out file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name

for fn in os.listdir('data/embl_splitted'):
    buffer_rec = []
    print fn
    fasta_seq = SeqIO.read('data/embl_splitted/' + fn, 'embl')
    for record in fasta_seq:
        print record