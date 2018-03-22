#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO

def FAspliter(sequence, seqs, outdir):
    """Well exclude some ids from a multifasta file
    """
    fasta_seq = SeqIO.parse(sequence, 'fasta')
    buffer_seqs = {}
    for record in fasta_seq:
        for seq in seqs:
            if seq in record.id:
                buffer_seqs.setdefault(seq, []).append(record)
    if outdir[-1] != "/":
        outdir += "/"
    for key, value in buffer_seqs.iteritems():
        SeqIO.write(value, outdir + key + ".fasta", "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()#pylint: disable=invalid-name
    parser.add_argument("-s", "--sequence", help="Sequence file (.fasta)", required=True)
    parser.add_argument("-l", "--seqs", help="Filter sequence that contains one of this seqs", nargs='+')
    parser.add_argument("-o", "--outdir", help="Output directory", default='')
    args = parser.parse_args()#pylint: disable=invalid-name
    FAspliter(args.sequence, args.seqs, args.outdir)
