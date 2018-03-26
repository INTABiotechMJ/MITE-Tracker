
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import OrderedDict
import time
import argparse

#defaults
parser = argparse.ArgumentParser()#pylint: disable=invalid-name
parser.add_argument("-d", "--detect", help="Fasta file", required=True)
parser.add_argument("-o", "--out", help="Out file", required=True)
args = parser.parse_args()#pylint: disable=invalid-name


start_time = time.time()

file_repbase = args.detect
#file_repbase = "f1.fa"
fasta_rec = SeqIO.parse(file_repbase, 'fasta')
fasta_index = SeqIO.index(file_repbase, 'fasta')
fasta_index = OrderedDict(sorted(fasta_index.items()))
elem_not_nested = []
for record in fasta_rec:
    record_1 = record.id
    seq_1, start_1, end_1, tir_len_1, tsd_len_1 = record_1.split("|")
    nested = False
    for record_2 in fasta_index:
        if record_1 == record_2:
            continue
        seq_2, start_2, end_2, tir_len_2, tsd_len_2 = record_2.split("|")
        if seq_1 == seq_2 and  \
            int(start_2) >= int(start_1)  and  \
            int(start_2) <= int(start_1) + int(tir_len_1)  and  \
            int(end_2) <= int(end_1)  and  \
            int(end_2) >= int(end_1) - int(tir_len_1):
            nested = True
            break
    if not nested:
        elem_not_nested.append(record)
print len(fasta_index)
print len(elem_not_nested)

print("--- %s seconds ---" % (time.time() - start_time))
SeqIO.write(elem_not_nested, args.out , "fasta")
#36029
#22763
#36029
#25960