
from Bio import SeqIO
from collections import OrderedDict
import time
start_time = time.time()

file_repbase = "../results/rice_wz_6_mcn_3/all.fasta"
#file_repbase = "f1.fa"
fasta_index = SeqIO.index(file_repbase, 'fasta')
fasta_index = OrderedDict(sorted(fasta_index.items()))
elem_not_nested = []
elem_nested = []
for record_1 in fasta_index:
    mid,seq_1, start_1, end_1, tsd, tir_len_1, fn = record_1.split("|")
    nested = False
    for record_2 in fasta_index:
        if record_1 == record_2:
            continue
        mid,seq_2, start_2, end_2, tsd, tir_len_2, fn = record_2.split("|")
        if seq_1 == seq_2 and  \
            int(start_2) >= int(start_1)  and  \
            int(start_2) <= int(start_1) + int(tir_len_1)  and  \
            int(end_2) <= int(end_1)  and  \
            int(end_2) >= int(end_1) - int(tir_len_1):
            #print record_1, record_2
            nested = True
            elem_nested.append(record_1)
            elem_nested.append(record_2)
            break
    if not nested:
        elem_not_nested.append(record_1)
print len(fasta_index)
print len(elem_not_nested)
print elem_nested
print("--- %s seconds ---" % (time.time() - start_time))