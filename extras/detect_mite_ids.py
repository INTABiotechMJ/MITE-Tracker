
from Bio import SeqIO
from collections import OrderedDict
import time
start_time = time.time()

file_repbase = "../detectmite/rice.miteSet.1.fasta"
#file_repbase = "f1.fa"
fasta_index = SeqIO.index(file_repbase, 'fasta')
fasta_index = OrderedDict(sorted(fasta_index.items()))
elem_not_nested = []
for record_1 in fasta_index:
    print record_1