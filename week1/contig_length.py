#!/bin/bash/env python
from Bio import SeqIO
import sys

file=sys.argv[1]
seq_lens = []
for record in SeqIO.parse(file,'fasta'):
	print('{}: length: {}'.format(record.id,len(record.seq)))
	seq_lens.append(len(record.seq))
print(seq_lens)
