from Bio import SeqIO
import sys
contigfile = sys.argv[1]
ref_file = sys.argv[2]

genome_length = 0

for record in SeqIO.parse(ref_file,'fasta'):
	genome_length+=len(record.seq)

n50_bp=0
n50=0
for record in SeqIO.parse(contigfile,'fasta'):
	if n50_bp>=.5*genome_length:
		break
	else:
		n50_bp+=len(record.seq)
		n50+=1
print('n50=',n50)
