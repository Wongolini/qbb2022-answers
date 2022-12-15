from Bio import SeqIO
from Bio import SeqUtils
import sys 

if __name__ == "__main__":
    fasta=sys.argv[1]
    nts = 0
    for record in SeqIO.parse(fasta,'fasta'):
        nts+=len(str(record.seq))
        print(record.id)
        print(SeqUtils.GC(record.seq))
    print(fasta)
    print(nts)