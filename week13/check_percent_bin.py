from Bio import SeqIO 
import sys 


def count_nts(query):
    contigs = 0
    for record in SeqIO.parse(query,'fasta'):
        contigs+=1
    return contigs

if __name__ == "__main__":
    bin = sys.argv[1]
    assembly = sys.argv[2]
    print(bin)
    print(100*(count_nts(bin)/count_nts(assembly)), '%')

