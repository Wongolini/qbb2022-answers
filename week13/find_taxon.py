import numpy as np
from Bio import SeqIO 
import sys

class TallyTaxon:
    def __init__(self,taxons):
        self.taxons=taxons # 2D array; col 1 = contig name; col 2 = strain

    def shape_taxon(self):
        taxon_list  = self.taxons[:,1]
        taxon_array = np.empty((len(taxon_list), 12),dtype=object)
        for i,t in enumerate(taxon_list):
            t = t.split(';')
            t = [str(_ or 'NA') for _ in t]
            taxon_array[i,:len(t)] = t
            taxon_array[i,len(t):] = ['NA']*(12-len(t))
        print(taxon_array[:,11])
        self.taxon_array = taxon_array

    
    def tally_taxon(self):
        taxon_headers = ['root','cellular organism',
                         'Super Kingdom','Group', 
                         'Phylum','Class','Order',
                         'Family','Genera','Species',
                         'Strain','subspecies']
        '''
        root;
        cellular organisms;
        (super kingdom)
        (unranked group?);
        (phylum)
        (class)
        (order)
        (family)
        (Genera)
        (species)
        (strain)
        (subsp)
        '''
        taxon_dict = {}
        for i,col in enumerate(self.taxon_array.T):
            members,member_counts = np.unique(col,return_counts=True)
            taxon_dict[taxon_headers[i]] = dict(zip(members, member_counts/len(self.taxon_array)))

        self.taxon_dict = taxon_dict

if __name__ == "__main__":
    fasta = sys.argv[1]
    kraken = sys.argv[2]
    k = np.loadtxt(kraken,delimiter='\t',dtype=object)
    indices = []
    for record in SeqIO.parse(fasta,'fasta'):
        if record.id in k[:,0]:
            i = np.where(record.id==k[:,0])[0][0]
            indices.append(i)
            #print(record.id,k[i,1][0])
    taxons = k[indices]
    T = TallyTaxon(taxons)
    T.shape_taxon()
    T.tally_taxon()
    print(T.taxon_dict)
    # add statistics separate names by semicolon


