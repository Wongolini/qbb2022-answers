from vcfParser import *
import sys 
'''
Create an analysis script that will import your VCF parser, 
load both the 1KGP genomes (random_snippet.vcf) and dbSNP 
(dbSNP_snippet.vcf) VCF files, and then annotate the 1KGP 
variants with the IDs from the dbSNP file, when appropriate. 
Your script should include the following features

 - Creates a dictionary of positions and IDs from the dbSNP file
 - Replaces the ID in each record from random_snippet.vcf with the correct label, if it exists, from your dbSNP dictionary
 - Finds and reports the number of random_snippet.vcf records that do not have a corresponding ID in your dbSNP ID dictionary

 vcf_parser outputs:
 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', {'AC': , 'AN',: ...}
'''
class VCF_Data:
    def __init__(self, vcf):
        # first thing to do is create separate lists
        self.vcf = vcf
        self.chrom = [_[0] for _ in self.vcf if _[0]!='CHROM']
        self.pos = [int(_[1]) for _ in self.vcf if _[1]!='POS']
        self.id = [_[2] for _ in self.vcf if _[2]!='ID']
        self.ref = [_[3] for _ in self.vcf if _[3] != 'REF']
        self.alt = [_[4] for _ in self.vcf if _[4] != 'ALT']
        self.qual = [_[5] for _ in self.vcf if _[5] != 'QUAL']
        self.filter = [_[6] for _ in self.vcf if _[6] != 'FILTER']  #
        self.info = [_[7] for _ in self.vcf if _[7] != 'INFO']    # info dictionary
        self.chrom_pos = set(zip(self.chrom, self.pos))
################################################################################
# Be careful knowing what the format of the parsed vcf file is. The first line
# in the list is the header line, which you clearly know. However, the INFO
# field isn't actually "INFO", but a dictionary of keys and descriptions for the
# fields in the info column. Same thing for the optional FORMAT field. Also, you
# aren't loading the FORMAT or sample columns here, so just be aware that they
# are being lost
################################################################################

    def wrap(self):
        # wrap chrom to ids
        # wrap ids to pos
        # {ids:{
        #        ('chrom#',pos),...,
        #               }
        # check ids are unique
        assert len(self.id) == len(set(self.id)) # according to docs, duplicate IDs are not allowed in VCF
        annotation_dict = dict(zip(self.chrom_pos, self.id))
        self.annotation_dict = annotation_dict
    
    def replace_ids(self, annotations):
        for i,(pos,chrom,id) in enumerate(zip(self.pos,self.chrom,self.id)):
            try:
                id = annotations[tuple([chrom,pos])]
                self.id[i] = id # replace newly annotated ids
            except:
                continue
################################################################################
# Just a note that "id" is a built-in python function that is getting over-
# written here
################################################################################
    
    def export_vcf(self):
        self.final_vcf = [self.chrom, self.pos, self.id,
                      self.ref, self.alt, self.qual,
                      self.filter, self.info]
    
class Annotate:
    def __init__(self, REF, UN):
        self.ref = REF      # dictionary of reference names with positions
        self.needs_map = UN # list of positions that require annotation
    
    def annotate(self):
        annotated = 0
        annotation_dict = {}
        for data in self.needs_map:
            if data in self.ref.keys():
                annotation_dict[data] = self.ref[data]
                annotated+=1
        unnanotated = len(self.needs_map)-annotated
        print('{} Annotated ids'.format(annotated))
        print('{} missing ids'.format(unnanotated))
        self.annotation_dict = annotation_dict

def main(vcf_1, vcf_2):
    d1 = parse_vcf(vcf_1)
    d2 = parse_vcf(vcf_2)

    V1 = VCF_Data(d1)
    V1.wrap()
    V2 = VCF_Data(d2)

    A = Annotate(V1.annotation_dict, V2.chrom_pos)
    A.annotate()
    V2.replace_ids(A.annotation_dict)
    V2.export_vcf()

    return V2.final_vcf

if __name__ == "__main__":
    vcf_1 = sys.argv[1] # annotated vcf
    vcf_2 = sys.argv[2] # vcf that requires annotation -- no ids.
    finalvcf=main(vcf_1, vcf_2)
    
################################################################################
# This is a really nice solution. I'm impressed by the class implementation.
# I'm a little surprised by the use of the Annotation class, since you've
# already built in the needed functionality into the "replace_ids" method. 
# All you're missing is a counter under the "except" statement to track missing
# IDs. Overall, excellent code! Also, just note that you were also supposed to
# comment the VCF parser code (but not the header parsing section).
################################################################################

