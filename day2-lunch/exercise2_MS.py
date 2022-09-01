import sys
import statistics
from bed_parser import parse_bed


def main(file):
    bed_data = parse_bed(file)
    gene_names = set([gn[3] for gn in bed_data])
    
    gene_exons_dict = {gn:0 for gn in gene_names}
    for line in bed_data:
        gn=line[3]
        gene_exons_dict[gn]+=line[9]
    median_exons = statistics.median(list(gene_exons_dict.values()))
    print("Median number of exons per gene: {}".format(median_exons))

if __name__ =="__main__":
    file = sys.argv[1]
    main(file)

################################################################################
# Very nice script. Excellent logic creating the dictionary full of zeros.
# While it is unneccessary in this case since there is only one line per gene
# in our bed file, this code can handle the situation of multiple lines for the
# same gene. My only real feedback is that by creating the gene_names set and
# the gene_exons_dict ahead of time, you are introducing 2 new for loops into
# your code. Since for loops are slow in python, you could easily speed up your
# code by doing something like this:
#
# gene_exons_dict = {}
# for line in bed_data:
#     gn=line[3]
#     if gn not in gene_exons_dict:
#         gene_exons_dict[gn] = 0
#     gene_exons_dict[gn] += line[9]
#
# or even easier:
#
# gene_exons_dict = {}
# for line in bed_data:
#     gn=line[3]
#     gene_exons_dict.set_default(gn, 0) # This only adds gn:0 if gn not in dict
#     gene_exons_dict[gn] += line[9]
#
################################################################################