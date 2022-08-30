import sys
import statistics
from bed_parser import parse_bed

field_names = ["chrom","chromStart","chromEnd",
                   "name","score","strand",
                   "thickStart","thickEnd",
                   "itemRGB","blockCount",
                   "blockSizes","blockStarts"]
# making field_names global

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