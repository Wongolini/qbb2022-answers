for x in chr11:1900000-2800000 chr14:100700000-100990000 chr15:23600000-25900000 chr20:58800000-58912000;

do
   medaka_variant -p -d -i nanopore/methylation.bam -f hg38.fa -r $x -o ${x%:*}_medaka_variant_out
done
