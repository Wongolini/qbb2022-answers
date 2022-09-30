for file in ./*.sam;
    do
        outfile=${file%".sam"}
        outfile=${outfile#"./"}.bam
        samtools sort -@ 4 -O BAM -o $outfile $file
    done
