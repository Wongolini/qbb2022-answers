for file in fastq/*;
    do
        filename=$(basename $file);
        samplename=${filename%".fastq"}
        R="@RG\tID:$samplename\tSM:$samplename"
        outfile=$samplename.sam
        bwa mem -t 4 -R $R reference_genome/sacCer3.fa $file > $outfile
    done
