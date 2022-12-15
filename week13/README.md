```
mkdir step0_workdir
mkdir step1_workdir

cd step0_workdir
for x in ../convert_kraken.py ../metagenomics_data/step0_givendata/KRAKEN/*; do python ../convert_kraken.py $x; done

# add Krona tools to ~/.profile
# export KRONA=/path/to/krona_tools/scripts/
mkdir krona_files
mv S* krona_files #move all the sample krona text files into ./krona_files

# run KRONA ImportText.pl on all files
arg=''
for x in krona_files/*; do 
    file=$(basename $x)
    sample=${file%_krona.txt}; 
    arg+=$x,$sample' '; 
done; 

perl $KRONA/ImportText.pl $arg -o step0.html -q
```

# Question 1: In your README, briefly comment on the trends you see in the gut microbiota throughout the first week.

1. Staphylococcaceae:
    a. Typically compose less than a quarter of the microbiome for all samples
    b. epidermis is present in all samples; with multiple strains in all samples but SRR492197.

2. Enterococci 
    a. Enterococci make up majority of the microbiome for all samples

3. Actinobacteria
    a. Sometimes are not present in samples.

# Question 2: In your README, comment on what metrics in the contigs could we use to group them together?
- --p1
    - p1 = probability cutoff for bin seeding. A higher p1-->higher specificity per bin
- --minProb 
    - minimum probability for binning consideration and controls sensitivity for bin creation.
- Can also use preset params --verysensitive or --sensitive or --specific, etc.
I do not expect the newborn's microbiome to be as diverse as the mother's or an older human. We could play with two extremes, --verysensitive, --specific, --superspecific and compare. We can then decide if superspecific is overfitting the microbiome, or if verysensitive is good enough to model the microbiome.

```
bwa index /Users/cmdb/qbb2022-answers/week13-homework/step0_workdir/bins_dir/assembly.fasta
bwa_mem.py assembly.fasta fastqs/

for x in ./*.sam; do 
    samtools sort $x -o ${x%.sam}.bam; 
done

conda activate metabat2
mkdir bins_dir
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
metabat2 -i assembly.fasta -a depth.txt -o bins_dir/bin

```
# Qusetion 3a:
I got 4 bins

# Question 3B
```
(base) [~/qbb2022-answers/week13-homework/step0_workdir $]for x in ./bins_dir/*; do python check_percent_bin.py $x assembly.fasta; done
./bins_dir/bin.1.fa
1.9254204240799415 %
./bins_dir/bin.2.fa
0.7799171338045333 %
./bins_dir/bin.3.fa
1.4867170363148916 %
./bins_dir/bin.4.fa
0.14623446258835 %
```

# Question 3C
for x in ./bins_dir/*; do python nt_size.py $x; done > bins_dir/basic_genome_stats.txt

./bins_dir/bin.1.fa
2,266,021
./bins_dir/bin.2.fa
587,549
./bins_dir/bin.3.fa
2,920,155
./bins_dir/bin.4.fa
2,862,852

All bins contain roughly the amount of nucleotides expected in a bacterial genome except bin 2. All the bins also contain similar GC content.
>> see "bins_dir/basic_genome_stats.txt to see GC content

# Question 3D
Run the bins on BUSCO or CheckM to see what taxon and taxon depth the bins can be predicted on. checkM and BUSCO both look for marker genes speific to a taxon group and can check the purity.

One could also check the level of clipped reads and discordant reads could also be signs of contmaination (for strains with overlapping repeat regions).

Check the GC content of each scaffold in the bin and see if they agree with each other.

```
for bin in bins_dir/*; do 
    binname=$(basename $bin); 
    python find_taxon.py $bin ../metagenomics_data/step0_givendata/KRAKEN/assembly.kraken > ${binname%.fa}_taxon.txt; 
done
```

# Question 4
    A.
        1. Bin 1 = 'Species': {'Staphylococcus aureus': 0.012658227848101266, 
                               'Staphylococcus epidermidis': 0.9873417721518988}, 
                    'Strain': {'Staphylococcus aureus subsp. aureus': 0.012658227848101266, 
                               'Staphylococcus epidermidis ATCC 12228': 0.3291139240506329, 
                               'Staphylococcus epidermidis RP62A': 0.6329113924050633}
            --> likely Staph epidermis.

        2. Bin 2 = Species: {'Enterococcus faecalis': 0.08, 
                              'Staphylococcus aureus': 0.2,                    
                              'Staphylococcus epidermidis': 0.24, 
                              'Staphylococcus haemolyticus': 0.4}, 
                   Strain: {'Enterococcus faecalis V583': 0.04, 
                              'Staphylococcus aureus subsp. aureus': 0.2, 
                              'Staphylococcus epidermidis ATCC 12228': 0.2, 
                              'Staphylococcus haemolyticus JCSC1435': 0.4}
            --> Decent chance Staph epidermis & haemolyticus with Enterococcus faecalis

        3. Bin 3  = Species: {'Anaerococcus prevotii': 0.01639344262295082, 
                              'Cutibacterium avidum': 0.21311475409836064, 
                              'Prevotella': 0.01639344262295082, 
                              'Staphylococcus aureus': 0.16393442622950818, 
                              'Staphylococcus epidermidis': 0.32786885245901637, 
                              'Staphylococcus haemolyticus': 0.03278688524590164}, 
                    'Strain': {'Anaerococcus prevotii DSM 20548': 0.01639344262295082, 
                               'Cutibacterium avidum 44067': 0.21311475409836064, 
                               'Prevotella dentalis': 0.01639344262295082, 
                               'Staphylococcus aureus M1': 0.03278688524590164, 
                               'Staphylococcus aureus subsp. aureus': 0.13114754098360656, 
                               'Staphylococcus epidermidis ATCC 12228': 0.14754098360655737, 
                               'Staphylococcus epidermidis RP62A': 0.14754098360655737, 
                               'Staphylococcus haemolyticus JCSC1435': 0.03278688524590164}
            --> Diverse mix of Staph with Cutibacterium

        4. Bin 4 = 'Species': {'Enterococcus faecalis': 1.0}, 
                    'Strain': {'Enterococcus faecalis D32': 0.16666666666666666, 
                               'Enterococcus faecalis OG1RF': 0.3333333333333333, 
                               'Enterococcus faecalis V583': 0.3333333333333333, 
                               'Enterococcus faecalis str. Symbioflor 1': 0.16666666666666666}
            --> Enterococcus faecalis
            
    B. Use checkM or BUSCO to see what percentage of marker genes for a taxon are detected.