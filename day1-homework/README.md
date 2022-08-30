Exercise 1.
Use -v command in awk to define bash variable as an awk variable.

Exercise 2.
awk required -v command to initialize bash variable as an awk variable
-F '\t' ensured delimited by tab
We used state 2 (flanking TSS) as putative promoter regions
C-->T is a pyrimadine to pyrimadine transition vs transversion. These mutations are more likely.
Also, an increase of T's lowers the gibbs free energy and can put the DNA in a more readily open configuration for DNA binding proteins.

Exercise 3.
First error requires tab delimited. We can fix this with awk command -v OFS='\t'
Second error is Error: Sorted input specified, but the file variants.bed has the following out of order record
chr21	5218156	5218157
fix this by piping the proper sort to the variants.bed file

There are 10293 variants
wc -l exercise3.out

There are 200 unique genes
cut -f 7 exercise3.out | sort | uniq | wc -l

10293/200 variants per gene average
51.465

Exercise 4.
Needed to change paths to correctly access data
grep command needed to have correct options, -Fxvf 
-F, --fixed-strings
             Interpret pattern as a set of fixed strings (i.e., force grep to behave as fgrep).
-f file, --file=file
             Read one or more newline separated patterns from file.  Empty pattern lines match every input line.  Newlines are not considered part of a
             pattern.  If file is empty, nothing is matched.
-v, --invert-match
             Selected lines are those not matching any of the specified patterns.

-x, --line-regexp
             Only input lines selected against an entire fixed string or regular expression are considered to be matching lines.


