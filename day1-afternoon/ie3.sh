inputfile=$1
nucoi=$2
grep -v "#" $inputfile | awk '{if ($4==$2) {print $5}}' | sort | uniq -c