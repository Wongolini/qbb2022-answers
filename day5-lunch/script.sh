python day5_lunch.py aau1043_dnm.csv
join -t "	" <(sort aau1043_parental_age.tsv -k1,1r) <(sort day5_lunch.txt -k1,1r) > aau1043_merged.tsv