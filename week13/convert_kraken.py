#!/usr/bin/env python

#USAGE: python convert_kraken.py KRAKEN_sample_input_file.kraken sample_name
#EXAMPLE: python make_krona_compatible.py KRAKEN/SRR492183.kraken SRR492183

import sys
import os 

input = sys.argv[1]
f = open(input)
filename = input.split('/')[-1]

sample = filename.strip('.kraken')
out = open(sample + "_krona.txt", "w")


for line in f:
    fields = line.strip('\r\n').split(';')
    out.write("\t".join(fields[1:]) + "\n")