#!/usr/bin/env python3
import sys

def parse_bed(fname):
    def numeric_list_parser(string_list):
        # special function to delimit lists into integer lists.
        # issue that the last value is , eg. [1,2,3,] --> ['1', '2', '3', ',']
        hold_list = string_list.rstrip(',').split(',')
        numeric_list = [int(x) for x in hold_list[:len(hold_list)]]
################################################################################
# Nice use of a function here!
# You don't actually need the [:len(hold_list)] as using a list as an iterator
# will always step through each item once with built-in bounds checking.
################################################################################
        return numeric_list
    try:
        fs = open(fname, 'r')
    except:
        raise FileNotFoundError("That file doesn’t appear to exist")
    bed = []
    field_types = [str, int, int, 
                   str, float, str,
                   int, int, set,
                   int, int, list]
    field_names = ["chrom","chromStart","chromEnd",
                   "name","score","strand",
                   "thickStart","thickEnd",
                   "itemRGB","blockCount","blockSizes",
                   "blockStarts"]
    # illegal field is blockCount and blockSizes
    # 0 chrom is string
    # 1 chromStart is int
    # 2 chromEnd is int
    # 3 name is str
    # 4 score is int 0-1000
    # 5 strand is str +/-
    # 6 thickStart int 
    # 7 thickEnd int
    # 8 itemRGB is a list, but set to a tuple?
    # 9 blockCount
    # 10 blockSizes (list, comma sep)
    # 11 blockStarts is a comma-separated list of block starts
    for i, line in enumerate(fs):
        if line.startswith("#"):
            continue
        fields = line.rstrip().split()
        fieldN = len(fields)
        if fieldN < 3 or fieldN > 12:
            print(f"Line {i} appears malformed", file=sys.stderr)
            continue
        
        for j in range(min(len(field_types), len(fields))):
            if j == 11 or j == 10: # field 11 is a comma delimited list. Parse into int list.
                try:
                    fields[j] = numeric_list_parser(fields[j]) #fix string lists with comma delimit

                except:
                    print('Broken field {}, index {}'.format(field_names[j], j))
                    continue
################################################################################
# This continue will simply go to the next value of j, not the next record. Not
# sure if this is intentional. You could have a boolean like "valid" and set it
# to False if you encounter a bad record, followed by break.
################################################################################
            elif j == 4:    #field 4 can sometimes have '.'. Ignore such cases as None.
                try:
                    fields[j] = field_types[j](fields[j])
                except:
                    if fields[j] == '.':
                        fields[j] = None
                    
            elif j == 8:    #field 8 is RGB. Cast to list.
                try:
                    RGB_hold = fields[j].split(',')
                    RGB_set = [int(x) for x in RGB_hold] #check int
                    fields[j] = set(RGB_set)
################################################################################
# I'm surprised that you didn't use your numeric_list_parse function here:
#
# fields[j] = set(numeric_list_parse(fields[j]))
################################################################################
                except:
                    print('Field 8 not composed of integers on line {}'.format(i))
                    fields[j] = None
            
            else:
                try:    # Generic exceptions
                    fields[j] = field_types[j](fields[j])
                except:
                    print('Unkown issue with field {}, index {}'.format(field_names[j],j))
                    print('File had value: {}'.format(fields[j]))
                    print('Expected {}'.format(field_types[j]))
  
        if len(fields[11]) == fields[9] and len(fields[11]) == len(fields[10]): # ensure proper values for fields 9,10,11
            bed.append(fields)
################################################################################
# This will only work if there are 12 fields in the bed. If there are less, it's
# going to throw an IndexOutOfRange error. Also you're not excluding records 
# with issues prior to this (except for fieldN < 3 | fieldN > 12). Finally, you
# need a check for fieldN == 10 | fieldN == 11, since those are both invalid.
################################################################################
        else:
            continue 

    fs.close()
    return bed

if __name__ == "__main__":
    fname = sys.argv[1]
    bed = parse_bed(fname)

################################################################################
# All in all, this is really good code. It's very close to right. I appreciate
# how readable it is, comments, when they appear, are informative, and the use 
# of a function for reused code is excellent.
################################################################################

 