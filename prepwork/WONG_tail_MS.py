#!/usr/bin/python3
import argparse
import sys 

def print_tail(file,numlines,maxlines,skip_header):
	line_num=0
	begin=maxlines-numlines # begin reading from this line
	print(begin)
	with open(file) as f:
		for line in f.readlines():
			if (str(line[0])=='#' and skip_header==True):
				continue
			else:
				line_num+=1
				if line_num in range(begin,maxlines):
					print(line.strip('\n'))
				else:
					pass

def print_head(file,numlines,maxlines,skip_header):
	# print_head and print_tail could be one function
	# but maybe this is more readable for grading
	begin=0
	line_num=0
	if numlines>maxlines:
		numlines=maxlines 
	else:
		pass 
	with open(file) as f:
		for line in f.readlines():
			if (str(line[0])=='#' and skip_header==True):
				continue
			elif line_num in range(begin,numlines):
				print(line.strip('\n'))
				line_num+=1
			else:
				break

def countlines(file,skip_header): # safety check
	i=0
	with open(file) as f:
		for line in f.readlines():
			if skip_header==True and str(line)[0]=='#':
					continue 
			else:
				i+=1
	return i 

def main(args):

	file=None  		  # file
	numlines=10		  # default always print 10 lines
	tail=False 		  # to read in reverse
	skip_header=False # Show lines beginning with '#' by default

	try:
		file = args.file
	except:
		sys.exit()
	try:
		numlines=args.numlines
	except:
		pass
	try:
		tail=args.tail
	except:
		pass
	try:
		skip_header=args.skip_headers
	except:
		pass
	
	maxlines=countlines(file,skip_header) # safety even if python doesn't care
	if tail==True:
		print_tail(file,numlines,maxlines,skip_header)
	else:
		print_head(file,numlines,maxlines,skip_header)
	

if __name__ == "__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('-f','--file',help='File',required=True, type=str)
	parser.add_argument('-n','--numlines',help='Number of lines to print',required=False,type=int)
	parser.add_argument('-t','--tail',help='Tail, print in reverse order',required=False,action='store_true')
	parser.add_argument('-s','--skip_headers',help='Skip header lines. By default script always prints headers',
						required=False,action='store_true')
	args=parser.parse_args()
	main(args)
	# code is sloppy and redundant but w/e
	# code allows you to choose between head or tail function 
	# -n option to print n lines
	# skip_header option allows you to skip header. By default prints headers


# Overall, this is an impressive script that goes above and beyond what was
# expected. I'm thrilled that you used argparse as well as try/except.
# When you check if the line starts with "#", you are checking line[0], which
# works here since each line should end with a newline character. But it would
# be safer to use line.startwith("#"), since this won't throw an OutOfBoundsError
# if the line is empty. Also, in the print_tail function, you don't need the
# else: pass statement since if the conditional is false, nothing will happen anyway.
# Those are just some tiny inconsequential things. As I said, this looks great. - Mike