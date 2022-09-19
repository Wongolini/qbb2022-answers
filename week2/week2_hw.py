#!/usr/bin/env python
'''
Write a script to perform global alignment between two sequences using a given scoring matrix and gap penalty. Your script will take four inputs:

A FASTA-style file containing two sequences to align
A text file containing the scoring matrix you’d like to use for this alignment
The penalty for gaps in your alignment
The filepath to write your alignment to

You’ll run your script twice:

Align the CTCF DNA transcript sequences from human and mouse using the HOXD70 scoring matrix and a gap penalty of 300.
Align the CTCF amino acid sequences from human and mouse using the BLOSUM62 scoring matrix and a gap penalty of 10.

NOTE: The DNA sequences are fairly long, and as such the DNA alignment may take a few minutes to run. 
We recommend testing your code with the protein alignment first (or even just a couple of small test sequences), 
and then running the DNA alignment when you’re confident it’s working.
'''
import numpy as np 
from Bio import SeqIO
import argparse
import matplotlib.pyplot as plt
from fasta import readFASTA

class ReadIn():
    def __init__(self,infile):
        self.infile = infile
        self.read_seq()
    
    def read_seq(self):
        # if there were more than 2 sequences, use a dictionary
        #for i,record in enumerate(SeqIO.parse(self.infile,'fasta')):
        #    if i == 0:
        #        self.query_id = record.id
        #        self.query_seq = record.seq
        #    if i == 1:
        #        self.subject_id = record.id
        #        self.subject_seq = record.seq
        sequences = readFASTA(self.infile)
        self.query_id = sequences[0][0]
        self.query_seq = sequences[0][1]
        self.subject_id = sequences[1][0]
        self.subject_seq = sequences[1][1]

class DNA_alignment(ReadIn):
    # HOXD70 scoring matrix
    # gap penalty of 10
    def __init__(self,infile,outfile):
        super().__init__(infile) 
        self.outfile = outfile
        self.gap_penalty = -300
        self.scoring_mat('needleman-wunsch/HOXD70.txt')      
        self.create_f_mat()
        self.align()
        self.calculate_score()

    def scoring_mat(self,score_file):
        header = []
        row_index = []
        matrix = []
        with open(score_file) as f:
            for i,line in enumerate(f.readlines()):
                line = line.split()
                if i == 0:
                    header = line
                else:
                    row_index.append(line[0])
                    matrix.append(line[1:])
        self.nucs = header 
        self.score_mat = np.array(matrix,dtype=np.int32)
        
    def create_f_mat(self):
        # n_rows = query length
        # n_cols = subject length
        self.f_mat = np.zeros((len(self.query_seq)+1,len(self.subject_seq)+1),dtype=np.int32)
        # query is rows index
        # subject is column index
        self.f_mat[:,0] = [j*self.gap_penalty for j in range(self.f_mat.shape[0])]
        self.f_mat[0] = [j*self.gap_penalty for j in range(self.f_mat.shape[1])]

    
    def align(self):
        for i,q in enumerate(self.query_seq):
            for j,s in enumerate(self.subject_seq):
                hoxd_i = self.nucs.index(q)
                hoxd_j = self.nucs.index(s)
                hoxd_score = self.f_mat[i,j]+self.score_mat[hoxd_i,hoxd_j]
                h = self.f_mat[i+1,j] + self.gap_penalty
                v = self.f_mat[i,j+1] + self.gap_penalty
                self.f_mat[i+1,j+1] = max(hoxd_score,h,v)

    
    def calculate_score(self):
        #traceback_mat = np.zeros((self.f_mat.shape[0], self.f_mat.shape[1]), dtype=np.int32)
        movement_dict = {'d':[0,0],'q':[0,1],'s':[1,0]}
        trace = []
        i = self.f_mat.shape[0]
        j = self.f_mat.shape[1]
        movement = [i,j]
        sub_f_mat = self.f_mat[movement[0]-2:movement[0], movement[1]-2:movement[1]] # slide tensor around f_mat

        while movement[0]!=1 or movement[1]!=1: # iterate upwards through rows

            # decide to align or gap
            # previous score is in (1,1) index
            # rows are query seq
            # columns are subject seq 
            diag = sub_f_mat[0][0] # move diag
            q_gap = sub_f_mat[1][0] # move up 
            s_gap = sub_f_mat[0][1] # move left
            if q_gap >= s_gap:
                choice = q_gap
                move_inds = [0,1]
                slide = 'q'
                
            else:
                choice = s_gap 
                move_inds = [1,0]
                slide = 's'
    
            if diag>=choice: # align
                # move diag
                movement[0] = movement[0] - 1
                movement[1] = movement[1] - 1
                slide = 'd'
                i = i-1
                j = j-1
            else: # do the gap
                movement[0] = movement[0] - move_inds[0]
                movement[1] = movement[1] - move_inds[1]
                i = i-movement_dict[slide][0]
                j = j-movement_dict[slide][1]

            sub_f_mat = self.f_mat[movement[0]-2:movement[0], movement[1]-2:movement[1]]
            trace.append(slide)
        trace = np.array(trace)
   
        seq1_ind = 0
        seq2_ind = 0
        aligned_seq1 = []
        aligned_seq2 = []
        for m in trace[::-1]:
            if m=='s':
                aligned_seq2.append('-')
                if seq1_ind < len(self.query_seq)-1:
                    aligned_seq1.append(self.query_seq[seq1_ind])
                    seq1_ind += 1
            elif m=='q':
                aligned_seq1.append('-')
                if seq2_ind < len(self.subject_seq)-1:
                    aligned_seq2.append(self.subject_seq[seq2_ind])
                    seq2_ind += 1
            else:
                if seq2_ind < len(self.subject_seq)-1:
                    aligned_seq2.append(self.subject_seq[seq2_ind])
                    seq2_ind += 1
                if seq1_ind < len(self.query_seq)-1:
                    aligned_seq1.append(self.query_seq[seq1_ind])
                    seq1_ind += 1
        
        with open('{}.fasta'.format(self.outfile),'w') as f:
            f.write('>{} ALIGNED\n'.format(self.query_id))
            f.write(''.join(aligned_seq1))
            f.write('\n>{} ALIGNED\n'.format(self.subject_id))
            f.write(''.join(aligned_seq2))
        chunks = np.arange(0,len(trace),50)
        
        print('score',self.f_mat[self.f_mat.shape[0]-1,self.f_mat.shape[1]-1])
        print('#'*55)
        print('S1: {}'.format(self.query_id))
        print('S2: {}'.format(self.subject_id))
        print('#'*55)
        print('\n')
        for i in range(1,len(chunks)):
            qc = ''.join(aligned_seq1[chunks[i-1]:chunks[i]])
            sc = ''.join(aligned_seq2[chunks[i-1]:chunks[i]])
            tc = trace[::-1][chunks[i-1]:chunks[i]]
            #if i <= len(trace):
            #    qc = ''.join(aligned_seq1[chunks[i-1]:chunks[i]])
            #    sc = ''.join(aligned_seq2[chunks[i-1]:chunks[i]])
            #    tc = trace[chunks[i-1]:chunks[i]]
            #else:
            #    qc = ''.join(aligned_seq1[chunks[i-1]:])
            #    sc = ''.join(aligned_seq2[chunks[i-1]:])
            #    tc = trace[chunks[i-1]:chunks[i-1]:]

            star = np.array((np.array(list(qc))==np.array(list(sc)))*1,dtype=object)
            lines = np.array((np.array(list(qc))==np.array(list(sc)))*1,dtype=object)
            gaps = np.where(tc=='s')[0].tolist()+np.where(tc=='q')[0].tolist()
       
            star[np.where(star==1)[0]] = ' ' # match
            star[np.where(star==0)[0]] = '*' # mismatch
            star[gaps] = ' ' # gap
            lines[np.where(lines==1)[0]] = '|'
            lines[np.where(lines==0)[0]] = '|'
            lines[gaps] = ' '
            print("  ",''.join(star))
            print('S1',qc,chunks[i])
            print('  ', ''.join(lines))
            print('S2' , sc)
            print('\n')

                
        '''
        (aligning -> gap in sequence 1 -> gap in sequence 2).
        (align --> q_gap --> s_gap)
        '''


class AA_alignment(DNA_alignment):
    # BLOSUM62 scoring matrix
    # gap penalty of 10
    def __init__(self,infile,outfile):
        self.infile = infile
        self.outfile = outfile
        self.read_seq()
        self.gap_penalty = -10
        self.scoring_mat('needleman-wunsch/BLOSUM62.txt')   
        self.create_f_mat()
        self.align()
        self.calculate_score()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-q','--query',required=True,help='Query sequence')
    parser.add_argument('-t','--type',required=True,help='Sequencing Type,(N=Nucleotide P=Protein)') # can accept N or P
    parser.add_argument('-o','--outfile',required=True,help='outfile name')
    args = parser.parse_args()
    if args.type=='N':
        D = DNA_alignment(args.query,args.outfile)
    elif args.type=='P':
        P = AA_alignment(args.query,args.outfile)
     