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
    def __init__(self,infile):
        super().__init__(infile) 
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
        self.f_mat = np.zeros((len(self.query_seq),len(self.subject_seq)),dtype=np.int32)
        # query is rows index
        # subject is column index
        self.f_mat[:,0] = [j*self.gap_penalty for j in range(self.f_mat.shape[0])]
        self.f_mat[0] = [j*self.gap_penalty for j in range(self.f_mat.shape[1])]

    
    def align(self):
        for i,q in enumerate(self.query_seq):
            if i == 0:
                continue
            for j,s in enumerate(self.subject_seq):
                if j == 0:
                    continue
                hoxd_i = self.nucs.index(q)
                hoxd_j = self.nucs.index(s)
                hoxd_score = +self.f_mat[i-1,j-1]+self.score_mat[hoxd_i,hoxd_j]
                h = self.f_mat[i,j-1] + self.gap_penalty
                v = self.f_mat[i-1,j] + self.gap_penalty
                self.f_mat[i,j] = max(hoxd_score,h,v)
    
    def calculate_score(self):
        traceback_mat = np.zeros((self.f_mat.shape[0], self.f_mat.shape[1]), dtype=np.int32)
        score_list = []
        movement_dict = {'d':[0,0], 's':[0,1],'q':[1,0]}
        trace = []
        i = self.f_mat.shape[0]
        j = self.f_mat.shape[1]
        movement = [i,j]
        sub_f_mat = self.f_mat[movement[0]-2:movement[0], movement[1]-2:movement[1]] # slide tensor around f_mat

        while i>=1 and j>=1: # iterate upwards through rows
            if movement[0]==0 or movement[1]==0:
                break 

            else:
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
                score_list.append(sub_f_mat[1,1])
                sub_f_mat = self.f_mat[movement[0]-2:movement[0], movement[1]-2:movement[1]]
                trace.append(slide)
                
            
        print('score',self.f_mat[self.f_mat.shape[0]-1,self.f_mat.shape[1]-1])
        seq1_ind = 0
        seq2_ind = 0
        aligned_seq1 = []
        aligned_seq2 = []
        for m in trace[::-1]:
            if m=='s':
                aligned_seq2.append('-')
                aligned_seq1.append(self.query_seq[seq1_ind])
                seq1_ind += 1
            if m=='q':
                aligned_seq1.append('-')
                aligned_seq2.append(self.subject_seq[seq2_ind])
                seq2_ind += 1
            else:
                aligned_seq1.append(self.query_seq[seq1_ind])
                aligned_seq2.append(self.query_seq[seq2_ind])
                seq1_ind+=1
                seq2_ind+=1
        with open('./aligned.fasta','w') as f:
            f.write('>{} ALIGNED\n'.format(self.query_id))
            f.write(''.join(aligned_seq1))
            f.write('\n>{} ALIGNED\n'.format(self.subject_id))
            f.write(''.join(aligned_seq2))
                
        '''
        (aligning -> gap in sequence 1 -> gap in sequence 2).
        (align --> q_gap --> s_gap)
        '''


class AA_alignment(DNA_alignment):
    # BLOSUM62 scoring matrix
    # gap penalty of 10
    def __init__(self,infile):
        self.infile = infile
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
    args = parser.parse_args()
    if args.type=='N':
        D = DNA_alignment(args.query)
    elif args.type=='P':
        P = AA_alignment(args.query)
     