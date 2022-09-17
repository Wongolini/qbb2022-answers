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

class ReadIn():
    def __init__(self,infile):
        self.infile = infile
        self.read_seq()
    
    def read_seq(self):
        # if there were more than 2 sequences, use a dictionary
        for i,record in enumerate(SeqIO.parse(self.infile,'fasta')):
            if i == 0:
                self.query_id = record.id
                self.query_seq = record.seq
            if i == 1:
                self.subject_id = record.id
                self.subject_seq = record.seq


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
        score_path = np.zeros((self.f_mat.shape[0]),dtype=np.int32)
        for i in range(self.f_mat.shape[0])[::-1]: # iterate upwards through rows
            row = self.f_mat[i]
            if i == self.f_mat.shape[0]-1:
                # next choice can be vertical or direct diagonal
                best_index = np.argmax(row) # index of the max 
                score_path[i] = best_index
            else:
                #best_index = min(abs(np.argsort(row)[::-1]-best_index)) # look for the smallest difference... should be 0 or 1
                prev_best_ind = score_path[i+1]
                if row[prev_best_ind-1] > row[prev_best_ind]:
                    score_path[i] = prev_best_ind-1
                elif row[prev_best_ind-1] <= row[prev_best_ind]:
                    score_path[i] = prev_best_ind
    
        self.score_path = score_path[::-1]
        self.score = np.sum([self.f_mat[i][score] for i,score in enumerate(score_path)])
        print('SCORE: ',self.score)
        plt.plot(self.score_path,range(self.score_path.shape[0]))
        #plt.xticks(list(range(len(self.subject_seq)))[::10],list(self.subject_seq)[::10],fontsize=5)
        #plt.yticks(list(range(len(self.query_seq))[::20]),list(self.query_seq)[::20],fontsize=5)
        plt.xlabel('SUBJECT index')
        plt.ylabel('QUERY index')
        plt.title('Alignment Score Path')
        plt.show()


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
     