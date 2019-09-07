# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 19:41:56 2018

@author: ark47
"""

from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import pandas as pd
import Excel_Writer
from itertools import groupby


def input_file(data_file):
    text_file = open(data_file, "r")
    lines = text_file.readlines()
    lines = [x.strip() for x in lines]
    text_file.close()
    return lines

def fasta_iter(fasta_name):
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

def aligner(new_file,SSR_list):
    #converts fasta file into a set of tuples for speed
    fiter = fasta_iter('AutoJobber_Logs/SSR_Containing_Genes.fa')
    SSR = 0
    #iterates through all of the sequences
    for ff in fiter:
        headerStr, seq = ff
        #uses pairwise2 align tool to generate alignments
        alignments = pairwise2.align.globalds(SSR_list[SSR], seq, blosum62,-10,-0.5, penalize_end_gaps = False, one_alignment_only = True)
        #writes alignments into a text document
        for a in alignments:
            new_file.write(pairwise2.format_alignment(*a))
        #increments the SSR index
        SSR += 1
    new_file.close()

def finder(SSR_list,SSR_locations):
    #converts the alignment file into a list
    lines = input_file('AutoJobber_Logs/Align products.txt')
    SSR_index = 0
    #iterates through every four lines in the alignment file, as those are the ones that contain the SSR
    for match_line in range(0, len(lines), 4):
        #assigns every line of the match data to a variable 
        match = str(lines[match_line])
        #declares the indexes of the lists
        letter_index = 0
        start_index = 0
        end_index = 0
        counting = False
        SSR = str(SSR_list[SSR_index])
        print(SSR)
        #iterates through the letters in the match
        for letter in range(len(match)):
            if match[letter] == SSR[letter_index]:
                if letter_index == len(SSR)-1:
                   counting = False
                   end_index = letter + 1
                   print('end')
                   print('Start: {}, End {}'.format(start_index,end_index))
                   SSR_locations.append('{}-{}\n'.format(start_index,end_index))
                if letter_index == 0:
                    counting = True
                    start_index = letter + 1
                    print('start')
                if counting == True:    
                    letter_index += 1
        SSR_index += 1      
            
           
   
    
def main(excel_file_name):
    #reads the SSRs from the excel file into a list
    df = pd.read_excel(excel_file_name, sheet_name='Sheet1')
    SSR_list = df['SSR']
    #opens text log file
    new_file = open('AutoJobber_Logs/Align products.txt','w+')
    SSR_locations = []
    print(SSR_list)
    
    #aligns SSRs against their genes and produces an alignment table
    aligner(new_file, SSR_list)
    #reads the alignment table and takes note of their 
    finder(SSR_list,SSR_locations)
    
    #converts the location data into a dataframe
    df1 = pd.DataFrame({'SSR Location': SSR_locations})
    #writes data into excel
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 3)
    new_file.close()
    
