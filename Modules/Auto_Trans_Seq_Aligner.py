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

    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

def aligner(new_file,SSR_list):
    
    #converts a fasta file into a tuple for faster reading
    fiter = fasta_iter('AutoJobber_Logs/Gene Translations.fa')
    
    SSR_index = 0
    #parses through each entry in the fasta file
    for ff in fiter:
        headerStr, seq = ff
        #if the SSR does not appear in the sequence, it will make a blank space where the alignment would be
        if SSR_list[SSR_index] == 'No Match':
            new_file.write('No Match\n\n\n\n')
        #aligns the sequences and writes the result to a text file
        else:        
            alignments = pairwise2.align.globalds(SSR_list[SSR_index], seq, blosum62,-10,-0.5, penalize_end_gaps = False, one_alignment_only = True)
            for a in alignments:
                new_file.write(pairwise2.format_alignment(*a))
        SSR_index += 1
    new_file.close()
    
def finder(SSR_list,SSR_locations):
    #writes the alignment file into a list
    lines = input_file('AutoJobber_Logs/Trans Align products.txt')
    SSR_index = 0
    for match_line in range(len(lines)):
        #reads every fourth line as thats where the tranlsated SSR is written
        if match_line % 4 == 0:
            match = str(lines[match_line])
            SSR_letter = 0
            start_index = 0
            end_index = 0
            counting = False
            SSR = str(SSR_list[SSR_index])
            print(SSR)
            #skips over the lines with no matches, and logs it
            if match == 'No Match':
                SSR_locations.append('No Match')    
            else:
                #scans through the line until it finds the SSR and records the range
                for letter in range(len(match)):                    
                    if match[letter] == SSR[SSR_letter]:
                        if SSR_letter == len(SSR)-1:
                           counting = False
                           end_index = letter + 1
                           print('end')
                           print('Start: {}, End {}'.format(start_index,end_index))
                           SSR_locations.append('{}-{}\n'.format(start_index,end_index))
                        if SSR_letter == 0:
                            counting = True
                            start_index = letter + 1
                            print('start')
                        if counting == True:    
                            SSR_letter += 1
            #increments the SSR index to keep the SSR and genes in sync                
            SSR_index += 1      
            
           
   
    
def main(excel_file_name):
    #reads the column in the excel sheet that contains all of the gene translations
    df = pd.read_excel(excel_file_name, sheet_name='Sheet1')
    SSR_list = df['Translated SSR']
    #creates temporary log to hold the align products
    new_file = open('AutoJobber_Logs/Trans Align products.txt','w+')
    SSR_locations = []
    print(SSR_list)
    #aligns SSR with translated sequence and writes it to a designated text file
    aligner(new_file, SSR_list)
    #if translated SSR exists then it will find its location in the sequence
    finder(SSR_list,SSR_locations)
    df1 = pd.DataFrame({'Translated SSR Location': SSR_locations})
    
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 7)
    new_file.close()
