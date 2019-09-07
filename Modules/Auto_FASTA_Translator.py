# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 09:45:16 2018

Will translate nucleotide sequences into the longest possible open reading 
frame and record them into an excel file alongside length, direction, 

@author: ark47
"""

import pandas as pd
from Bio.Seq import Seq
from itertools import groupby
import Excel_Writer

def strand_dictionary(strand):
    #translates the working frame syntax into a readable one
    strand_translate =	{
            -1: "3'5",
            1: "5'3",
            }
    return strand_translate[strand]

def fasta_iter(fasta_name):
    
    fh = open(fasta_name)

    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

def dna_translator(temp_file_name,Translated_Protein,Protein_Length,Seq_Direction,Frame_Info,Complete_Translated_Protein):
    #opens the fasta file and converts it into arrays of strings grouped by their title and sequence
    fiter = fasta_iter(temp_file_name)
    
    #opens new file that the translations go into
    Trans_fasta = open('AutoJobber_Logs/Gene Translations.fa', 'w+')
    
    #groups the gene title and sequence into the same loop
    for ff in fiter:
        headerStr, seq = ff
        
        #turns the sequence into a Biopython Seq object
        seq = Seq(str(seq))

        #sets table to a standard codon one
        table = 1
        longest_trans = ""
        longest_frame = 0
        longest_strand = 1
        
        #iterates through the 1-3 frames of the forwards and reverse complement stands
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                #cuts the translation off at the stop codon
                for pro in nuc[frame:].translate(table).split("*"):
                    output = []
                    process = False
                   
                    #starts reading at the first M and will read until the end of the translation
                    for i in range(len(str(pro))):
                        if (list(str(pro))[i] == 'M'):
                            process = True
                        if process:
                            #no need to worry about stop codons since they were accounted for
                            output.append(str(pro)[i])
                    
                    #converts the list into a string        
                    trans_seq = "".join(str(x) for x in output)       
                    
                    #only records the longest protein
                    if len(trans_seq) > len(longest_trans):
                        longest_trans = trans_seq
                        longest_frame = frame
                        longest_strand = strand
                        
        #records the complete translation, as the former block only keeps the ORF
        if longest_strand == 1:             
            longest_complete_trans = seq[longest_frame:].translate(table)
        elif longest_strand == -1:
            rev_comp_seq = seq.reverse_complement()
            longest_complete_trans = rev_comp_seq[longest_frame:].translate(table)
        
        #the complete translation is saved in a fasta format with the STOP codons being represented with 'Z' as the aligner only takes letter inputs
        Trans_fasta.write(">{} Translation\n{}\n".format(headerStr,str(longest_complete_trans).replace('*', 'Z')))
        #records info about the protein into lists for later input
        Translated_Protein.append(longest_trans)
        Protein_Length.append(len(longest_trans))
        Seq_Direction.append(str(strand_dictionary(longest_strand)))
        Frame_Info.append(longest_frame+1)
        Complete_Translated_Protein.append(longest_complete_trans)
        print (">%s Translation\n%s\tlength %i, strand %s, frame %i\n" % (headerStr, longest_trans, len(longest_trans), strand_dictionary(longest_strand), longest_frame+1))
    return Translated_Protein,Frame_Info                        


def main(excel_file_name):
    #calls on the text file that has the gene data
    temp_file_name = 'AutoJobber_Logs/SSR_Containing_Genes.fa'
    temp_file = open(temp_file_name,'w+')
    df = pd.read_excel(excel_file_name, sheet_name='Sheet1')
    Translated_Protein = []
    Complete_Translated_Protein = []
    Seq_Direction = []
    Frame_Info = []
    Protein_Length = []
    
    #writes genes into a text file in fasta format        
    for i in df.index:
        temp_file.write('{}\n'.format(df['Gene'][i]))
    temp_file.close()  
    
    #finds the protein translation with the longest ORF and records data about it
    dna_translator(temp_file_name,Translated_Protein,Protein_Length,Seq_Direction,Frame_Info,Complete_Translated_Protein)
    
    print('Writing into excel file')
    
    #formats lists as pandas dataframes
    df1 = pd.DataFrame({'Translation': Translated_Protein})
    df2 = pd.DataFrame({'Protein Length': Protein_Length})
    df3 = pd.DataFrame({'Direction': Seq_Direction})
    df4 = pd.DataFrame({'Reading Frame': Frame_Info})
    df5 = pd.DataFrame({'Completed Translation': Complete_Translated_Protein})
 
    #uses Excel_Writer module to write into the specified file
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 4)
    Excel_Writer.append_df_to_excel(excel_file_name, df5,startcol = 5)
    Excel_Writer.append_df_to_excel(excel_file_name, df2,startcol = 8)
    Excel_Writer.append_df_to_excel(excel_file_name, df3,startcol = 9)
    Excel_Writer.append_df_to_excel(excel_file_name, df4,startcol = 10)
    print('done')


    