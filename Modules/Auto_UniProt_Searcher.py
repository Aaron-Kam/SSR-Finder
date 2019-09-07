# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 20:13:27 2018

@author: ark47
"""

from Bio import ExPASy
from Bio import SeqIO
import Excel_Writer
import pandas as pd

def input_file(data_file):
    text_file = open(data_file, "r")
    lines = text_file.readlines()
    lines = [x.strip() for x in lines]
    text_file.close()
    return lines

def accession_finder(accessions,hit_id_list):
 
    for input_id in hit_id_list:
        bar_number = 0
        read_id = False
        accession = ''
        if input_id == '':
            accession = 'No Result'
        else:
            for letter in range(len(input_id)):
                if bar_number == 3:
                    read_id = True
                if input_id[letter] == '|' or input_id[letter] == '.':
                    bar_number += 1
                if bar_number == 4:
                    read_id = False
                if read_id == True:
                    accession = accession + input_id[letter]         
        
        print(accession)
        accessions.append(accession)
        
def uniprot_searcher(accessions,hit_id_list, Top_hit, Keyword_list,Organism_list):
    seq_num = 0
    for accession in accessions:
        if accession == 'No Result':
            Keyword_list.append('No Match')
            Organism_list.append('No Match')
            Top_hit.append('No Match')
        else:
            with ExPASy.get_sprot_raw(accession) as handle:
                seq_record = SeqIO.read(handle, "swiss")
            Top_hit.append('>{}\n{}\n'.format(hit_id_list[seq_num],seq_record.seq))
            Keyword_list.append(seq_record.annotations["keywords"])
            Organism_list.append(seq_record.annotations["organism"])
            seq_num += 1
            
            print(seq_record.id)
            print(seq_record.name)
            print(seq_record.description)
            print(repr(seq_record.seq))
            print("Length %i" % len(seq_record))
            print(seq_record.annotations["keywords"])
        
def main(excel_file_name):
    accessions = []
    Top_hit = []
    Keyword_list = []
    Organism_list = []
    
    hit_id_list = input_file("AutoJobber_Logs/SwissProt_Top_BLAST_Hits.txt")
    accession_finder(accessions,hit_id_list)
    uniprot_searcher(accessions, hit_id_list, Top_hit, Keyword_list,Organism_list)
    
    df1 = pd.DataFrame({'Top Hit Keywords (SwissProt database)': Keyword_list})
    df2 = pd.DataFrame({'Top Hit and Identities (SwissProt database)': Top_hit})
    df3 = pd.DataFrame({'Organism (SwissProt database)': Organism_list})
    
    Excel_Writer.append_df_to_excel(excel_file_name, df2,startcol = 15)
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 16)    
    Excel_Writer.append_df_to_excel(excel_file_name, df3,startcol = 17)

