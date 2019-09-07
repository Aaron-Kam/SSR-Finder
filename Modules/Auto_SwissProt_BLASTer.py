# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 10:40:30 2018

@author: ark47
"""

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO
from Bio.Blast import NCBIXML
import Excel_Writer
import pandas as pd


def blaster(query_file,temp_top_hits,save_file_name,Top_hit,Top_three):
    
    #opens the complete xml file for writting
    save_file = open(save_file_name, "w+")
    #keeps track of the sequence number
    seq_num = 1
    ERROR_COUNT = 0
    ERROR_LOCATIONS = []
    #BLASTs all sequences against the SwissProt Database
    for seq_record in SeqIO.parse(query_file, "fasta"):
        #opens the xml file for single searches
        temp_xml_file_name = "AutoJobber_Logs/temp_swissprot_results.xml"
        temp_xml_file = open(temp_xml_file_name, "w+")
        #writes the query into a file for blasting
        blast_query_file = open('AutoJobber_Logs/swissprot_query.fa', 'w+')
        blast_query_file.write(str(seq_record.seq))
        blast_query_file.close()
        
        #prints the blasted sequence's information
        print('BLASTing Sequence {}...'.format(seq_num))
        print(seq_record.id)
        print(repr(seq_record.seq))
        print('Length:{}\n'.format(len(seq_record)))
        
        #actual blast commandline
        result_handle = NcbiblastxCommandline(cmd='blastx', db = "swissprot/swissprot", 
                                     outfmt=5, out= temp_xml_file_name, query = 'AutoJobber_Logs/swissprot_query.fa')
        #iterates the sequence number
        seq_num += 1
        stdout, stderr = result_handle()
        temp_xml_file.close()
        #opens the single blast result file
        temp_xml_file = open(temp_xml_file_name, "r")
        BLAST_xml_string = temp_xml_file.read()
        #saves the single result into a collective xml file
        save_file.write(BLAST_xml_string)
        if BLAST_xml_string == "":
            print('THIS IS BLANK')
            ERROR_COUNT += 1
            ERROR_LOCATIONS.append(seq_num)
        temp_xml_file.close()
        
    save_file.close()
    
    save_file = open(save_file_name,'r')
    
    #reads through the xml file
    records = NCBIXML.parse(save_file)
    
    #parses through the records in the file
    for i, record in enumerate(records):
        hit_rank = 0
        #makes a line break after each line, except the first time through
        if i > 0:
            temp_top_hits.write('\n')
        #Sets parameter to true in case there are no top 3 hits
        Top_three_hit_bool = True
        
        
        temp_results = []
        print(i)
        #adds a slot for the alignment in case there are no matches
        Top_hit.append('No Match')
        Top_three.append('No Match')
        #extracts info about the BLAST results
        for align in record.alignments:
            for hsp in align.hsps:
                hit_rank += 1   
                #saves specific info about the first hit
                if hit_rank == 1:
                    print(align.hit_id)
                    Top_hit[i] = '{}'.format(align.hit_id)
                    temp_top_hits.write('{}'.format(align.hit_id))
                    
                print('Hit no. {}'.format(hit_rank))
                print('{} Length: {}, e-value: {}'.format(align.hit_id, align.length,hsp.expect))
                temp_results.append('Hit no.{}\n{}'.format(hit_rank,align.hit_id))
                #only iterates through the 
                break
            
            #saves the top 3 hits into a single element in a list
            if hit_rank == 3:
                triplet = ''
                for j in range(len(temp_results)):
                    triplet = triplet + str(temp_results[j]) + '\n'
                
                Top_three[i] = triplet
                Top_three_hit_bool = False
                break
        if Top_three_hit_bool:
            Top_three[i] = 'Less than three matches'
            
        print('')
        
    temp_top_hits.close()    
    save_file.close()
    return Top_hit,Top_three
         
def main(excel_file_name):
    print('Initating SwissProt BLAST')
    temp_xml_results = "AutoJobber_Logs/SwissProt_BLAST_Results.xml"
    top_hits_fasta = "AutoJobber_Logs/SwissProt_Top_BLAST_Hits.txt"
    Top_hit = []
    Top_three = []
    #erases the file of the same length creates file if one of the same name doesnt exist
    temp_top_hits = open(top_hits_fasta, "w+")
    
    #BLASTs sequences
    blaster("AutoJobber_Logs/SSR_Containing_Genes.fa",temp_top_hits,temp_xml_results,Top_hit,Top_three)
    print('done BLASTing')
    
    #writes data into the excel file
    df1 = pd.DataFrame({'Top Three Hits (SwissProt database)': Top_three})
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 14)

