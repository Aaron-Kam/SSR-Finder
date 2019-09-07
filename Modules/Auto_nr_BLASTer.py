# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 21:31:51 2018

@author: Aaron Kam

Will BLAST a series of sequences 
"""

from Bio.Blast import NCBIWWW
from itertools import groupby
from Bio import SeqIO
from Bio.Blast import NCBIXML
import Excel_Writer
import pandas as pd


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

def blaster(query_file,temp_top_hits,save_file_name,Top_hit,Top_three):
    save_file = open(save_file_name, "w+")
    seq_num = 1
    
    #BLASTs all Genes against nr protein database and saves it in a xml file
    for seq_record in SeqIO.parse(query_file, "fasta"):
        print('BLASTing Sequence {}...'.format(seq_num))
        print(seq_record.id)
        print(repr(seq_record.seq))
        print('Length:{}\n'.format(len(seq_record)))
        result_handle = NCBIWWW.qblast("blastx","nr",seq_record.seq)
        seq_num = seq_num + 1
        save_file.write(result_handle.read())
        
    save_file.close()
    save_file = open(save_file_name,'r')
    records = NCBIXML.parse(save_file)
    
    #exracts information for each of the BLAST searches
    for i, record in enumerate(records):
        hit_rank = 0
        print("Input {}".format(i+1))
        temp_results = []
        for align in record.alignments:
            for hsp in align.hsps:
                hit_rank = hit_rank +1
                #Records the top hit
                if hit_rank == 1:
                    Top_hit.append('>{}\n'.format(align.title))
                    temp_top_hits.write('>{}\n'.format(align.title))
                print('Hit no. {}'.format(hit_rank))
                print('{} Length: {}, e-value: {}'.format(align.hit_id, align.length,hsp.expect))
                temp_results.append('Hit no.{}\n{}'.format(hit_rank,align.title))
                break
            #Records Top three hits as single element in a list
            if hit_rank == 3:
                triplet = ''
                for j in range(len(temp_results)):
                    triplet = triplet + str(temp_results[j]) + '\n'
                Top_three.append(triplet)
                break
        print('')
    temp_top_hits.close()    
    save_file.close()

    result_handle.close()
    return Top_hit,Top_three
         
def main(excel_file_name):
    print('Initating nr BLAST')
    temp_xml_results = "AutoJobber_Logs/BLAST_Results.xml"
    top_hits_fasta = "AutoJobber_Logs/Top_BLAST_Hits.fa"
    Top_hit = []
    Top_three = []
    
    temp_top_hits = open(top_hits_fasta, "w+")
    
    #calls on blaster function
    blaster("AutoJobber_Logs/SSR_Containing_Genes.fa",temp_top_hits,temp_xml_results,Top_hit,Top_three)
    print('done BLASTing')
    
    #Records info in Excel via my Excel_Writer module and pandas
    df1 = pd.DataFrame({'Top Three Hits': Top_three})
    df2 = pd.DataFrame({'Top Hit and Identities': Top_hit})
    
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 12)
    Excel_Writer.append_df_to_excel(excel_file_name, df2,startcol = 13)
    
    