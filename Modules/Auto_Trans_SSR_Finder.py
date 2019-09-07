# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 08:49:16 2018

Uses BLASTx to translate and check if an SSR is contained within the longest ORF of the Protein

@author: ark47
"""
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from itertools import groupby
import pandas as pd
import Excel_Writer

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


def fasta_list_maker(input_file_name):

    fiter = fasta_iter(input_file_name)
    fasta_list = []
    for ff in fiter:
        headerStr, seq = ff
        fasta_list.append('>{}\n{}'.format(headerStr,seq))
    print(fasta_list)    
    return fasta_list

def main(excel_file_name):
    Translated_SSR_list = []
    SSR_fasta_list = fasta_list_maker('AutoJobber_Logs/SSR_fasta.fa')
    subject_fasta_list = fasta_list_maker('AutoJobber_Logs/Gene Translations.fa')
    save_file_name ='AutoJobber_Logs/test_blaster.xml'
    save_file = open(save_file_name, 'w+')   
    num_line = 0
    
    for i in range(len(subject_fasta_list)):
        query_file = open('AutoJobber_Logs/BLAST query.fa','w+')
        query_file.write(SSR_fasta_list[i])
        print(SSR_fasta_list[i])
        query_file.close()
        
        subject_file = open('AutoJobber_Logs/BLAST subject.fa', 'w+')
        subject_file.write(subject_fasta_list[i])
        print(subject_fasta_list[i])
        subject_file.close()
        
        blastx_cline = NcbiblastxCommandline(cmd='blastx', out= save_file_name, outfmt=5, query='AutoJobber_Logs/BLAST query.fa', 
        subject ='AutoJobber_Logs/BLAST subject.fa')
        stdout, stderr = blastx_cline()
        
        save_file = open(save_file_name, 'r')
        num_line += 1
        subject_seq = 'No Match'
        record = NCBIXML.read(save_file)
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print(alignment)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)
                subject_seq = hsp.sbjct
                break
            break
        Translated_SSR_list.append(subject_seq)
    save_file.close()
    print(Translated_SSR_list)
    if len(Translated_SSR_list) == num_line:
        print('All SSRs accounted for')
    
    df1 = pd.DataFrame({'Translated SSR': Translated_SSR_list})
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 6)
