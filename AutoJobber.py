# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 19:48:42 2018

@author: ark47

Takes a transcriptome and an SSR list and
analyzes their relationship
"""

from Runtime_Timer import Timer
from Modules import Auto_RegEx
from Modules import Auto_FASTA_Translator
from Modules import Auto_nr_BLASTer
from Modules import Auto_Seq_Aligner
from Modules import Auto_SwissProt_BLASTer
from Modules import Auto_UniProt_Searcher
from Modules import Auto_Trans_SSR_Finder
from Modules import Auto_Trans_Seq_Aligner
from Modules import Auto_FASTA_Converter
from Modules import Auto_Hydrophobicity_Plotter


Timer.Start()
#name of the transcriptome
transcriptome = 'akoa.fa'
#name of the SSR files seperated by comma
#e.g. ['10(2)_SSR.txt', '20(2)_repeats.txt']
input_list = ['cleaned_test_SSR.txt']

input_file_path = 'AutoJobber_Logs/SSR Files/'
output_file_path = 'AutoJobber_Logs/Output Files/test/'

for SSR_file_name in input_list:
    #Current code assumes SSR files follow the above format
    excel_file_name ='{} Compiled Data.xlsx'.format(str(SSR_file_name)[:5])
    
    #Combines the file path and name into a single variable
    SSR_file = input_file_path + SSR_file_name
    excel_file = output_file_path + excel_file_name
    
    #Finds SSR within fasta file and records the gene and title into an excel sheet
    Auto_RegEx.main(SSR_file,transcriptome,excel_file)
    
    #Converts SSR list into a fasta file for later use
    Auto_FASTA_Converter.main(excel_file)
    
    #Translates Genes from Regex into Amino acid sequence with the longest ORF, and records it inside excel
    Auto_FASTA_Translator.main(excel_file)
    
    #Uses BLASTx to translate the SSR and tests if it is within the Translated Gene
    Auto_Trans_SSR_Finder.main(excel_file)
    
    #Finds location of SSR inside gene and logs it in excel
    Auto_Seq_Aligner.main(excel_file)
    
    #Finds location of translated SSR within Gene and logs it in excel
    #note: it will not consider copies of the translated SSR within the gene
    Auto_Trans_Seq_Aligner.main(excel_file)
    
    #BLASTs SSR containing genes against the nr and SwissProt databases
    #Auto_nr_BLASTer.main(excel_file)
    Auto_SwissProt_BLASTer.main(excel_file)
    
    #Takes SwissProt accession number from BLAST and runs it against the UniProt.org site
    #Records keywords associated with the entry, including domain and function
    Auto_UniProt_Searcher.main(excel_file)
    
    #Plots the Hydrophobicity of the existing translated SSRs on the 
    Auto_Hydrophobicity_Plotter.main(excel_file,str(SSR_file_name)[:5],output_file_path)
    
    
Timer.End()

print('done')