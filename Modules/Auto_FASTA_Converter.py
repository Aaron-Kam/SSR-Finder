# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 09:57:01 2018

@author: ark47
"""

import pandas as pd

def main(excel_file):
    #converts SSR file into a list
    
    SSR_fasta = open('AutoJobber_Logs/SSR_fasta.fa','w+')
    index = 0
    df = pd.read_excel(excel_file, sheet_name='Sheet1')
    SSR_list = df['SSR']
    #writes SSRs into a fasta file
    for SSR in SSR_list:
        index += 1
        print(SSR)
        SSR_fasta.write('>SSR no {}\n{}\n'.format(index,SSR))
    SSR_fasta.close()

