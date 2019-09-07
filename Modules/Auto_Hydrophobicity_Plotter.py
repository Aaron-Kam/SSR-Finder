# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:02:41 2018

@author: ark47
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Excel_Writer


def Hydrophobicity_dict(SSRs, SSR_Hydro_Values):
   
   #Kyte and Doolittle hydrophobicity dictionary
   kd = {"A": 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 
         'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 
          'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2, '_': "STOP"} 
   
   #iterates through SSRs and only reads ones that exists
   for SSR in SSRs:
       ind_list = []
       if SSR == 'No Match':
           SSR_Hydro_Values.append('')
       else:
           #creates 2D arraylist that contains the hydrophobicity values for each SSR
           for letter in range(len(str(SSR))):
               if str(SSR)[letter] != 'Z':   
                   ind_list.append(kd[str(SSR)[letter]])
           
           SSR_Hydro_Values.append(ind_list) 


def main(excel_file_name,SSR_type,output_file_path):
    #reads the translated SSRs from Excel
    df = pd.read_excel(excel_file_name, sheet_name='Sheet1')
    Trans_SSRs = df['Translated SSR']
    
    #reads raw SSRs for naming's sake
    SSRs = df['SSR']
    SSR_Hydro_Values = []
    
    #converts proteins into their hydrophobicity values
    Hydrophobicity_dict(Trans_SSRs,SSR_Hydro_Values)
    
    index_count = 0
    
    #iterates through the individual SSR's hydrophobicity values and skips over nonexistent ones
    for SSR in SSR_Hydro_Values:
        if SSR != '':
            #SSR means the translated Hydrophobicity values
            y = SSR
            #x is in increments of 1, starting 
            x = np.arange(1,len(y) + 1, 1)
            
            #creates graph
            fig = plt.figure()
            hyd = fig.add_subplot(111)
            
            #creates file title
            SSR_length = int(len(str(SSRs[index_count]))/2)
            title = str(SSRs[index_count])[:SSR_length]
            #assigns attributes to the graph
            plt.title("Hydrophobicity of {}".format(title))
            plt.grid(True)
            plt.plot(x,y)
            plt.xlabel('Position', fontsize=12)
            plt.ylabel('Value', fontsize=12)
            #annotates the graph with the x,y coordinates
            for xy in zip(x,y):
                hyd.annotate('(%s, %s)' % xy, xy=xy, textcoords='data')
            fig.savefig('{}{} {}.jpeg'.format(output_file_path,SSR_type, title))
            fig.clf()
        index_count += 1
    
    #writes data into excel file    
    df1 = pd.DataFrame({'Hydrophobicity Values': SSR_Hydro_Values})    
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 11)   
