from itertools import groupby
from string import Template
import pandas as pd
import Excel_Writer

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


def file_parser_maker(search_seq, subject_file,SSR_list,Gene_list,Length_list):
    #opens the fasta file
    fiter = fasta_iter(subject_file)
    
    #makes it so the linebreak is maintained inside the excel sheet
    s = Template("\">$name\n$gene\"" + "\n")
    a = Template('>$name\n$gene')
    
    print('parsing for {}'.format(search_seq))
    #Creates a variable for the ID and sequence for each of the genes
    for ff in fiter:
        headerStr, seq = ff

        #searches for the SSR in the current gene
        if search_seq in seq:
            
            #only writes if the gene length is larger than 500
            if len(seq) > 500:
                if 'N' in seq:
                    print('Match found, but contains indeterminate bases')
                else:
                    #saves SSR only if it has a match with a large enough gene
                    SSR_list.append(search_seq)
                    #saves gene info
                    Gene_list.append('{}'.format(a.substitute(name=headerStr, gene=seq)))
                    Length_list.append('{}'.format(len(seq)))
                    print('Valid match found')
            else:
                print('Match found, but of insufficient length')
    return SSR_list,Gene_list


def main(SSR_file,subject_file,excel_file_name):
    
    print('Initiating search')
    
    #converts SSR file into a list
    search_list = input_file(SSR_file)
    #initializes lists
    SSR_list = []
    Gene_list = []
    Length_list = []
    #iterates through the SSRs and logs every time it appears in the subject file
    for i in range(len(search_list)):
        file_parser_maker(str(search_list[i]), subject_file,SSR_list,Gene_list, Length_list)
    
    print('Writing into excel file')
    
    #saves lists into pandas dataframes
    df1 = pd.DataFrame({'SSR': SSR_list})
    df2 = pd.DataFrame({'Gene': Gene_list})
    df3 = pd.DataFrame({'Gene Length': Length_list})
    
    #logs the dataframes into Excel
    Excel_Writer.append_df_to_excel(excel_file_name, df1,startcol = 0)
    Excel_Writer.append_df_to_excel(excel_file_name, df2,startcol = 1)
    Excel_Writer.append_df_to_excel(excel_file_name, df3,startcol = 2)

    
    print('done searching')

