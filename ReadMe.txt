Auto_SSR_Analysis v2.1

*********************************************
WARNING: CODE IS UNPOLISHED AND NOT OPTIMIZED
USE SOMETHING ELSE FIRST
*********************************************

Instructions:
    
    Dependencies:
        -Biopython
            *to download type 'pip install biopython' in the command prompt without quotes
            * when python is installed
        -ncbi-blast-2.8.1+ suite
            *can be downloaded here ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
        -pandas
            *should be automatically downloaded with python
        -swissprot prot db
            *can be downloaded here ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz
            *place the swissprot db inside a folder called swissprot, and place that folder inside the modules directory
        
    To run in commandline:
        1)install python
            *can be downloaded here https://www.python.org/downloads/
        2)navigate to the Auto_SSR_Analysis folder where AutoJobber.py is
        3)type "python AutoJobber.py" into the commandline without
          quotations
    
    To configure the program:
        
        File placement:
            -place transcriptome file inside the same directory as AutoJobber.py
            -place SSR text files in the directory SSR Files
                *SSRs should be delimmited by line
            -the graphs and spreadsheets should be inside the Output Files Directory
        
        Changing files to use
            -Open AutoJobber.py in your favorite text editor
            -Read the comments and follow their directions
            -To omit the use of a module, place a # in front of it
                *to use it again, delete the #
        
        Important Caveats:
            -The Uniprot Searcher module is reliant on the output of the Swissprot Blast module. Only
                run the Uniprot Searcher if Swissprot has been run successfully
            -The Swissprot Blast module requires an installation of the swissprot db which can be found on
                the NCBI FTP server. Place the swissprot protein db inside a folder called swissprot, 
                inside the modules folder.
                
