## SCRIPT FOR USING BLAST THROUGH BIOPYTHON

### Dependencies

> $ pip3 install biopython 

or you can download windows package from http://biopython.org/wiki/Download

## HOWTO

1. Paste at **organisms_<name_of_project>.txt** file "taxonomy name" (search query) 
2. At **main.py** find variable **MAM_GIDS** where add yours protein **Name** and **GID**
3. Run script:

> $ python3 main.py project=<name_of_project>

> $ python3 main.py project=gastro (example)


**<name_of_project>** min length 3 chars

Script will save alignments at directory saved_data_**<name_of_project>**
Script will save results at directory **results** with name **<name_of_project>.tsv**

You can open .tsv file at MS EXCEL or another program for edit received data