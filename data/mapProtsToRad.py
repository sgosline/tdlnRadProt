'''
Downloads radiomic data
maps patient identifiers to PDC
gets protein data
'''


from osfclient import cli
import pandas as pd
import os

##adding this here to prep for the download of data
token = "C6a7xYr87esl61NiQaXOHCYZykMVwiE0kLscQkD62fXmfMVgj1A1fFgJQkDG7ZuyNWeEc5"
image_dat = '5x4n2/NewImageData02142017.xlsx'
prot_dat = '5x4n2/S054_HNSCC_imputed_0920.tsv'

loc_image_dat = 'NewImageData02142017.xlsx'
loc_prot_dat = 'S054_HNSCC_imputed_0920.tsv'

def main():
    ##download the files
    ##unclear how to do this, cheating right now
    image = pd.read_excel(loc_image_dat)
    prot = pd.read_tsv(loc_prot_dat)
    ##get the file mapping


if __name__=='__main__':
    main()
