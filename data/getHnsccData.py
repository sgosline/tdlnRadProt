#!/usr/local/bin/python
'''
Basic CLI to import CPTAC proteomic data
'''

import cptac
import argparse


def main():
    cptac.download(dataset='hnscc')

    dat = cptac.Hnscc()
    df = dat.get_proteomics()
    df.to_csv(path_or_buf='file.tsv', sep="\t")

if __name__=='__main__':
    main()
