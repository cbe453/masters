#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import numpy as np
from scipy import stats
from collections import Counter

def parseData(cdsPath):
    cdsArray = []
    for cds in SeqIO.parse(cdsPath, "fasta"):
        cdsArray.append(len(cds.seq))

    return(cdsArray)

def main(args):
    data1 = parseData(args.input1)
    data2 = parseData(args.input2)
    sortedLogData1 = sorted(np.log10(data1))
    sortedLogData2 = sorted(np.log10(data2))
 
    results = stats.ks_2samp(sortedLogData1, sortedLogData2)
    print(results)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Kolmogorov-Smirnov Test between CDFs of gene lengths')
    parser.add_argument('-1', type=str, dest='input1', help='Path to CDS file from gene prediction tool 1.')
    parser.add_argument('-2', type=str, dest='input2', help='Path to CDS file from gene prediction tool 2.')
    args = parser.parse_args()


    exit(main(args))
