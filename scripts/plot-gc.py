#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main(args):
    print(args)
    dc1 = pd.read_csv(r"X:/Connor/masters/results/gc-content/dc1/dc1-gc-raw.tsv", header=None, sep='\t', encoding='utf-8')
    tsth20 = pd.read_csv(r"X:\Connor\masters\results\gc-content\tsth20\tsth20-gc-raw.tsv", header=None, sep='\t', encoding='utf-8')
    reesei = pd.read_csv(r"X:\Connor\masters\results\gc-content\t-reesei\t-reesei-gc-raw.tsv", header=None, sep='\t', encoding='utf-8')
    harzianum = pd.read_csv(r"X:\Connor\masters\results\gc-content\t-harzianum\t-harzianum-gc-raw.tsv", header=None, sep='\t', encoding='utf-8')
    virens = pd.read_csv(r"X:\Connor\masters\results\gc-content\t-virens\t-virens-gc-raw.tsv", header=None, sep='\t', encoding='utf-8')
    
    plt.plot(dc1.iloc[:, [1]], dc1.iloc[:, [0]])
    plt.plot(tsth20.iloc[:, [1]], tsth20.iloc[:, [0]])
    plt.plot(reesei.iloc[:, [1]], reesei.iloc[:, [0]])
    plt.plot(harzianum.iloc[:, [1]], harzianum.iloc[:, [0]])
    plt.plot(virens.iloc[:, [1]], virens.iloc[:, [0]])
    plt.axvline(x=0.28, color='black', linestyle='--')
    plt.legend(['DC1', 'Tsth20', 'T. reesei', 'T. harzianum', 'T. virens'])
    plt.xlabel('Percent GC')
    plt.ylabel('Frequency')
    plt.savefig(args.out, format='pdf')
    plt.show()
    
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script to plot GC content from windows identified by isochore.')
    #parser.add_argument('-g', type=str, dest='gc', help='path to file containing raw GC values from isochore')
    parser.add_argument('-o', type=str, dest='out', help='path to desired output file')
    args = parser.parse_args()
    main(args)
