#!/usr/bin/env/ python

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import random
import seaborn as sns

# Nice color palette for colorblind individuals
#CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
#                  '#f781bf', '#a65628', '#984ea3',
#                  '#999999', '#e41a1c', '#dede00']

def histPlot(args):
    print(args)

    braker = parseData(args.braker)
    refseq = parseData(args.refseq)
    genemark = parseData(args.genemark)

    tools = ['Braker2', 'RefSeq', 'GeneMark']
    cds = []
    cds.append(braker)
    cds.append(refseq)
    cds.append(genemark)
    pvalues = []
    
    for i in range(0, len(cds)):
        for j in range(i+1, len(cds)):
            print(tools[i] + ' ' + tools[j])
            pvalues.append(stats.ks_2samp(cds[i], cds[j]).pvalue)

    print(pvalues)
    #fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    #axes = [ax1, ax2, ax3]

    #for i in range(0, len(cds)):
    #    axes[i].hist(cds[i], bins=1000)
    #    axes[i].set_xlim(0,5000)
    #    axes[i].set_ylim(0, 500)
    #    axes[i].set_xlabel(tools[i])

    #ax2.set_yticks([])
    size = 75
    sns.histplot(cds[0], bins=range(0, max(cds[0]) + size, size), kde=True, color='#e41a1c', label=tools[0])
    sns.histplot(cds[1], bins=range(0, max(cds[0]) + size, size), kde=True, color='#377eb8', label=tools[1])
    sns.histplot(cds[2], bins=range(0, max(cds[0]) + size, size), kde=True, color='#dede00', label=tools[2])
    #plt.ylim(0, 600)
    plt.xlim(0, 5500)
    plt.xlabel('Gene Length (bp)')
    plt.tight_layout()
    plt.legend()
    plt.savefig(args.output)
    plt.show()
    
def cumDist(args):
    tools = ['Braker2', 'RefSeq', 'GeneMark']
    for i in range(0, len(args.input)):
        data = parseData(args.input[i])
        sortedData = np.sort(data)
        finalSum = sum(data)
        curSum = 0
        i = 0
        x = []
        y = []
        for cdsLen in sortedData:
            curSum = cdsLen + curSum
            plotVal = (curSum / finalSum)
            i += 1
            y.append(plotVal)
            x.append(i)
        #sortedData = np.sort(finalData)
        print(x[-5::1])
        print(y[-5::1])
        print(curSum)
        print(finalSum)
        #plt.step(finalData, np.arange(0, len(finalData)))
        #sns.histplot(finalData, kde=True, color='#dede00', label=tools[2])
        plt.plot(y, x)
    plt.show()
        
def parseData(cdsPath):
    cdsArray = []
    for cds in SeqIO.parse(cdsPath, "fasta"):
        cdsArray.append(len(cds.seq))

    return(cdsArray)

def main(args):
    cumDist(args)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analysis of gene lengths from gene finding tools')
    parser.add_argument('-i', type=str, dest='input', nargs='+', help='List of CDS sequence files to be plotted.')
    parser.add_argument('-o', type=str, dest='output', help='Output file name')
    args = parser.parse_args()
    
    
    exit(main(args))
