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
from collections import Counter
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

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
    tools = ['Braker2', 'GeneMark', 'RefSeq']
    #tools = ['Braker2', 'GeneMark']
    fig, ax = plt.subplots()
    #axins = ax.inset_axes([0.2, 0.13, 0.51, 0.7])

    for i in range(0, len(args.input)):
        if (len(tools) > 2):
            next
        data = parseData(args.input[i])
        sortedData = sorted(np.log10(data))
        finalSum = sum(data)
        finalCount = len(data)
        x = []
        p = []
        counter = Counter(sortedData)

        for key in counter.keys():
            x.append(key)
            p.append(counter[key] / finalCount)

        y = np.cumsum(p)
        
        ax.plot(x, y, label=tools[i])
        #axins.plot(x, y, label=tools[i])
        #sns.histplot(finalData, kde=True, color='#dede00', label=tools[2])
        #plt.step(x, y, label = tools[i])

    #mark_inset(ax, axins, loc1=1, loc2=4, fc='none', ec='0.5')
    #axins.set_xlabel('')
    #axins.set_xlim(0, 2500)
    #axins.set_ylim(0,1)
    plt.ylabel('Proportion of total genes (%)')
    plt.xlabel('Log10 Gene Length')
    plt.legend(loc = 'center right')
    plt.savefig(args.output, format = 'pdf')
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
