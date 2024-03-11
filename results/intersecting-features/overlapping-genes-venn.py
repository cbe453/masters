#!/usr/bin/env/ python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2

def main(args):
    print(args.counts)
    countTable = pd.read_csv(args.counts, sep='\t', header=None)
    labels = ['Braker', 'GeneMark']

    v = venn2(subsets = (countTable.iloc[:, 1]), set_labels = labels)
    v.hide_zeroes()
    labels = v.set_labels
    print(labels[0][0])
    
    plt.title("DC1")
    #plt.show()
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analysis of overlaping genes based on counts from bedtools intersect tool')
    parser.add_argument('-d', type=str, dest='counts', help='Path to the file containing counts for genes for each caller along with counts for the bedtools intersect command.')
    args = parser.parse_args()
    exit(main(args))
