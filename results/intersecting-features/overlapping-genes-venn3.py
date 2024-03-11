#!/usr/bin/env/ python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn3
from matplotlib_venn import venn3_unweighted

def main(args):
    print("sup")
    print(args.counts)
    countTable = pd.read_csv(args.counts, sep='\t', header=None)
    venn3_unweighted(subsets = (countTable.iloc[:, 1]), set_labels = ('RefSeq', 'Braker2', 'GeneMark'))
    plt.title('T. virens')
    plt.show()
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analysis of overlaping genes based on counts from bedtools intersect tool')
    parser.add_argument('-d', type=str, dest='counts', help='Path to the file containing counts for genes for each caller along with counts for the bedtools intersect command.')
    args = parser.parse_args()
    exit(main(args))
