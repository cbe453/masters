#!/usr/bin/env python
import argparse
import numpy as np

def ctDC1():
    totalLength = 38616239
    lowGcLength = 2064693
    lowGcFraction = ((lowGcLength / totalLength) * 100)
    print(str(lowGcFraction))
    totalGenes = 19858 + 42
    
def main(args):
    print(args)
    output = ''
    ctDC1()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Chi squared test on regions \
    gene predictions in low GC content regions..')
    parser.add_argument('-o', type=str, dest='output', help='List of spaced GFF files for analysis', required=True))
    args = parser.parse_args()
    main(args)
