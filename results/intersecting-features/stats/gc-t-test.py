#!/usr/bin/env python

import argparse
from BCBio import GFF
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def countGenes(data, toolCount):
    gffHandle = open(data)
    normalGC = []
    lowGC = []
    
    # Iterate over contigs from supplied GFF file
    for contig in GFF.parse(gffHandle):
        geneCount = 0
        genes = 0
        previousFlag = None
        feature = None
        currentFlag = None
        skip = None

        # Iterate over features (regions) from each contig
        for feature in contig.features:
            genes = int(feature.qualifiers['geneCount'][0])
            rsmiList = ["'AUGUSTUS'", 'GeneMark']
            #refseqList = ['AUGUSTUS', 'GeneMark', 'RefSeq']
            currentFlag = feature.qualifiers['lowGC'][0]
            skip = 'False'

            # Check if all tools in region agree.
            if int(feature.qualifiers['toolCount'][0]) == toolCount:
                skip = 'True'
                continue
            # Check if all gene finers in region agree.
            elif int(feature.qualifiers['toolCount'][0]) == (toolCount -1 ) and ('isochore' not in str(feature.qualifiers['tools'])):
                skip = 'True'
                continue
            # Check GC content flags from GFF file in next chunks.
            elif (currentFlag == 'False'):
                if previousFlag == 'True':
                    lowGC.append(geneCount)
                    previousFlag = currentFlag
                    geneCount = genes
                else:
                    geneCount += genes
                    previousFlag = currentFlag
            elif (currentFlag == 'True'):
                if previousFlag == 'False':
                    normalGC.append(geneCount)
                    previousFlag = currentFlag
                    geneCount = genes
                else:
                    geneCount += genes
                    previousFlag = currentFlag

        # End condition where final region is a complete agreement case from
        # from above.
        if (skip == 'True'):
            if (previousFlag == 'True'):
                lowGC.append(geneCount)
            else:
                normalGC.append(geneCount)
        # If final region has disagreement, append gene counts appropriately.
        else:
            if (currentFlag == previousFlag):
                if (currentFlag == 'True'):
                    lowGC.append(geneCount)
                else:
                    normalGC.append(geneCount)
            else:
                if (currentFlag == 'True'):
                    normalGC.append(geneCount)
                    lowGC.append(genes)
                else:
                    lowGC.append(geneCount)
                    normalGC.append(genes)

            #print(contig)
            #print(feature)
            #print(str(geneCount))
        
                    
    #print('Normal genes: ' + str(normalGC))
    #print('Anomalous genes: ' + str(lowGC))
    return(normalGC, lowGC)

def dc1Chi2():
    normal, low = (countGenes('/Users/cbe453/Desktop/masters/masters/masters/results/intersecting-features/dc1/genomic-regions/gc-stats.gff', 3))
    totalLength = 38616239
    lowGcLength = 2064693
    lowGcFraction = ((lowGcLength / totalLength))
    totalGenes = (sum(normal) + sum(low))
    expLowGcGenes = int(totalGenes * lowGcFraction)
    expNormGcGenes = totalGenes - expLowGcGenes

    print('Normal Genes: ' + str(sum(normal)))
    print('% assembly low GC: ' + str(lowGcFraction))
    print('Total Gene: ' + str(totalGenes))
    print('Expected Low GC genes: ' + str(expLowGcGenes))
    print('Expected Normal GC genes: ' + str(expNormGcGenes))

    #chisquare
    
def main(args):
    
    dc1Chi2()
    #for gff in args.gffFiles:
        #normal, low = Ttest(gff, args.toolCount)
        


        #df = pd.DataFrame(normal)
        #print(df.describe())
        #print('Variance normal: ' + str(np.var(normal)))
        #print('Mean normal: ' + str(sum(normal) / len(normal)))
        #print('Variance low: ' + str(np.var(low)))
        #print('Mean low: ' + str(sum(low) / len(low)))
        #print(stats.ttest_ind(normal, low, alternative='two-sided'))
        #fig, (ax1, ax2) = plt.subplots(2, 1)
        #ax1.hist(normal)
        #ax1.set_title('# of genes in normal GC content regions')
        #ax2.hist(low)
        #ax2.set_title('# of genes in low GC content regions')
        #fig.tight_layout()
        #plt.show()
        #distTest(normal)
        #distTest(low)
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run stats on regions \
    identified by find-genomic-regions.py.')
    #parser.add_argument('-f', type=str, dest='gffFiles', nargs='+', \
    #                    help='List of spaced GFF files for analysis', required=True)
    parser.add_argument('-n', type=int, dest='toolCount', help='Number of tools considered in input to region identification. (Gene finders, gc content, repeat regions, etc.)', required=True)
    args = parser.parse_args()
    main(args)
