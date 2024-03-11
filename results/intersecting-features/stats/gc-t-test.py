#!/usr/bin/env python

import argparse
from BCBio import GFF
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def Ttest(data, toolCount):
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
        
                    
    print('Normal genes: ' + str(normalGC))
    print('Anomalous genes: ' + str(lowGC))
    return(normalGC, lowGC)

# Code taken from StackExchange to fit distribution of gene counts in regions
# to known distributions. Results are kind of funky. 
def distTest(data):
    dist_names = ["norm", "exponweib", "weibull_max", "weibull_min", "pareto", "genextreme"]
    dist_results = []
    params = {}
    for dist_name in dist_names:
        dist = getattr(stats, dist_name)
        param = dist.fit(data)

        params[dist_name] = param
        # Applying the Kolmogorov-Smirnov test
        D, p = stats.kstest(data, dist_name, args=param)
        print("p value for "+dist_name+" = "+str(p))
        dist_results.append((dist_name, p))

    # select the best fitted distribution
    best_dist, best_p = (max(dist_results, key=lambda item: item[1]))
    # store the name of the best fit and its p value

    print("Best fitting distribution: "+str(best_dist))
    print("Best p value: "+ str(best_p))
    print("Parameters for the best fit: "+ str(params[best_dist]))

    return best_dist, best_p, params[best_dist]

def main(args):
    for gff in args.gffFiles:
        normal, low = Ttest(gff, args.toolCount)
        df = pd.DataFrame(normal)
        print(df.describe())
        print('Variance normal: ' + str(np.var(normal)))
        print('Mean normal: ' + str(sum(normal) / len(normal)))
        print('Variance low: ' + str(np.var(low)))
        print('Mean low: ' + str(sum(low) / len(low)))
        print(stats.ttest_ind(normal, low, alternative='two-sided'))
        fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.hist(normal)
        ax1.set_title('# of genes in normal GC content regions')
        ax2.hist(low)
        ax2.set_title('# of genes in low GC content regions')
        fig.tight_layout()
        plt.show()
        #distTest(normal)
        #distTest(low)
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run T-tests on regions \
    identified by find-genomic-regions.py.')
    parser.add_argument('-f', type=str, dest='gffFiles', nargs='+', \
                        help='List of spaced GFF files for analysis', required=True)
    parser.add_argument('-n', type=int, dest='toolCount', help='Number of tools considered in input to region identification. (Gene finders, gc content, repeat regions, etc.)', required=True)
    args = parser.parse_args()
    main(args)
