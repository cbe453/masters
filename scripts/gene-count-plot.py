#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def main():

    DC1 = [8546, 8637, 11353, 11353]
    Tsth20 = [8784, 8858, 12362, 12362]
    Tree = [9659, 10175, 9203, 9203]
    Tharz = [8314, 8385, 12164, 12164]
    Tvir = [7801, 7863, 11877, 11877]

    brakerGene = [8546, 8784, 9659, 8314, 7801]
    brakerTran = [8637, 8858, 10175, 8385, 7863]
    genemarkGene = [11353, 12362, 9203, 12164, 11877]
    genemarkTran = [11353, 12362, 9203, 12164, 11877]
    labels = ['DC1', 'Tsth20', 'T. reesei', 'T. harzianum', 'T. virens']
    
    index = np.arange(5)
    barWidth = 0.2

    plt.bar(index - 0.3, brakerGene, barWidth, alpha=0.7, label='Braker Genes')
    plt.bar(index - 0.1, brakerTran, barWidth, alpha=0.7, label='Braker Transcripts')
    plt.bar(index + 0.1, genemarkGene, barWidth, alpha=0.7, label='GeneMark Genes')
    plt.bar(index + 0.3, genemarkTran, barWidth, alpha=0.7, label='GeneMark Transcripts')
    plt.ylabel('Gene/Transcript Counts')
    plt.xlabel('Assembly')
    plt.xticks(index, labels)
    plt.legend(loc='upper center', ncol=2)
    plt.ylim(0, 15000)
    plt.tight_layout()
    #plt.show()
    plt.savefig('/Users/cbe453/Desktop/masters/masters/working-thesis/figures/genecounts.pdf', format='pdf')
    
if __name__ == '__main__':
    main()
