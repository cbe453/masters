#!/usr/env python
import matplotlib.pyplot as plt
import numpy as np

print("Hello world!")

def basic_counts():
    barWidth = 0.25
    fungiNames = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T.virens"]
    toolNames = ["Braker2", "GeneMark", "RefSeq"]
    
    geneData = [
        [8546, 11353, 0],
        [8784, 12362, 0],
        [9659, 9196, 9109],
        [8314, 12164, 14269],
        [7801, 11866, 12405]
    ]
    cdsData = [[8637, 11353, 0],
               [8858, 12362, 0],
               [10175, 9196, 9118],
               [8385, 12164, 14090],
               [7863, 11866, 12406]]

    dc1Genes = [8546, 11353, 0]
    dc1CDS = [8637, 11353, 0]
    tsth20Genes = [8784, 12362, 0]
    tsth20CDS = [8858, 12362, 0]
    treeseiGenes = [9659, 9196, 9109]
    treeseiCDS = [10175, 9196, 9118]
    tharizianumGenes = [8314, 12164, 14269]
    tharizianumCDS = [8385, 12164, 14090]
    tvirensGenes = [7801, 11866, 12405]
    tvirensCDS = [7863, 11866, 12406]

    pos1 = np.arange(len(dc1Genes))
    pos2 = pos1 + barWidth
    pos3 = pos2 + barWidth

    fig, ax = plt.subplots(3,2)
    fig.delaxes(ax[2][1])
    ax[0,0].bar(pos1, dc1Genes, width=barWidth, edgecolor='white', label='Genes')
    ax[0,0].bar(pos2, dc1CDS, width=barWidth, edgecolor='white', label='CDS')
    #ax[0,0].set_xlabel('Tools', fontweight='bold')
    ax[0,0].set_xticks(pos1 + barWidth)
    ax[0,0].set_xticklabels(['Braker2', 'GeneMark', 'RefSeq'])
    ax[0,0].set_ylim(0,15000)

    ax[0,1].bar(pos1, tsth20Genes, width=barWidth, edgecolor='white', label='Genes')
    ax[0,1].bar(pos2, tsth20CDS, width=barWidth, edgecolor='white', label='CDS')
    #ax[0,0].set_xlabel('Tools', fontweight='bold')
    ax[0,1].set_xticks(pos1 + barWidth)
    ax[0,1].set_xticklabels(['Braker2', 'GeneMark', 'RefSeq'])
    ax[0,1].set_ylim(0,15000)

    ax[1,0].bar(pos1, treeseiGenes, width=barWidth, edgecolor='white', label='Genes')
    ax[1,0].bar(pos2, treeseiCDS, width=barWidth, edgecolor='white', label='CDS')
    #ax[0,0].set_xlabel('Tools', fontweight='bold')
    ax[1,0].set_xticks(pos1 + barWidth)
    ax[1,0].set_xticklabels(['Braker2', 'GeneMark', 'RefSeq'])
    ax[1,0].set_ylim(0,15000)

    ax[1,1].bar(pos1, tharizianumGenes, width=barWidth, edgecolor='white', label='Genes')
    ax[1,1].bar(pos2, tharizianumCDS, width=barWidth, edgecolor='white', label='CDS')
    #ax[0,0].set_xlabel('Tools', fontweight='bold')
    ax[1,1].set_xticks(pos1 + barWidth)
    ax[1,1].set_xticklabels(['Braker2', 'GeneMark', 'RefSeq'])
    ax[1,1].set_ylim(0,15000)

    ax[2,0].bar(pos1, tvirensGenes, width=barWidth, edgecolor='white', label='Genes')
    ax[2,0].bar(pos2, tvirensCDS, width=barWidth, edgecolor='white', label='CDS')
    #ax[0,0].set_xlabel('Tools', fontweight='bold')
    ax[2,0].set_xticks(pos1 + barWidth)
    ax[2,0].set_xticklabels(['Braker2', 'GeneMark', 'RefSeq'])
    ax[2,0].set_ylim(0,15000)

    ax[0,0].legend()
    plt.show()

def main():
    ### gene/CDS counts
    basic_counts()

if __name__ == "__main__":
    print("This is a script for generating result plots.")
    main()



