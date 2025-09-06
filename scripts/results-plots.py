#!/usr/env python
import matplotlib.pyplot as plt
import numpy as np

print("Hello world!")

def blast_counts():
    # Data: [Fungus, Reference, Braker2, GeneMark, RefSeq]
    data = [
        ["DC1", "T. atroviride", 5902, 4679, None],
        ["DC1", "F. graminarium", 4955, 4114, None],
        ["DC1", "S. cerevisiae", 2105, 1850, None],
        ["Tsth20", "T. atroviride", 6065, 4626, None],
        ["Tsth20", "F. graminarium", 5365, 4191, None],
        ["Tsth20", "S. cerevisiae", 2211, 1869, None],
        ["T. reesei", "T. atroviride", 5072, 5174, 4989],
        ["T. reesei", "F. graminarium", 4577, 4685, 4529],
        ["T. reesei", "S. cerevisiae", 2055, 2114, 2022],
        ["T. harzianum", "T. atroviride", 6363, 4611, 6835],
        ["T. harzianum", "F. graminarium", 5659, 4198, 5982],
        ["T. harzianum", "S. cerevisiae", 2424, 1963, 2560],
        ["T. virens", "T. atroviride", 6256, 4437, 6415],
        ["T. virens", "F. graminarium", 5568, 4075, 5664],
        ["T. virens", "S. cerevisiae", 2318, 1861, 2352]
    ]

    refs_list = ["T. atroviride", "F. graminarium", "S. cerevisiae"]
    width = 0.25

    for ref in refs_list:
        ref_data = [row for row in data if row[1] == ref]
        fungi = [row[0] for row in ref_data]
        braker2 = [row[2] for row in ref_data]
        genemark = [row[3] for row in ref_data]
        refseq = [row[4] if row[4] is not None else 0 for row in ref_data]
        x = np.arange(len(fungi))

        fig, ax = plt.subplots(figsize=(8, 6))
        rects1 = ax.bar(x - width, braker2, width, label='Braker2')
        rects2 = ax.bar(x, genemark, width, label='GeneMark')
        rects3 = ax.bar(x + width, refseq, width, label='RefSeq')

        ax.set_xticks(x)
        ax.set_xticklabels(fungi, rotation=45, ha='right', fontsize=14)
        ax.set_title(f'BLAST hits to {ref}', fontsize=16)
        ax.bar_label(rects1, padding=3)
        ax.bar_label(rects2, padding=3)
        ax.bar_label(rects3, padding=3)
        ax.set_ylabel('Counts (#)', fontsize=16)
        ax.legend(fontsize=14)

        # Set y-limits to allow space for the legend
        ymax = max(max(braker2), max(genemark), max(refseq))
        ax.set_ylim(0, ymax * 1.15)

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.show()

def basic_counts():
    width = 0.25
    fungiNames = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T.virens"]
    toolNames = ["Braker2", "GeneMark", "RefSeq"]
    x = np.arange(len(fungiNames))
    multiplier = 0

    #geneData = {
    #    'Braker2': (8546, 8784, 9659, 8314, 7801),
    #    'GeneMark': (11353, 12362, 9196, 12164, 11866),
    #    'RefSeq': (0, 0, 9109, 14269, 12405)
    #}

    cdsData = {
        'Braker2': (8637, 8858, 10175, 8385, 7863),
        'GeneMark': (11353, 12362, 9196, 12164, 11866),
        'RefSeq': (0, 0, 9118, 14090, 12406)
    }

    fig, ax = plt.subplots(layout='constrained')
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
    plt.rc('font', **font)

    for attribute, measurement in cdsData.items():
        print(attribute, measurement)
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3)
        multiplier += 1
    
    ax.set_ylabel('Counts (#)')
    ax.set_title('Coding sequences predicted by different gene finding tools')
    ax.set_xticks(x + width, fungiNames)
    ax.legend(loc='upper left', ncols=3)
    ax.set_ylim(0, 15500)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    ax.set_title('Number of coding sequences predicted by different gene finding tools in Trichoderma assemblies', fontsize=16)
    ax.set_ylabel('Counts (#)', fontsize=16)
    ax.legend(loc='upper left', ncols=3, fontsize=16)

    plt.show()

    #dc1Genes = [8546, 11353, 0] 
    #dc1CDS = [8637, 11353, 0]
    #tsth20Genes = [8784, 12362, 0]
    #tsth20CDS = [8858, 12362, 0]
    #treeseiGenes = [9659, 9196, 9109]
    #treeseiCDS = [10175, 9196, 9118]
    #tharizianumGenes = [8314, 12164, 14269]
    #tharizianumCDS = [8385, 12164, 14090]
    #tvirensGenes = [7801, 11866, 12405]
    #tvirensCDS = [7863, 11866, 12406]
    #'DC1': (546, 11353, 0),
    #'Tsth20': (8784, 12362, 0),
    #'T. reesei': (9659, 9196, 9109),
    #'T. harzianum': (8314, 12164, 14269),
    #'T. virens': (7801, 11866, 12405)
    #'DC1': (8637, 11353, 0),
    #'Tsth20': (8858, 12362, 0),
    #'T. reesei': (10175, 9196, 9118),
    #'T. harzianum': (8385, 12164, 14090),
    #'T. virens': (7863, 11866, 12406)

def main():
    ### gene/CDS counts
    #basic_counts()
    ### BLAST counts
    blast_counts()

if __name__ == "__main__":
    print("This is a script for generating result plots.")
    main()



