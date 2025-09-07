#!/usr/env python
import matplotlib.pyplot as plt
import numpy as np

print("Hello world!")

def blast_total_counts():
    # Data: [Fungus, DC1, Tsth20, T. reesei, T. harzianum, T. virens, Average]
    data = [
        ["T. atroviride", 11552, 11080, 10601, 11081, 11078],
        ["F. graminarium", 10327, 10429, 10064, 10434, 10490],
        ["S. cerevisiae", 3537, 3517, 3445, 3509, 3500]
    ]

    fungi = [row[0] for row in data]
    braker2 = [row[1] for row in data]
    genemark = [row[2] for row in data]
    treesei = [row[3] for row in data]
    tharzianum = [row[4] for row in data]
    tvirens = [row[5] for row in data]

    width = 0.13  # slightly narrower bars
    gap = 0.03    # space between bars
    x = np.arange(len(fungi))

    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - 2*width - 2*gap, braker2, width, label='DC1')
    rects2 = ax.bar(x - width - gap, genemark, width, label='Tsth20')
    rects3 = ax.bar(x, treesei, width, label='T. reesei')
    rects4 = ax.bar(x + width + gap, tharzianum, width, label='T. harzianum')
    rects5 = ax.bar(x + 2*width + 2*gap, tvirens, width, label='T. virens')

    ax.set_xticks(x)
    ax.set_xticklabels(fungi, fontsize=14)
    ax.set_ylabel('Counts (#)', fontsize=16)
    ax.set_title('Total tblastn hits for selected Trichoderma assemblies ', fontsize=16)
    ax.legend(fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    ax.bar_label(rects4, padding=3)
    ax.bar_label(rects5, padding=3)
    plt.tight_layout()
    plt.show() 

def blast_region_counts():
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
    width = 0.22
    gap = 0.04

    for ref in refs_list:
        ref_data = [row for row in data if row[1] == ref]
        fungi = [row[0] for row in ref_data]
        braker2 = [row[2] for row in ref_data]
        genemark = [row[3] for row in ref_data]
        refseq = [row[4] if row[4] is not None else 0 for row in ref_data]
        x = np.arange(len(fungi))

        fig, ax = plt.subplots(figsize=(8, 6))
        rects1 = ax.bar(x - width - gap, braker2, width, label='Braker2')
        rects2 = ax.bar(x, genemark, width, label='GeneMark')
        rects3 = ax.bar(x + width + gap, refseq, width, label='RefSeq')

        ax.set_xticks(x)
        ax.set_xticklabels(fungi, rotation=45, ha='right', fontsize=14)
        ax.set_title(f'BLAST hits to {ref}', fontsize=16)
        ax.bar_label(rects1, padding=3)
        ax.bar_label(rects2, padding=3)
        ax.bar_label(rects3, padding=3)
        ax.set_ylabel('Counts (#)', fontsize=16)
        ax.legend(fontsize=14)

        ymax = max(max(braker2), max(genemark), max(refseq))
        ax.set_ylim(0, ymax * 1.15)

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.tight_layout()
        plt.show()

def basic_counts():
    width = 0.22
    gap = 0.04
    fungiNames = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T.virens"]
    toolNames = ["Braker2", "GeneMark", "RefSeq"]
    x = np.arange(len(fungiNames))

    geneData = {
        'Braker2': (8546, 8784, 9659, 8314, 7801),
        'GeneMark': (11353, 12362, 9196, 12164, 11866),
        'RefSeq': (0, 0, 9109, 14269, 12405)
    }

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

    offsets = [-width-gap, 0, width+gap]
    for i, (attribute, measurement) in enumerate(geneData.items()):
        print(attribute, measurement)
        rects = ax.bar(x + offsets[i], measurement, width, label=attribute)
        ax.bar_label(rects, padding=3)
    
    ax.set_ylabel('Counts (#)')
    ax.set_title('Coding sequences predicted by different gene finding tools')
    ax.set_xticks(x)
    ax.set_xticklabels(fungiNames)
    ax.legend(loc='upper left', ncols=3)
    ax.set_ylim(0, 15500)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    ax.set_title('Number of genes predicted', fontsize=16)
    ax.set_ylabel('Counts (#)', fontsize=16)
    ax.legend(loc='upper left', ncols=3, fontsize=16)

    plt.show()

def main():
    ### gene/CDS counts
    basic_counts()
    ### BLAST counts
    #blast_region_counts()
    #blast_total_counts()

if __name__ == "__main__":
    print("This is a script for generating result plots.")
    main()



