#!/usr/env python
import matplotlib.pyplot as plt
import numpy as np

print("Hello world!")

def at_counts():
    # Data for AT-rich regions
    assemblies = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T. virens"]
    full_support = [11, 2, 25, 26, 8]
    partial_support = [0, 0, 18, 43, 11]
    singletons = [20, 9, 54, 68, 0]
    # no_genes = [42, 13, 194, 265, 49]  # Ignored in the first plot

    # Create a bar plot (excluding "No Genes")
    x = np.arange(len(assemblies))
    width = 0.18
    gap = 0.08

    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - width - gap, full_support, width, label='Full Support')
    rects2 = ax.bar(x, partial_support, width, label='Partial Support')
    rects3 = ax.bar(x + width + gap, singletons, width, label='Singletons')

    ax.set_xticks(x)
    ax.set_xticklabels(assemblies, fontsize=14)
    ax.set_ylabel('Number of Genes', fontsize=16)
    ax.set_title('Regions of Gene Predictions in AT-rich Sequence', fontsize=16)
    ax.legend(loc='upper right', fontsize=12)
    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.show()

    # Second bar plot: Assembly counts for Braker2, GeneMark, RefSeq
    assemblies = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T. virens"]
    braker2 = [31, 11, 39, 81, 21]
    genemark = [11, 2, 48, 30, 8]
    refseq = [0, 0, 107, 154, 20]  # N/A replaced with 0 for plotting

    width = 0.22
    gap = 0.04
    x = np.arange(len(assemblies))

    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - width - gap, braker2, width, label='Braker2')
    rects2 = ax.bar(x, genemark, width, label='GeneMark')
    rects3 = ax.bar(x + width + gap, refseq, width, label='RefSeq')

    ax.set_xticks(x)
    ax.set_xticklabels(assemblies, fontsize=14)
    ax.set_ylabel('Counts (#)', fontsize=16)
    ax.set_title('Counts of genes in AT-rich sequence by Assembly and Tool', fontsize=16)
    ax.legend(fontsize=14)
    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.show()


def busco_counts():
    import matplotlib.pyplot as plt

    # Data for each tool
    strains_braker2 = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T. virens"]
    complete_braker2 = [4319, 4309, 4269, 4309, 4309]
    fragmented_braker2 = [2, 6, 27, 6, 7]
    missing_braker2 = [2, 8, 27, 8, 7]

    strains_genemark = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T. virens"]
    complete_genemark = [4313, 4304, 4297, 4301, 4301]
    fragmented_genemark = [7, 10, 11, 11, 11]
    missing_genemark = [3, 9, 15, 11, 11]

    # RefSeq only for T. reesei, T. harzianum, T. virens
    strains_refseq = ["T. reesei", "T. harzianum", "T. virens"]
    complete_refseq = [4102, 4244, 4174]
    fragmented_refseq = [159, 52, 116]
    missing_refseq = [62, 27, 33]

    width = 0.22
    gap = 0.08

    # Calculate positions with gaps
    x = np.arange(len(strains_braker2))
    x_braker2 = x - width - gap
    x_genemark = x
    x_refseq = x[2:] + width + gap

    # Plot 1: Complete BUSCOs
    fig, ax = plt.subplots(figsize=(8, 6))
    rects1 = ax.bar(x_braker2, complete_braker2, width, label='Braker2')
    rects2 = ax.bar(x_genemark, complete_genemark, width, label='GeneMark')
    rects3 = ax.bar(x_refseq, complete_refseq, width, label='RefSeq')

    ax.set_xticks(x)
    ax.set_xticklabels(strains_braker2, fontsize=14)
    ax.set_ylabel('Complete BUSCOs', fontsize=16)
    ax.set_title('Complete BUSCO counts by tool and strain', fontsize=16)
    ax.legend(loc='upper right', ncol=3, fontsize=14)
    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    ax.set_ylim(0, 5000)
    # Add transparent horizontal line at 4323
    ax.axhline(4323, color='black', linestyle='--', linewidth=2, alpha=0.3)
    plt.tight_layout()
    plt.show()

    # Plot 2: Fragmented and Missing BUSCOs
    fig, ax = plt.subplots(figsize=(8, 6))
    # Fragmented
    rects1 = ax.bar(x_braker2, fragmented_braker2, width, label='Braker2 - Fragmented')
    rects2 = ax.bar(x_genemark, fragmented_genemark, width, label='GeneMark - Fragmented')
    # Use correct order for RefSeq fragmented values
    rects3 = ax.bar(x_refseq, fragmented_refseq, width, label='RefSeq - Fragmented')
    # Missing
    rects4 = ax.bar(x_braker2, missing_braker2, width, bottom=fragmented_braker2, label='Braker2 - Missing')
    rects5 = ax.bar(x_genemark, missing_genemark, width, bottom=fragmented_genemark, label='GeneMark - Missing')
    # Use correct order for RefSeq missing values
    rects6 = ax.bar(x_refseq, missing_refseq, width, bottom=fragmented_refseq, label='RefSeq - Missing')

    ax.set_xticks(x)
    ax.set_xticklabels(strains_braker2, fontsize=14)
    ax.set_ylabel('Fragmented & Missing BUSCOs', fontsize=16)
    ax.set_title('Fragmented and Missing BUSCO counts by tool and strain', fontsize=16)
    ax.legend(loc='upper left', fontsize=12, ncol=1)
    plt.tight_layout()
    plt.show()

def interproscan_counts():
    # Data: [Assembly, Braker2 (%), GeneMark (%), RefSeq (%)]
    assemblies = ["DC1", "Tsth20", "T. reesei", "T. harzianum", "T. virens"]
    braker2 = [73.73, 73.26, 72.38, 73.79, 74.68]
    genemark = [74.12, 74.10, 76.01, 74.49, 74.76]
    refseq = [0, 0, 76.44, 66.07, 73.18]  # N/A replaced with 0 for plotting

    width = 0.22
    gap = 0.04
    x = np.arange(len(assemblies))

    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - width - gap, braker2, width, label='Braker2')
    rects2 = ax.bar(x, genemark, width, label='GeneMark')
    rects3 = ax.bar(x + width + gap, refseq, width, label='RefSeq')

    ax.set_xticks(x)
    ax.set_xticklabels(assemblies, fontsize=14)
    ax.set_ylabel('Pfam matches (%)', fontsize=16)
    ax.set_title('InterProScan annotation percentages by assembly and tool', fontsize=16)
    ax.legend(fontsize=14)
    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    ax.bar_label(rects3, padding=3)
    plt.ylim(0, 100)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    plt.show()


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
    #basic_counts()
    ### BLAST counts
    #blast_region_counts()
    #blast_total_counts()
    #interproscan_counts()
    #busco_counts()
    at_counts()

if __name__ == "__main__":
    print("This is a script for generating result plots.")
    main()



