import re

genomes = {
    "dc1": "C:\\Users\\cbe45\\OneDrive\\Desktop\\Masters\\masters\\results\\IGV\\dc1-igv\\dc1-genome.nextpolish.fmt.fasta",
    "tsth20": "C:\\Users\\cbe45\\OneDrive\\Desktop\\Masters\\masters\\results\\IGV\\tsth20-igv\\tsth20-genome.nextpolish.fmt.fasta",
    "tReesei": "C:\\Users\\cbe45\\OneDrive\\Desktop\\Masters\\masters\\results\\IGV\\t-reesei-igv\\GCF_000167675.1_v2.0_genomic.fna",
    "tHarz": "C:\\Users\\cbe45\\OneDrive\\Desktop\\Masters\\masters\\results\\IGV\\t-harzianum-igv\\GCF_003025095.1_Triha_v1.0_genomic-fmt.fna",
    "tVirens": "C:\\Users\\cbe45\\OneDrive\\Desktop\\Masters\\masters\\results\\IGV\\t-virens-igv\\GCF_000170995.1_TRIVI_v2.0_genomic-fmt.fna"
}

gffFiles = {
    "dc1": "X:\\Connor\\masters\\results\\intersecting-features\\dc1\\gc-content\\dc1-intersecting-gc-genes.gff",
    "tsth20": "X:\\Connor\\masters\\results\\intersecting-features\\tsth20\\gc-content\\tsth20-intersecting-gc-genes.gff",
    "tReesei": "X:\\Connor\\masters\\results\\intersecting-features\\t-reesei\\gc-content\\t-reesei-intersecting-gc-genes.gff",
    "tHarz": "X:\\Connor\\masters\\results\\intersecting-features\\t-harzianum\\gc-content\\t-harzianum-intersecting-gc-genes.gff",
    "tVirens": "X:\\Connor\\masters\\results\\intersecting-features\\t-virens\\gc-content\\t-virens-intersecting-gc-genes.gff"
}

def calculate_genome_length(fasta_file):
    """Calculate total length of sequences in a FASTA file"""
    total_length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length

def calculate_lowgc_stats(gff_file):
    """Calculate total length and count of regions with lowGC=True in GFF file"""
    total_length = 0
    count = 0
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                attributes = fields[8]
                if 'lowGC=True' in attributes:
                    start = int(fields[3])
                    end = int(fields[4])
                    total_length += (end - start + 1)
                    count += 1
    return total_length, count

# Calculate genome lengths and lowGC proportions
# Print table header
print(f"{'Genome':<15} {'Total Length (bp)':>20} {'Low GC Count':>15} {'Low GC Length (bp)':>20} {'Proportion (%)':>15}")
print("-" * 90)

# Calculate genome lengths and lowGC proportions
for genome_name in genomes:
    genome_length = calculate_genome_length(genomes[genome_name])
    lowgc_length, lowgc_count = calculate_lowgc_stats(gffFiles[genome_name])
    proportion = (lowgc_length / genome_length) * 100 if genome_length > 0 else 0
    
    print(f"{genome_name:<15} {genome_length:>20,} {lowgc_count:>15} {lowgc_length:>20,} {proportion:>14.2f}%")


