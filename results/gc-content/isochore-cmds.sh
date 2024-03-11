# extract contig IDs from assembly
grep '^>' contigs.fa > contig-ids.txt

# example isochore command
for contig in `cat dc1-contig-ids.txt`; do echo $contig; outfile=`echo $contig | sed 's/$/.out/'`; isochore -sequence ../../IGV/dc1-igv/dc1-genome.nextpolish.fmt.fasta:$contig -outfile $outfile -graph png -goutfile $contig -window 250 -shift 50; done

# parse output from isochore into GFF format
for file in `cat dc1-contig-ids.txt`; do awk -v infile=$file -v id=1 '{if ($2 < 0.28) {print infile, "isochore\tAT-rich", ($1 - 124), $1, ".\t-\t.", "ID=" id++ ";percent=" $2}}' OFS='\t' $file.out ; done | sort -k 1 > .gff
