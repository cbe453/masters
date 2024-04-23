#!/usr/bin/env python

import argparse
from Bio import SeqIO

def main(args):
    
    fileTypes = ['perfect', 'perfect.spacer']
    
    
    gffEntries = []
    for repType in ('./tir/', './tdr/'):
        for fileType in fileTypes:
            count = 0
            fasta = repType + fileType + '.fasta'
            for record in SeqIO.parse(fasta, 'fasta'):
                splitStr = record.id.split(':')
                count += 1
                lineID = (fileType + '.' + str(count))
                outStr = str("%s\tGRF\tTIR\t%s\t%s\t.\t.\t.\t%s;type=%s:match=%s\n" % (splitStr[0], splitStr[1], splitStr[2], lineID, fileType, splitStr[3]))
                gffEntries.append(outStr)

        print(str(len(gffEntries)))

    count = 0
    for record in SeqIO.parse('./mite/candidate.fasta', 'fasta'):
        splitStr = record.id.split(':')
        count += 1
        lineID = (fileType + '.' + str(count))
        outStr = str("%s\tGRF\tTIR\t%s\t%s\t\t.\t.\t%s;type=mite:match=%s\n" % (splitStr[0], splitStr[1], splitStr[2], lineID, splitStr[3]))
        gffEntries.append(outStr)

    print(str(len(gffEntries)))

    with open(args.out, 'w') as handle:
        for entry in gffEntries:
            handle.write(entry)
        
    return(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to convert fasta output from Generic Repeat Finder to GFF format for further analysis.')
    parser.add_argument('-d', dest='basedir', help='Directory containing the MITE, TIR and TDR output directories from GRF.')
    parser.add_argument('-o', dest='out', help='Name of resulting GFF output.')
    args = parser.parse_args()
    main(args)
