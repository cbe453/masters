#!/usr/env python

import argparse
from BCBio import GFF 

def parseIPS(ips):
    ipsDict = {}
    ipsDict['declare'] = []
    ipsHandle = open(args.ipsGFF)
    for record in GFF.parse(ipsHandle):
        for item in record.features:
            if ('source' in item.qualifiers):
                if str(item.qualifiers['source'][0]) == 'Pfam':
                    geneSource = item.qualifiers['Target'][0]
                    geneNum = str(geneSource.split(' ')[0].split('_')[1])
                    match = item.qualifiers['ID'][0]
                    if (geneNum in ipsDict):
                        ipsDict[geneNum].append(item)
                    else:
                        ipsDict[geneNum] = [item]
    ipsHandle.close()
    return(ipsDict)
    
def parseGeneMark(args, ipsDict):
    gmHandle = open(args.geneGFF)
    outHandle = open(args.outGFF, 'w')
    for record in GFF.parse(gmHandle):
        for item in record.features:
            #print(item)
            ipsID = str(item.qualifiers['ID'][0].split('_')[1])
            if ((ipsID in ipsDict) and(item.type == 'gene')):
                ipsLine = ipsDict[ipsID]
                for ipsLine in ipsDict[ipsID]:
                    if (item.strand == 1):
                        newStart = item.location.start + ipsLine.location.start
                        newEnd = item.location.start + ipsLine.location.end
                        strand = '+'
                    elif (item.strand == -1):
                        newEnd = item.location.end - ipsLine.location.start
                        newStart = item.location.end - ipsLine.location.end
                        strand = '-'
                    
                        contig = record.id.split(' ')[0].split('_')[0]
                        source = 'InterProScan'
                        ipsType = ipsLine.qualifiers['source'][0]
                        ogIpsID = ipsLine.qualifiers['ID'][0]
                        match = ipsLine.qualifiers['Name'][0]
                        ogGene = item.qualifiers['ID'][0]
                        mRNA = ipsLine.qualifiers['Target'][0]
                        score = ipsLine.qualifiers['score'][0]
                        desc = ipsLine.qualifiers['signature_desc'][0]

                        outStr = (contig + '\t' + source + '\t' + ipsType + '\t' \
                                  + str(newStart) + '\t' + str(newEnd) + '\t' \
                                  + '.' + '\t' + strand + '\t' + '.\tID=' + ogIpsID + ';score=' \
                                  + str(score) + ';match=' + match + ';desc=\'' + desc \
                                  + '\';mRNA=' + mRNA + ';ogGene=' + ogGene + ';ipsID=' + ipsID + '\n')
                        #print(outStr)
                        outHandle.write(outStr)
    gmHandle.close()
    outHandle.close()
        
def main(args):
    print(args)
    ipsDict = parseIPS(args.ipsGFF)
    if (args.finder == 'GeneMark'):
        parseGeneMark(args, ipsDict)
    else:
        print('Argument to -f unknown! Please supply one of Braker2, GeneMark or RefSeq.')
    return(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge interproscan Pfam hits with genes in predicted genes')
    parser.add_argument('-f', type=str, dest='finder', help='Name of the gene finder used for supplied predictions. Can be one of [Braker2, GeneMark, RefSeq].', required=True)
    parser.add_argument('-i', type=str, dest='ipsGFF', help='GFF file output by InterProScan', required=True)
    parser.add_argument('-g', type=str, dest='geneGFF', help='GFF with predicted genes.', required=True)
    parser.add_argument('-o', type=str, dest='outGFF', help='Path and anme of the newly mappend IPS GFF file', required=True)
    args = parser.parse_args()
    main(args)
