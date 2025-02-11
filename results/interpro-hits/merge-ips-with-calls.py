#!/usr/env python

import argparse
import gffutils
from BCBio import GFF

def parseIPS(tool, ips):
    ipsDict = {}
    ipsDict['declare'] = []
    ipsHandle = open(args.ipsGFF)
    for record in GFF.parse(ipsHandle):
        for item in record.features:
            if ('source' in item.qualifiers):
                if str(item.qualifiers['source'][0]) == 'Pfam':
                    geneSource = item.qualifiers['Target'][0]
                    if (tool == 'GeneMark'):
                        geneNum = str(geneSource.split(' ')[0].split('_')[1])
                    elif (tool == 'Braker2'):
                        geneNum = str(geneSource.split(' ')[0].split('.')[0])
                    else:
                        geneNum = str(geneSource.split('.')[0].split('_')[1])
                        
                    match = item.qualifiers['ID'][0]
                    if (geneNum in ipsDict):
                        ipsDict[geneNum].append(item)
                    else:
                        ipsDict[geneNum] = [item]

    del ipsDict['declare']
    ipsHandle.close()
    return(ipsDict)


def parseRefSeq(args, ipsDict):
    db = gffutils.create_db(args.geneGFF, dbfn='t-reesei-refseq.db', force=True, keep_order=True,
                            merge_strategy='merge', sort_attribute_values=True)

    db = gffutils.FeatureDB('t-reesei-refseq.db', keep_order=True)
    outHandle = open(args.outGFF, 'w')
    
    for key in ipsDict.keys():
        cdsKey = 'cds-XP_' + key + '.1'
        try:
            cds = db[cdsKey]
            #print(cds.attributes['Parent'])
        except:
        
            print('No CDS found for protein XP_' + key)

        rnaKey = cds.attributes['Parent'][0]
        rna = db[rnaKey]
        #print(rna.attributes['Parent'][0])
        #print('Start: ' + str(rna.start) + '  End: ' + str(rna.end) + '   Strand: ' + str(rna.strand))
        #print(len(ipsDict[key]))

        for pfam in ipsDict[key]:
            #print(str(pfam.location.start) + '\t' + str(pfam.location.end))
            if rna.strand == '+':
                newStart = rna.start + pfam.location.start
                newEnd = rna.start + pfam.location.end
                strand = '+'
            else:
                newStart = rna.end - pfam.location.end
                newEnd = rna.end - pfam.location.start
                strand = '-'

            contig = rna.seqid
            source = 'InterProScan'
            ipsType = pfam.qualifiers['source'][0]
            ogIpsID = pfam.qualifiers['ID'][0]
            match = pfam.qualifiers['Name'][0]
            ogGene = rna.attributes['Parent'][0]
            mRNA = cdsKey
            score = pfam.qualifiers['score'][0]
            desc = pfam.qualifiers['signature_desc'][0]
            ipsID = 'XP_' + key
        
            outStr = (contig + '\t' + source + '\t' + ipsType + '\t' \
                      + str(newStart) + '\t' + str(newEnd) + '\t' \
                      + '.' + '\t' + strand + '\t' + '.\tID=' + ogIpsID + ';score=' \
                      + str(score) + ';match=' + match + ';desc=\'' + desc \
                      + '\';mRNA=' + mRNA + ';ogGene=' + ogGene + ';ipsID=' + ipsID + '\n')
            #print(outStr)
            outHandle.write(outStr)
    outHandle.close()

    
def parseGeneMark(args, ipsDict):
    tool = args.finder
    gmHandle = open(args.geneGFF)
    outHandle = open(args.outGFF, 'w')
    for record in GFF.parse(gmHandle):
        for item in record.features:
            if (tool == 'GeneMark'):
                ipsID = str(item.qualifiers['ID'][0].split('_')[1])
            elif (tool == 'Braker2'):
                ipsID = str(item.qualifiers['ID'][0].split('.')[0])
            else:
                print('Something went wrong...')
    
            if ((ipsID in ipsDict) and (item.type == 'gene')):
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

                    #DC1 / Tsth20
                    #contig = record.id.split(' ')[0].split('_')[0]
                    #RefSeq assemblies
                    contig = record.id
                        
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
    #print(args)
    if args.finder not in ['Braker2', 'GeneMark', 'RefSeq']:
        print('Argument to -f unknown. Please supply one of Braker2, GeneMark or RefSeq.')
        return(1)
    elif (args.finder == 'RefSeq'):
        ipsDict = parseIPS(args.finder, args.ipsGFF)
        parseRefSeq(args, ipsDict)
    else:
        print('Creating new GFF file with InterProScan hits from ' + args.finder + ' mapped to original gene predictions.')
        ipsDict = parseIPS(args.finder, args.ipsGFF)
        print(len(ipsDict.keys()))
        parseGeneMark(args, ipsDict)

    return(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge interproscan Pfam hits with genes in predicted genes')
    parser.add_argument('-f', type=str, dest='finder', help='Name of the gene finder used for supplied predictions. Can be one of [Braker2, GeneMark, RefSeq].', required=True)
    parser.add_argument('-i', type=str, dest='ipsGFF', help='GFF file output by InterProScan', required=True)
    parser.add_argument('-g', type=str, dest='geneGFF', help='GFF with predicted genes.', required=True)
    parser.add_argument('-o', type=str, dest='outGFF', help='Path and anme of the newly mappend IPS GFF file', required=True)
    args = parser.parse_args()
    main(args)
