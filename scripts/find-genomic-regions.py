#!/User/bin/env python
# Author: Connor Burbridge
# Affiliation: University of Saskatchewan, Computer Science.
# Supervisors: Dave Schneider and Tony Kusalik.
# Purpose: program to identify 'regions' of overlapping features from one or multiple GFF
# files for MSc.

import argparse
from BCBio import GFF


# Class definition for a region. Will likely need to be expanded on
# for future analysis
class region:
    def __init__(self, regionID, start, end, features, programs):
        self.region = ('region' + str(regionID))
        self.startPos = start
        self.endPos = end
        self.featureList = features
        self.tools = programs

    def printRegion(self):
        itemStr = ''
        for item in self.featureList:
            #print('\t' + item.feature + '\t' + str(item.startPos) + '\t' + str(item.endPos))
            itemStr += (';' + item.feature + ':' + item.startPos)
        print(str(self.region) + '\t' + str(self.startPos) + '\t' + str(self.endPos) + itemStr)
            
            
# Class definition for a feature from a GFF file. Will likely need to
# be expanded on for future analysis
class feature:
    def __init__(self, featureID, contig, start, end, strand, program, gcFlag):
        self.feature = featureID
        self.contig = contig
        self.startPos = start
        self.endPos = end
        self.strand = strand
        self.tool = program
        self.lowGC = gcFlag
        #self.overlap = None
        
    def printFeature(self):
        print(self.feature + '\t' + str(self.startPos) + '\t' + str(self.endPos) \
              + '\t' + str(self.strand))
    

# Original testing. Kept in case something goes awry
def testFindRegion():
    feature1 = feature('feature1', 1, 10)
    feature2 = feature('feature2', 5, 15)
    feature3 = feature('feature3', 7, 20)
    feature4 = feature('feature4', 35, 45)
    feature5 = feature('feature5', 1, 45)
    feature6 = feature('feature6', 75, 90)
    feature7 = feature('feature7', 78, 87)
    feature8 = feature('feature8', 1, 45)
    feature9 = feature('feature9', 100, 125)
    featureList = [feature3, feature1, feature9, feature5, feature6, feature7, \
                   feature2, feature4, feature8]
    sortedFeatures = sorted(featureList, key=lambda x: x.startPos)
    item = None
    currentStart = 0
    currentEnd = -1
    overlapType = None
    regionID = 0
    featureList = []
    regionList = []
    finalFeature = sortedFeatures[-1]
    finalStart = finalFeature.startPos


# Method to write results from findRegion to the supplied output file.
# Output format is currently 'propietary'.
# Input: A dictionary of contigs, with contigs being the keys, pointing
# to a list of regions for each contig
# Returns: 0
def writeRegions(regions, outfile):
    with open(outfile, 'w') as outHandle:
        for key in sorted(regions.keys()):
            regionStr = ''
            for region in regions[key]:
                finalStr = ''
                regionLowGC = False
                geneCount = 0
                for feature in region.featureList:
                    if feature.lowGC == True:
                        regionLowGC = True
                        continue
                    geneCount += 1
                    finalStr += (str(feature.feature) + ':' +
                                 str(feature.startPos) + '-' + str(feature.endPos)
                                 + ':' + str(feature.strand) + ';')
                outHandle.write(str(key) + '\t' +
                                'regionFinder\tregion\t' +
                                str(region.startPos) + '\t' +
                                str(region.endPos) + '\t.\t.\t.\t' +
                                'ID=' + str(region.region) + ';' +
                                finalStr + 'tools=' +
                                str(region.tools) + ';geneCount=' +
                                str(geneCount) + ';' + 'lowGC=' +
                                str(regionLowGC) + ';toolCount=' +
                                str(len(region.tools)) + '\n')
            
    return(0)
    
# Method to identify regions based on a list of features from each
# contig.
# Input: a list of features for an individual contig -> features
# Returns: a list of regions of overlapping features for the contig
# being considered -> regionList
def findRegion(features):
    item = None
    currentStart = 0
    currentEnd = -1
    overlapType = None
    regionID = 0
    featureList = []
    sortedFeatures = sorted(features, key=lambda x: x.startPos)
    regionList = []
    programSet = set([])
    firstFeature = sortedFeatures[0]
    finalFeature = sortedFeatures[-1]
    firstStart = firstFeature.startPos
    finalStart = finalFeature.startPos
    first = True

    # Iterate over sorted features for one contig, and identify
    # regions. Lots of unnecessary overlap type detection for the
    # approach that I am taking. Can certainly be simplified later...
    for item in sortedFeatures:
        # a check to start the region identification process, if
        # the end feature start position has not been ecnountered.
        if (currentStart == 0):
            currentStart = item.startPos
            currentEnd = item.endPos
            featureList.append(item)
            programSet.add(item.tool)
        elif (item.startPos == finalStart):
            if (item == finalFeature):
                if first == True:
                    regionID += 1
                currentEnd = item.endPos
                featureList.append(item)
                programSet.add(item.tool)
                regionList.append(region(regionID, currentStart, currentEnd, featureList, \
                                         programSet))
            elif (item.startPos > currentEnd):
                first = False
                regionID += 1
                regionList.append(region(regionID, currentStart, currentEnd, featureList, \
                                         programSet))
                regionID += 1
                featureList = []
                programSet = set([])
                featureList.append(item)
                programSet.add(item.tool)
                currentStart = item.startPos
                currentEnd = item.endPos
            else:
                featureList.append(item)
        # check to see if a new region should be started.
        elif (item.startPos > currentEnd):
            if item.feature == 'gene-TRIVIDRAFT_163286':
                print('it closed the region')
                for feature in featureList:
                    feature.printFeature()
            regionID += 1
            newRegion = region(regionID, currentStart, currentEnd, featureList, programSet)
            currentStart = item.startPos
            currentEnd = item.endPos
            featureList = []
            programSet = set([])
            featureList.append(item)
            programSet.add(str(item.tool))
            regionList.append(newRegion)
        else:
            if item.endPos > currentEnd:
                currentEnd = item.endPos
            featureList.append(item)
            programSet.add(item.tool)
    return(regionList)


# Method to parse supplied GFF files to create a list of features for each
# contig. Each contig ID is a key used to access a list of features
# from a dictionary.
# Input: a list of GFF files -> gffFiles
# Returns: a dictionary with contig IDs as keys, with each
# key pointing to a list of features -> contigs
def parseGFF(gffFiles):
    contigs = {}
    for currentGFF in gffFiles:
        gffHandle = open(currentGFF)
        print(currentGFF)
        # pull record from GFF. each record is high level, referring
        # to the parent gene. program can be expanded to include
        # features such as start and stop codons as well as mRNA/CDS
        # sequences for each record.
        for rec in GFF.parse(gffHandle):
            features = []
            lowGC = False
            if rec.id not in contigs:
                contigs[rec.id] = []
            for item in rec.features:
                if item.qualifiers['source'][0] == 'isochore':
                    lowGC = True
                else:
                    lowGC = False
                features.append(feature(item.id, rec.id, item.location.start, \
                                        item.location.end, item.location.strand, \
                                        item.qualifiers['source'][0], lowGC))
            contigs[rec.id] += features
    gffHandle.close()
    return(contigs)

# Main function
def main(args):
    # parse GFF features provided by input from command line
    contigFeatures = parseGFF(args.gffFiles)
    contigRegions = {}

    # iterate over each contig parsed by parseGFF, and determine a list
    # of regions for each contig's features
    for contig in contigFeatures.keys():
        finalRegions = []
        if contig not in contigRegions:
            contigRegions[contig] = []
        finalRegions += (findRegion(contigFeatures[contig]))
        contigRegions[contig] = finalRegions

    #print(args.outFile)
    writeRegions(contigRegions, args.outFile)
        
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse and identify overlapping regions from GFF files. Contig names MUST match.')
    parser.add_argument('-f', type=str, dest='gffFiles', nargs='+', help='List of spaced GFF files for analysis', required=True)
    parser.add_argument('-o', type=str, dest='outFile', help='Name of desired output file.', required=True)
    args = parser.parse_args()
    main(args)
    
