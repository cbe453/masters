#!/User/bin/env python
# Author: Connor Burbridge
# Affiliation: University of Saskatchewan, Computer Science.
# Supervisors: Dave Schneider and Tony Kusalik.
# Purpose: program to identify 'regions' of overlapping features from one or multiple GFF
# files for MSc.

import argparse
#from BCBio.GFF import GFFExaminer
from BCBio import GFF


# Class definition for a region. Will likely need to be expanded on
# for future analysis
class region:
    def __init__(self, regionID, start, end, features):
        self.region = regionID
        self.startPos = start
        self.endPos = end
        self.featureList = features

    def printRegion(self):
        print(str(self.region) + '\t' + str(self.startPos) + '\t' + str(self.endPos))
        #print(len(self.featureList))
        for item in self.featureList:
            #print('\t' + item.feature + '\t' + str(item.startPos) + '\t' + str(item.endPos))
            item.printFeature()

            
# Class definition for a feature from a GFF file. Will likely need to
# be expanded on for future analysis
class feature:
    def __init__(self, featureID, contig, start, end, strand):
        self.feature = featureID
        self.contig = contig
        self.startPos = start
        self.endPos = end
        self.strand = strand
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
    print('\t' + 'Feature count: ' + str(len(sortedFeatures)))
    regionList = []
    firstFeature = sortedFeatures[0]
    finalFeature = sortedFeatures[-1]
    firstStart = firstFeature.startPos
    finalStart = finalFeature.startPos

    # iterate over sorted features for one contig, and identify
    # regions. lots of unnecessary overlap type detection for the
    # approach that I am taking. can certainly be simplified later...
    for item in sortedFeatures:
        # end feature condition. completes the final region for a contig
        if item.startPos == finalStart:
            regionID += 1
            currentStart = item.startPos
            currentEnd = item.endPos
            featureList.append(item)
            newRegion = region(regionID, currentStart, currentEnd, featureList)
            regionList.append(newRegion)
        # a check to start the region identification process, if
        # the end feature start position has not been ecnountered.
        elif (currentStart == 0):
            currentStart = item.startPos
            currentEnd = item.endPos
            featureList.append(item)
        # check to see if a new region should be started relative to the 
        # previous region
        elif (item.startPos > currentEnd):
            regionID += 1
            newRegion = region(regionID, currentStart, currentEnd, featureList)
            currentStart = item.startPos
            currentEnd = item.endPos
            featureList = []
            featureList.append(item)
            regionList.append(newRegion)
        # check if complete overlap. not really necessary at this point...
        elif (item.startPos == currentStart) and (item.endPos == currentEnd):
            overlapType = 'complete overlap'
            item.overlap = overlapType
            featureList.append(item)
        # if start position of feature is within range of currentStart
        # and currentEnd.
        elif (item.startPos >= currentStart) and (item.startPos <= currentEnd):
            if item.endPos > currentEnd:
                currentEnd = item.endPos
                overlapType = 'right overlap'
                item.overlap = overlapType
            elif item.endPos <= currentEnd:
                overlapType = 'within'
                item.overlap = overlapType
            featureList.append(item)
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
        # pull record from GFF. each record is high level, referring
        # to the parent gene. program can be expanded to include
        # features such as start and stop codons as well as mRNA/CDS
        # sequences for each record.
        for rec in GFF.parse(gffHandle):
            features = []
            if rec.id not in contigs:
                contigs[rec.id] = []
            for item in rec.features:
                features.append(feature(item.id, rec.id, item.location.start, \
                                        item.location.end, item.location.strand))
            contigs[rec.id] += features
    gffHandle.close()
    return(contigs)

# Main function
def main(args):
    print(args.gffFiles)
    # parse GFF features provided by input from command line
    contigFeatures = parseGFF(args.gffFiles)
    contigRegions = {}

    # iterate over each contig parsed by parseGFF, and determine a list
    # of regions for each contig's features
    for contig in contigFeatures.keys():
        print(contig)
        finalRegions = []
        if contig not in contigRegions:
            contigRegions[contig] = []
        finalRegions += (findRegion(contigFeatures[contig]))
        contigRegions[contig] = finalRegions
        print('\t' + 'Region count: ' + str(len(contigRegions[contig])))
            
    return(0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse and identify overlapping regions from GFF files. Contig names MUST match.')
    parser.add_argument('-f', type=str, dest='gffFiles', nargs='+', help='List of spaced GFF files for analysis')
    args = parser.parse_args()
    main(args)
    
