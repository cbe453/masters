import plotly.graph_objects as go
import argparse
from BCBio import GFF

def agree(region):
    numFeatures = int(region.qualifiers['geneCount'][0])
    initialStart = 0
    initialStop = 0
    for i in range(1, numFeatures+1):
        key = 'ftr' + str(i)
        #print(region.qualifiers[key])
        feature = region.qualifiers[key][0]
        if initialStart == 0:
            initialStart = int(feature.split(':')[1].split('-')[0])
            initialStop = int(feature.split(':')[1].split('-')[1])
            continue
        curStart = int(feature.split(':')[1].split('-')[0])
        curStop = int(feature.split(':')[1].split('-')[1])

        if (initialStart != curStart) or (initialStop != curStop):
            return(False)

    return(True)
        #print('start: ' + str(start))
        #print('stop: ' + str(stop))
        
def parseGFF(gffFile):
    contigs = {}
    gffHandle = open(gffFile)

    single = 0
    partialSupport = 0
    partialSupportAgree = 0
    partialSupportDisagree = 0
    fullSupport = 0
    fullSupportAgree = 0
    fullSupportDisagree = 0
    uncertainRegions = 0
    totalRegions = 0
    
    for contig in GFF.parse(gffHandle):
        for region in contig.features:
            #print(feature)
            #change for DC1 and Tsth20, no refseq
            toolCount = int(region.qualifiers['toolCount'][0])
            geneCount = int(region.qualifiers['geneCount'][0])
            totalRegions += 1
            
            if toolCount == 3:
                fullSupport += 1
                if(agree(region)):
                    fullSupportAgree += 1
                else:
                    fullSupportDisagree += 1
            if toolCount == 2:
                fullSupport += 1
                if(agree(region)):
                    fullSupportAgree += 1
                else:
                    fullSupportDisagree += 1
            elif toolCount == 1:
                single += 1

            if geneCount > 3:
                uncertainRegions += 1

    print('Total number of regions: ' + str(totalRegions))
    print('Fully supported genes: ' + str(fullSupport))
    print('Fully supported genes with agreeing models: ' + str(fullSupportAgree))
    print('Fully supported genes with disagreeing models: ' + str(fullSupportDisagree))
    print('Partially supported genes: ' + str(partialSupport))
    print('Partially supported genes with agreeing models: ' + str(partialSupportAgree))
    print('Partially supported genes with disagreeing models: ' + str(partialSupportDisagree))
    print('Singletons: ' + str(single))
    print('Regions with more than 3 predicted models: ' + str(uncertainRegions))
    

    fig = go.Figure(data=[go.Sankey(
        node = dict(
            line = dict(color = "black", width = 0.3),
            label = ["Regions: " + str(totalRegions), "Full Support: " + str(fullSupport), \
                     "Partial Support: " + str(partialSupport), \
                     "Models Agree: " + str(fullSupportAgree) , \
                     "Models Disagree: " + str(fullSupportDisagree), \
                     "Models Agree: " + str(partialSupportAgree), \
                     "Models Disagree: " + str(partialSupportDisagree),\
                     "Singletons: " + str(single)]
        ),
        link = dict(
            arrowlen = 15,
            source = [0, 0, 1, 1, 2, 2, 0], # indices correspond to labels, eg A1, A2, A1, B1, ...
            target = [1, 2, 3, 4, 5, 6, 7],
            value = [fullSupport, partialSupport, fullSupportAgree, fullSupportDisagree, \
                     partialSupportAgree, partialSupportDisagree, single]
        ))])
    
    fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
    fig.show()
    gffHandle.close()

    return
            
def main(args):
    print(args)
    parseGFF(args.gffFile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process regions and report interesting stats.')
    parser.add_argument('-f', type=str, dest='gffFile', help='GFF file containing regions', required=True)
    parser.add_argument('-o', type=str, dest='outFile', help='Name of desired output file.', required=True)
    args = parser.parse_args()
    main(args)
