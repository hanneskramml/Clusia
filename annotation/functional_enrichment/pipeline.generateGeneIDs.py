import sys, re
import gffutils


print("GFF file: {}".format(sys.argv[1]))
print("Prefix for gene IDs: {}".format(sys.argv[2]))
print("Output: {}".format(sys.argv[3]))

inGFF = sys.argv[1]
preGene = sys.argv[2]
outGFF = sys.argv[3]

currentChr = None
currentGene = None
currentRna = None
geneCounter = 0
rnaCounter = 0
utrCounter = 0
cdsCounter = 0
exonCounter = 0


def transform(feature):
    global currentChr, currentGene, currentRna, geneCounter, rnaCounter, utrCounter, cdsCounter, exonCounter

    if feature.attributes['ID'].__len__() > 1:
        print("Warning: Multiple IDs in feature {}".format(feature.attributes['ID']))

    if feature.seqid[3:5] != currentChr:
        currentChr = feature.seqid[3:5]
        geneCounter = 0

    if feature.featuretype == 'gene':
        rnaCounter = 0
        geneCounter += 1
        currentGene = '{}{}.g{}'.format(preGene, feature.seqid[3:5], geneCounter)  #"{:04d}".format(counter)
        feature.attributes['ID'][0] = currentGene

    elif feature.featuretype == 'mRNA':
        utrCounter = 0
        cdsCounter = 0
        exonCounter = 0
        rnaCounter += 1

        currentRna = '{}{}.g{}.t{}'.format(preGene, feature.seqid[3:5], geneCounter, rnaCounter)
        feature.attributes['ID'][0] = currentRna
        feature.attributes['Parent'][0] = currentGene

    else:
        feature.attributes['Parent'][0] = currentRna

        if feature.featuretype in ['five_prime_UTR', 'three_prime_UTR']:
            utrCounter += 1
            feature.attributes['ID'][0] = '{}.utr{}'.format(currentRna, utrCounter)

        elif feature.featuretype == 'exon':
            exonCounter += 1
            feature.attributes['ID'][0] = '{}.exon{}'.format(currentRna, exonCounter)

        elif feature.featuretype == 'CDS':
            cdsCounter += 1
            feature.attributes['ID'][0] = '{}.cds{}'.format(currentRna, cdsCounter)

        else:
            print("Warning: Unknown feature type {} in {}".format(feature.featuretype, feature.attributes['ID']))

    return feature


print("Processing...")
db = gffutils.create_db(inGFF, ':memory:', force=True, transform=transform)

print("Writing file...")
with open(outGFF, "w") as file:
    for feature in db.all_features():
        file.write(feature.__str__() + "\n")

print("Done!")
