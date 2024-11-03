import sys, re
import gffutils


print("GFF file: {}".format(sys.argv[1]))
print("List of gene IDs: {}".format(sys.argv[2]))
print("Output: {}".format(sys.argv[3]))

inGFF = sys.argv[1]
geneList = sys.argv[2]
outGFF = sys.argv[3]

print("Parsing & filtering GFF...")
db = gffutils.create_db(inGFF, ':memory:', force=True, merge_strategy="create_unique", keep_order=True)

print("Writing output...")
with open(geneList, "r") as genes, open(outGFF, "w") as file:
    for geneID in genes:
        gene = db[geneID.strip()]
        file.write(gene.__str__() + "\n")

        for mRNA in db.children(gene, featuretype='mRNA', order_by='start'):
            file.write(mRNA.__str__() + "\n")

            for feature in db.children(mRNA, order_by='start'):
                file.write(feature.__str__() + "\n")

        file.write("\n")
