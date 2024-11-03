import sys, re
import gffutils


inGFF = sys.argv[1]
outBED = sys.argv[2]

db = gffutils.create_db(inGFF, ':memory:', force=True, keep_order=True, verbose=True, merge_strategy="create_unique")

with open(outBED, "w") as file:
    for feature in db.create_introns():
        file.write("{}\t{}\t{}\t{}\t0\t{}\n".format(feature.seqid, feature.start-1, feature.stop-1, feature.attributes['Parent'][0], feature.strand))

# print("Updating DB...")
# introns = list(db.create_introns())
# def intron_id(f):
#     return '_'.join(f['ID'])
#
# db.update(
#     introns,
#     id_spec={'intron': [intron_id]},
#     verbose=True,
#     merge_strategy="create_unique"
# )
