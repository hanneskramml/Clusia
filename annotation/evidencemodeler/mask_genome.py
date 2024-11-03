import sys
import funannotate.library as lib

genome = sys.argv[1]
gaps_bed = sys.argv[2]
repeats_bed = sys.argv[3]
repeats_gff = sys.argv[4]

lib.checkMasklowMem(genome, repeats_bed, gaps_bed, 1, tmpdir='tmp')
lib.bed2gff3(repeats_bed, repeats_gff)
