awk '{print $1, substr($2,4,8) "_1" substr($7,10,11)};/#contigs_kept/{exit}' Cmultiflora.curated.fid.reassignments.tsv > Cmultiflora.curated.haplotigs.renamed.tsv
# delete comment lines from Cmultiflora.curated.haplotigs.renamed.tsv
awk 'NR==FNR{ a[$1]=$2; next } /^>/{ id=a[substr($0, 2)]; if (id!=""){ print ">" id; next } } 1' Cmultiflora.curated.haplotigs.renamed.tsv Cmultiflora.curated.haplotigs.fna > Cmultiflora.curated.haplotigs.renamed.fna
sed 's/^>CMU\(.*\)/>\1/' Cmultiflora.haplotigs.fna Cmultiflora.curated.haplotigs.renamed.fna > Cmultiflora.curated.haplotigs.merged.fid.fna
