experiment=human_animal_subsets
prefix=mink.Netherlands.n3750
prefix=mink.Denmark.n10512
dates=data/metadata/${experiment}/${prefix}.date_file.tsv
aln=data/alignments/${experiment}/${prefix}.audacity_only.v8_masked.aln.fasta
tree=data/trees/${experiment}/${prefix}.audacity_only.v8_masked.tree
outgroup=NC_045512.2

iqtree -s ${aln} \
	-T 10 \
	--date ${dates} \
	-o $outgroup \
	-te ${tree} \
	-m GTR+G \
	--date-ci 100 \
	--date-no-outgroup

