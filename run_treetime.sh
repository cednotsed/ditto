aln_path=data/alignments/human_animal_subsets/V5/
aln_path=data/alignments/temporal

for aln in $(ls $aln_path/*.dedup.fasta)
do
	prefix=$(echo $aln|sed "s|$aln_path||g"|sed "s|.audacity_only.v8_masked.unambiguous.dedup.fasta||g")
	#echo $prefix

	aln_file=$aln
	tree_file=$(echo $aln|sed "s|alignments|trees|g"| sed "s|.dedup.fasta|.dedup.tree|g")
	date_file=$(echo $aln|sed "s|alignments|metadata|g"| sed "s|.dedup.fasta|.dedup.dates_only.csv|g"| sed "s|audacity_only.v8_masked.||g")
	echo $tree_file
	echo $date_file
	echo $aln_file
	res_dir=results/mutation_rates/${prefix}.unambiguous.dedup

	echo $aln_file
	echo $tree_file
	echo $date_file

	treetime \
	--tree ${tree_file} \
	--aln ${aln_file} \
	--dates ${date_file} \
	--clock-filter 3 \
	--covariation \
	--confidence \
	--coalescent skyline \
	--relax 1.0 0 \
	--outdir ${res_dir}

#treetime clock \
#--tree ${tree_file} \
#--dates ${date_file} \
#--clock-filter 3 \
#--aln ${aln_file} \
#--relax 1.0 0 \
#--outdir ${res_dir}

done
