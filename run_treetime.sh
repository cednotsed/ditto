aln_path=data/alignments/human_animal_subsets/V5/

for i in $(ls data/alignments/human_animal_subsets/V5/*.dedup.fasta)
do
	prefix=$(echo $i|sed "s|$aln_path||g"|sed "s|.audacity_only.v8_masked.unambiguous.dedup.fasta||g")
	echo $prefix

	aln_file=data/alignments/human_animal_subsets/V5/${prefix}.audacity_only.v8_masked.unambiguous.dedup.fasta
	tree_file=data/trees/human_animal_subsets/V5/${prefix}.audacity_only.v8_masked.unambiguous.dedup.tree
	date_file=data/metadata/human_animal_subsets/V5/${prefix}.unambiguous.dedup.dates_only.csv
	res_dir=results/human_animal_subsets/V5/dating_out/${prefix}.unambiguous.dedup

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
