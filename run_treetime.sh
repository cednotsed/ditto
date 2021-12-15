prefix=mink.Denmark.n10512
prefix=mink.Netherlands.n3750
prefix=mink.USA.n35777
aln_file=data/alignments/human_animal_subsets/${prefix}.audacity_only.v8_masked.aln.unambiguous.dedup.fasta
tree_file=data/trees/human_animal_subsets/${prefix}.audacity_only.v8_masked.aln.unambiguous.dedup.tree
date_file=data/metadata/human_animal_subsets/${prefix}.unambiguous.dedup.dates_only.csv
res_dir=results/human_animal_subsets/dating_out/${prefix}.unambiguous.dedup

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
--outdir ${res_dir}

treetime clock \
--tree ${tree_file} \
--dates ${date_file} \
--clock-filter 3 \
--aln ${aln_file} \
--covariation \
--outdir ${res_dir}

