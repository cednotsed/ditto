aln_path=data/alignments/human_animal_subsets/V5/
for i in $(ls data/alignments/human_animal_subsets/V5/*.dedup.fasta)
do
	prefix=$(echo $i|sed "s|$aln_path||g"|sed "s|.audacity_only.v8_masked.unambiguous.dedup.fasta||g")
	echo $prefix
done
