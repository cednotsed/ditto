threads=10
#alignment_dir=../data/alignments
#alignment_dir=../data/alignments/human_animal_subsets/V5
alignment_dir=../data/alignments/temporal
#alignment_dir=../data/alignments/all_animals
#alignment_dir=../data/alignments/mini_animal_trees
#alignment_dir=../data/alignments/matched
#alignment_dir=../data/alignments/mink_paper
steps=audacity_only.v8_masked.unambiguous.dedup
#steps=audacity_only.v8_masked

for input in ${alignment_dir}/*${steps}.fasta
do
	output=$(echo $input| sed 's/.fasta/.tree/g'| sed 's/alignments/trees/g')
	echo $output

	augur tree \
		--nthreads ${threads} \
		--alignment ${input} \
		--method iqtree \
		--substitution-model GTR+G \
		--output ${output}
#		--tree-builder-args "-mem 16G"
done
