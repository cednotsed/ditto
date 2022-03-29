threads=10
#alignment_dir=../data/alignments/all_animals
alignment_dir=../data/alignments/human_animal_subsets/V1
steps=audacity_only.v8_masked

for input in ${alignment_dir}/*${steps}.fasta
do
	output=$(echo $input| sed 's/.fasta//g'| sed 's/alignments/trees/g')
	echo $output

	iqtree2 \
		-nt ${threads} \
		-s ${input} \
		-m GTR+G \
		--prefix ${output} \
		-B 1000 \
		-alrt 1000 \
		-bnni
done
