alignment_dir=../data/alignments/human_animal_subsets/V5
alignment_dir=../data/alignments/temporal

for input in ${alignment_dir}/*v8_masked.unambiguous.fasta
do
	output=$(echo $input|sed 's/v8_masked.unambiguous.fasta/v8_masked.unambiguous.dedup.fasta/g')
	echo $output
	seqkit rmdup -s < ${input} > ${output}

done
