alignment_dir=../data/alignments/human_animal_subsets

for input in ${alignment_dir}/*v8_masked.aln.unambiguous.fasta
do
	output=$(echo $input|sed 's/v8_masked.aln.unambiguous.fasta/v8_masked.aln.unambiguous.dedup.fasta/g')
	echo $output
	seqkit rmdup -s < ${input} > ${output}

done
