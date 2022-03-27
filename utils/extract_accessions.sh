for i in accessions/hiv/*accessions_only.csv
do
	echo $i
	output=$(echo $i| sed 's/.accessions_only.csv/.audacity_only.v8_masked.fasta/g' | sed 's/accessions/aln/g')
	echo $output
	seqtk subseq mmsa_2021-12-08/2021-12-08_masked.fa ${i} > ${output}
done
