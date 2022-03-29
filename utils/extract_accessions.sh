mmsa_date=2022-03-26
prefix=V6
prefix=all_animals

for i in accessions/$prefix/*accessions_only.csv
do
	echo $i
	output=$(echo $i| sed 's/.accessions_only.csv/.audacity_only.v8_masked.fasta/g' | sed 's/accessions/aln/g')
	echo $output
	seqtk subseq mmsa_$mmsa_date/${mmsa_date}_masked.fa ${i} > ${output}
done
