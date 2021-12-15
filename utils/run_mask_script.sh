mask_vcf=../data/to_mask/problematic_sites_sarsCov2.v8.vcf
reference_id=NC_045512.2
#alignment_dir=../data/alignments
alignment_dir=../data/alignments/human_animal_subsets
#alignment_dir=../data/alignments/mini_animal_trees
#alignment_dir=../data/alignments/matched
 alignment_dir=../data/alignments/mink_paper
#input=../data/alignments/deer_n76_201121.aln
#output=../data/alignments/deer_n76_201121.v8_masked.aln
#
#python mask_alignment_using_vcf.py \
#	-i ${input} \
#	-o ${output} \
#	-v ${mask_vcf} \
#	-r ${reference_id} \
#	--both

for input in ${alignment_dir}/*audacity_only.aln
do
	output=$(echo $input|sed 's/audacity_only.aln/audacity_only.v8_masked.aln.fasta/g')
	echo $output

	python mask_alignment_using_vcf.py \
   		-i ${input} \
   		-o ${output} \
   		-v ${mask_vcf} \
   		-r ${reference_id} \
   		--both
done
