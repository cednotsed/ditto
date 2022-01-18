reference=../data/genomes/Wuhan-Hu-1_NC_045512.2.fasta
threads=10
#genome_dir=../data/genomes
genome_dir=../data/genomes/human_animal_subsets/V2
#genome_dir=../data/genomes/mini_animal_trees
#genome_dir=../data/genomes/matched
#genome_dir=../data/genomes/mink_paper

#genomes=../data/genomes/deer_n76_201121.fasta
#alignment_out=../data/alignments/deer_n76_201121.aln
#
#augur align \
#	--nthreads ${threads} \
#	--sequences ${genomes} \
#	--output ${alignment_out} \
#	--reference-sequence ${reference} \
#	--method mafft

for genomes in ${genome_dir}/mink.*audacity_only.fasta
do
	echo $genomes
	alignment_out=$(echo $genomes|sed 's/.fasta/.aln/g'| sed 's/genomes/alignments/g')
	augur align \
	   --nthreads ${threads} \
	   --sequences ${genomes} \
	   --output ${alignment_out} \
	   --reference-sequence ${reference} \
	   --method mafft
done

