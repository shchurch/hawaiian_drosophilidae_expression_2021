# SHC 2021
# These commands were used to annotate the assembled transcriptomes 
#  using local BLAST against custom Drosophila databases 

#	### BUILD DATABASES
#		# nt
#		update_blastdb.pl nt --decompress
#		
#		# nr
#		update_blastdb.pl nr --decompress
#		
#		# Dmel transcripts
#		makeblastdb -in ../Hawaiian_Sequences/Other_Drosophila_Genomes_Sep2018/dmel-all-transcript-r6.29.fasta  -dbtype nucl -title dmel-transcript -out dmel-transcript
#		
#		# Dmel translation
#		makeblastdb -in ../Hawaiian_Sequences/Other_Drosophila_Genomes_Sep2018/dmel-all-translation-r6.29.fasta  -dbtype prot -title dmel-translation -out dmel-translation
#		
#		# Dmel Dvir and Dgri transcript
#		cd ../Hawaiian_Sequences/Other_Drosophila_Genomes_Sep2018/
#		cat dvir-all-transcript-r1.07.fasta dgri-all-transcript-r1.05.fasta dmel-all-transcript-r6.29.fasta > concatenated_transcripts_dmel_dgri_dvir.fasta
#		cd ../../blast_trinity_Feb2021
#		
#		makeblastdb -in ../Hawaiian_Sequences/Other_Drosophila_Genomes_Sep2018/concatenated_transcripts_dmel_dgri_dvir.fasta  -dbtype nucl -title dros-transcript -out dros-transcript
#		
#		# Dmel Dvir and Dgri translation
#		cd ../Hawaiian_Sequences/Other_Drosophila_Genomes_Sep2018/
#		cat dvir-all-translation-r1.07.fasta dgri-all-translation-r1.05.fasta dmel-all-translation-r6.29.fasta > concatenated_translation_dmel_dgri_dvir.fasta
#		cd ../../blast_trinity_Feb2021
#		
#		makeblastdb -in ../Hawaiian_Sequences/Other_Drosophila_Genomes_Sep2018/concatenated_translation_dmel_dgri_dvir.fasta  -dbtype prot -title dros-translation -out dros-translation


# SPLIT TRANSCRIPTOMES

 ASSEMBLY=$'008D'
# ASSEMBLY=$'106A'
# ASSEMBLY=$'025A'
# ASSEMBLY=$'055A'
# ASSEMBLY=$'040C'
# ASSEMBLY=$'16_1'
# ASSEMBLY=$'002D'
# ASSEMBLY=$'CFB'
# ASSEMBLY=$'088B'
# ASSEMBLY=$'020A'
# ASSEMBLY=$'029A'
# ASSEMBLY=$'043D'

# cp ../assembly_phylogeny_expression_Hawaiian_Drosophila_June2021/agalma/reports/${ASSEMBLY}/*.assembly.fa ${ASSEMBLY}_assembly.fa

awk -v var=${ASSEMBLY} 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf(var"_%d.fa",++i);} print >> file; n_seq++; next;} { print >> file; }' < ${ASSEMBLY}_assembly.fa

NUM_FILES=$(ls ${ASSEMBLY}_[0-9]*fa | wc -l)


# BLAST

### blastn against dmel transcripts
touch script_${ASSEMBLY}_dmel-transcript_array.sh
echo '#!/bin/bash' >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -J dmel-transcript " >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -p shared,unrestricted,serial_requeue" >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -n 4 " >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -N 1 " >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -t 1-00:00 " >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH --mem 12G " >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -o dmel-transcript%A_%a.out # Standard output" >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "#SBATCH -e dmel-transcript%A_%a.err # Standard error" >>script_${ASSEMBLY}_dmel-transcript_array.sh
echo "blastn -evalue 1e-3 -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -db dmel-transcript -max_target_seqs 1 -query "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}.fa -out "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}_dmel-transcript.tsv -num_threads 4" >>script_${ASSEMBLY}_dmel-transcript_array.sh

sbatch --array=1-$NUM_FILES script_${ASSEMBLY}_dmel-transcript_array.sh

### blastn against dmel translated proteins
touch script_${ASSEMBLY}_dmel-translation_array.sh
echo '#!/bin/bash' >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -J dmel-translation " >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -p shared,unrestricted,serial_requeue" >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -n 4 " >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -N 1 " >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -t 1-00:00 " >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH --mem 12G " >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -o dmel-translation%A_%a.out # Standard output" >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "#SBATCH -e dmel-translation%A_%a.err # Standard error" >>script_${ASSEMBLY}_dmel-translation_array.sh
echo "blastx -evalue 1e-3 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -db dmel-translation -max_target_seqs 1 -query "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}.fa -out "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}_dmel-translation.tsv -num_threads 4" >>script_${ASSEMBLY}_dmel-translation_array.sh

sbatch --array=1-$NUM_FILES script_${ASSEMBLY}_dmel-translation_array.sh

### blastn against dmel, dvir, and dgri transcripts
touch script_${ASSEMBLY}_dros-transcript_array.sh
echo '#!/bin/bash' >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -J dros-transcript" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -p shared,unrestricted,serial_requeue" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -n 4" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -N 1" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -t 1-00:00" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH --mem 12G" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -o dros-transcript%A_%a.out # Standard output" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "#SBATCH -e dros-transcript%A_%a.err # Standard error" >>script_${ASSEMBLY}_dros-transcript_array.sh
echo "blastn -evalue 1e-3 -word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -db dros-transcript -max_target_seqs 1 -query "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}.fa -out "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}dros-transcript.tsv -num_threads 4" >>script_${ASSEMBLY}_dros-transcript_array.sh

sbatch --array=1-$NUM_FILES script_${ASSEMBLY}_dros-transcript_array.sh

### blastn against dmel, dvir, and dgri translation
touch script_${ASSEMBLY}_dros-translation_array.sh
echo '#!/bin/bash' >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -J dros-translation" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -p shared,unrestricted,serial_requeue" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -n 4" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -N 1" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -t 1-00:00" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH --mem 12G" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -o dros-translation%A_%a.out # Standard output" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "#SBATCH -e dros-translation%A_%a.err # Standard error" >>script_${ASSEMBLY}_dros-translation_array.sh
echo "blastx -evalue 1e-3 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' -db dros-translation -max_target_seqs 1 -query "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}.fa -out "${ASSEMBLY}"_\${SLURM_ARRAY_TASK_ID}_dros-translation.tsv -num_threads 4" >>script_${ASSEMBLY}_dros-translation_array.sh

sbatch --array=1-$NUM_FILES script_${ASSEMBLY}_dros-translation_array.sh

### move files to out directory
	mkdir ${ASSEMBLY}_dmel-transcript_out
	mv ${ASSEMBLY}_*dmel-transcript* ${ASSEMBLY}_dmel-transcript_out/
	cd ${ASSEMBLY}_dmel-transcript_out
	cat ${ASSEMBLY}*tsv >> ${ASSEMBLY}_dmel-transcript.tsv
	cd ..
	
	mkdir ${ASSEMBLY}_dmel-translation_out
	mv ${ASSEMBLY}_*dmel-translation* ${ASSEMBLY}_dmel-translation_out/
	cd ${ASSEMBLY}_dmel-translation_out
	cat ${ASSEMBLY}*tsv >> ${ASSEMBLY}_dmel-translation.tsv
	cd ..
	
	mkdir ${ASSEMBLY}_dros-transcript_out
	mv ${ASSEMBLY}_*dros-transcript* ${ASSEMBLY}_dros-transcript_out/
	cd ${ASSEMBLY}_dros-transcript_out
	cat ${ASSEMBLY}*tsv >> ${ASSEMBLY}_dros-transcript.tsv
	cd ..
	
	mkdir ${ASSEMBLY}_dros-translation_out
	mv ${ASSEMBLY}_*dros-translation* ${ASSEMBLY}_dros-translation_out/
	cd ${ASSEMBLY}_dros-translation_out
	cat ${ASSEMBLY}*tsv >> ${ASSEMBLY}_dros-translation.tsv
	cd ..
	
	mv ${ASSEMBLY}_[0-9]*fa ${ASSEMBLY}_dros-translation_out/