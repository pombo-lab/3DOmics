#!/bin/bash
#$ -N metagenome_pipeline
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -l m_mem_free=12G
#$ -l data,os=any

"""
====================================================================================================================
Bash script name metagenome.pipeline.sh
Script written by Alexander Kukalev
Ana Pombo lab, Berlin Institute of Medical systems Biology, 
Max-Delbrueck Center for Molecular Medizine, Berlin, Germany
Date: August 25, 2023

This pipeline was developed as part of 3DOmics project. It takes any number of fastq raw sequencing files and compute the reads and coverage for different bacterial species, as well as human and chicken contaminations to determine quality of DNA extraction and quality metrics.

It requires the following software package to be installed: bedtools, samtools, bowtie2, python, python pandas and python numpy libraries.
This can be done with the following conda command:

conda create -n 3domics_pombo -c bioconda bedtools samtools bowtie2 python pandas numpy

Usage example:
$./metagenome.pipeline.sh input_folder name_of_the_dataset

=====================================================================================================================
"""


export PATH="$HOME/.guix-profile/bin${PATH:+:}:/data/pombo/Sasha/Common_files:/data/pombo/Sasha/Common_files/fastqscreen:$PATH"
export GUIX_LOCPATH="$HOME/.guix-profile/lib/locale"
export PYTHONPATH="$HOME/.guix-profile/lib/python3.9/site-packages${PYTHONPATH:+:}$PYTHONPATH"
FASTA_FILE=../metagenome_fasta_HsGgMoPaSrXr/3DOmics.metagenome.fa
BOWTIE_INDEX=../metagenome_indices_HsGgMoPaSrXr/3DOmics_indices

FOLDER=$1
NAME=$2

#=========================================================================================
# Mapping section

FASTQ_FILES="${FOLDER}*.fastq.gz"
for f in $FASTQ_FILES
do
	echo
	echo "Processing" $f
	bowtie2 -p 12 -x $BOWTIE_INDEX $f -S $f.sam && samtools view -F 4 -bS $f.sam > $f.bam #keep multimappers, remove only unmapped
	echo "Mapping successfuly done for" $f
	samtools sort $f.bam -o $f.sorted.bam && samtools index $f.sorted.bam && samtools rmdup -s $f.sorted.bam $f.rmdup.bam && samtools index $f.rmdup.bam
	echo "PCR duplicates removed for" $f.rmdup.bam
done

# Reads and coverage section
# Create bed-file with all chromosomes from fasta-file

samtools faidx $FASTA_FILE 
awk -v OFS='\t' '{print $0, "0", $2}' $FASTA_FILE.fai > ${NAME}.chrom.sizes
cut -f1,6,7 ${NAME}.chrom.sizes > ${NAME}.chromosomes.bed
cp ${NAME}.chromosomes.bed ${NAME}.reads.chromosomes.bed
cp ${NAME}.chromosomes.bed ${NAME}.coverage.chromosomes.bed

RMDUP_FILES="${FOLDER}*.rmdup.bam"
for f in $RMDUP_FILES
do
	bedtools coverage -a ${NAME}.chromosomes.bed -b $f | tee >(awk '{print $4}' > ${NAME}.reads.bed) >(awk '{print $5}' > ${NAME}.coverage.bed)
	paste ${NAME}.reads.chromosomes.bed ${NAME}.reads.bed > ${NAME}.merged.bed && mv ${NAME}.merged.bed ${NAME}.reads.chromosomes.bed
	paste ${NAME}.coverage.chromosomes.bed ${NAME}.coverage.bed > ${NAME}.merged.bed && mv ${NAME}.merged.bed ${NAME}.coverage.chromosomes.bed
done

echo -e "chrom\tstart\tstop" $RMDUP_FILES > ${NAME}.list_of_files.txt
tr ' ' \\t < ${NAME}.list_of_files.txt > ${NAME}.list_of_files.tab.txt && mv ${NAME}.list_of_files.tab.txt ${NAME}.list_of_files.txt
cat ${NAME}.list_of_files.txt ${NAME}.reads.chromosomes.bed > ${NAME}.reads.chromosomes.header.bed && mv ${NAME}.reads.chromosomes.header.bed ${NAME}.reads.all.contigs.table
cat ${NAME}.list_of_files.txt ${NAME}.coverage.chromosomes.bed > ${NAME}.coverage.chromosomes.header.bed && mv ${NAME}.coverage.chromosomes.header.bed ${NAME}.coverage.all.contigs.table
rm ${NAME}.list_of_files.txt ${NAME}.reads.bed ${NAME}.coverage.bed ${NAME}.chromosomes.bed ${NAME}.reads.chromosomes.bed ${NAME}.coverage.chromosomes.bed ${NAME}.chrom.sizes $FASTA_FILE.fai

python3 metagenome.pipeline.py ${NAME}.reads.all.contigs.table ${NAME}.coverage.all.contigs.table $NAME 


