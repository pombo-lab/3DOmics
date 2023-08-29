# MetaGAM QC analysis pipeline
This pipeline was developed in Ana Pombo lab as part of work for 3DOmics consortia. It takes any number of fastq raw sequencing files and compute the reads and coverage for different bacterial species, as well as human and chicken contaminations to determine quality of DNA extraction and quality metrics.

It requires the following software package to be installed: bedtools, samtools, bowtie2, python, python pandas and python numpy libraries.
This can be done with the following conda command:

conda create -n 3domics_pombo -c bioconda bedtools samtools bowtie2 python pandas numpy

Usage example:
$./metagenome.pipeline.sh input_folder name_of_the_dataset

The pre-formatted fasta metagenome containing human, chicken and bacterial genomes, as well as matching bowties indexes are available for download from the links below:

http://bimsbstatic.mdc-berlin.de/pombo/sasha/3DOmics_pipeline/metagenome.fasta.zip

http://bimsbstatic.mdc-berlin.de/pombo/sasha/3DOmics_pipeline/metagenome.indices.zip
