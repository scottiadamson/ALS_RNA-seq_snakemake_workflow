#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=star_index_%j.log


#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#unzip hg38.fa.gz
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gff3.gz
#unzip gencode.v43.basic.annotation.gff3.gz

ref_dir="/gpfs/commons/groups/nygcfaculty/knowles_phatnani/ALS_compartment_RNA_seq/references"

module load star/2.7.10b
star --runThreadN 8 --runMode genomeGenerate --genomeDir $ref_dir/star_index --genomeFastaFiles $ref_dir/hg38_canonical_only.fa --sjdbGTFfile $ref_dir/gencode.v43.basic.annotation.gff3 

