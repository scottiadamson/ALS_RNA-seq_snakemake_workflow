#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --output=kallisto_index_%j.log

source ~/.bashrc
mamba activate ALS_analysis_env

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz
#gunzip gencode.v43.transcripts.fa.gz

base_dir="/gpfs/commons/groups/nygcfaculty/knowles_phatnani/ALS_compartment_RNA_seq/references"
kallisto index -i $base_dir/kallisto.index $base_dir/gencode.v43.transcripts.fa 

