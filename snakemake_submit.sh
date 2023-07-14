#!/bin/bash
#SBATCH --job-name=ALS_RNA-seq_process
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=5G
#SBATCH --time=48:00:00
#SBATCH --output=ALS_RNA-seq_process_%j.log

source ~/.bashrc
mamba activate ALS_analysis_env

base_dir="/gpfs/commons/groups/nygcfaculty/knowles_phatnani/ALS_compartment_RNA_seq"
slurm_config_folder="/gpfs/commons/groups/nygcfaculty/knowles_phatnani/ALS_compartment_RNA_seq/config/slurm_config"
snakemake --snakefile $base_dir/Snakefile --profile $slurm_config_folder --nolock -k

