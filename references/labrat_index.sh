#!/bin/bash
#SBATCH --job-name=labrat_index
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=labrat_index_%j.log

source /gpfs/commons/home/sadamson/.bashrc
mamba activate /gpfs/commons/home/sadamson/dev/mambaforge/envs/labrat_env

ref_dir="/gpfs/commons/groups/nygcfaculty/knowles_phatnani/ALS_compartment_RNA_seq/references"

LABRAT.py --mode makeTFfasta --gff $ref_dir/gencode.v43.basic.annotation.gff3 --genomefasta $ref_dir/hg38.fa --lasttwoexons --librarytype RNAseq
