# snakemake pipeline for mapping and pre-processing ALS RNA-seq #

Scott Adamson  
David Knowles / Tuuli Lappalainen labs

### Installation: ###
First install this snakemake pipeline:
```
git clone git@github.com:scottiadamson/ALS_RNA-seq_snakemake_workflow.git
cd ALS_RNA-seq_snakemake_workflow
```
Then create a mamba environment with dependencies:
```
mamba create -n ALS_analysis_env
mamba activate ALS_analysis_env
mamba install --file requirements_updated.txt
``` 

### Preparing references ###
Download the reference files if you haven't already:
```
cd references
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz
gunzip gencode.v43.transcripts.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gff3.gz
gunzip gencode.v43.basic.annotation.gff3.gz
```
Submit all of the reference building stuff  
Be sure to edit the paths in the .sh files to the reference folder in your snakemake directory
```
sbatch kallisto_index.sh
sbatch labrat_index.sh
sbatch star_index.sh
```

### Adjust config/config.yaml ###
Edit the config.yaml file to have the paths to this snakemake file, as well as the samples of interest.  

Then run it:
```
cd ../ # only if you're still in the references folder
sbatch snakemake_submit.sh
```

### Figuring out what's going on ###
If you just want to see the steps that the snakemake pipeline is running:
```
snakemake --snakefile Snakefile --dry-run -p -k
```


