# snakemake pipeline for mapping and pre-processing ALS RNA-seq #

Scott Adamson  
David Knowles / Tuuli Lappalainen labs

### Installation: ###
First install this snakemake pipeline:
```
git clone git@github.com:scottiadamson/ALS_RNA-seq_snakemake_workflow.git
cd ALS_RNA-seq_snakemake_workflow.git 
```
Then create a mamba environment with dependencies:
```
mamba create -n ALS_analysis_env
mamba activate ALS_analysis_env
mamba install --file requirements_updated.txt
``` 

### Preparing references ###
Submit all of the reference building stuff  
Be sure to edit the paths in the .sh files to the reference folder in your snakemake directory
```
cd references
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


