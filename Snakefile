from os.path import join
import yaml, sys, os, csv

snakemake_dir = '/gpfs/commons/groups/nygcfaculty/knowles_phatnani/ALS_compartment_RNA_seq/'
configfile: join(snakemake_dir, 'config/config.yaml')

fastq_dir = config['fastq_dir']

def scheduler_time_hours(hour):
    return hour*60

def scheduler_memory(memory_gb):
    return memory_gb*1000

rule all:
    input:
        expand(join(snakemake_dir, "fastqc_out/{sample}_trimmed_R1_fastqc.html"), sample = config['samples']),
        expand(join(snakemake_dir, "fastqc_out/{sample}_trimmed_R2_fastqc.html"), sample = config['samples']),
        expand(join(snakemake_dir, "junction_counts/dedup/{sample}_junctions.out"), sample = config['samples']),
        expand(join(snakemake_dir, "junction_counts/normal/{sample}_junctions.out"), sample = config['samples']),
        expand(join(snakemake_dir, "kallisto_out/normal/{sample}/abundance.tsv"), sample = config['samples']),
        expand(join(snakemake_dir, "bam_stats/deduped/{sample}_stats.txt"), sample = config['samples']),
        expand(join(snakemake_dir, "bam_stats/normal/{sample}_stats.txt"), sample = config['samples']),
        expand(join(snakemake_dir, "bam_stats/deduped/{sample}.metrics.tsv"), sample = config['samples']),
        expand(join(snakemake_dir, "bam_stats/normal/{sample}.metrics.tsv"), sample = config['samples']),
        expand(join(snakemake_dir, "labrat_quant/{sample}/quant.sf"), sample = config['samples'])

ruleorder: combine_files > trim_adapter

rule combine_files:
    output:
        R1 = temp(join(snakemake_dir, "fastq/{sample}_R1.fastq.gz")),
        R2 = temp(join(snakemake_dir, "fastq/{sample}_R2.fastq.gz"))
    params:
        R1 = join(fastq_dir, "Sample_{sample}/fastq/{sample}*_001.R1.fastq.gz"), 
        R2 = join(fastq_dir, "Sample_{sample}/fastq/{sample}*_001.R2.fastq.gz")
    shell:
        """
        cat {params.R1} > {output.R1}
        cat {params.R2} > {output.R2}
        """

rule trim_adapter:
    input:
        R1 = join(snakemake_dir, "fastq/{sample}_R1.fastq.gz"),
        R2 = join(snakemake_dir, "fastq/{sample}_R2.fastq.gz")
    output:
        R1 = temp(join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz")),
        R2 = temp(join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz"))
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        cutadapt -a AGATCGGAAGAG -a CTGTCTCTTATACACATCT -A CTCTTCCGATCT -A TGTCTCTTATACACAT -m 40 -j {threads} -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        """

rule fastqc_R1:
    input:
        join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz")
    output:
        join(snakemake_dir, "fastqc_out/{sample}_trimmed_R1_fastqc.zip"),
        join(snakemake_dir, "fastqc_out/{sample}_trimmed_R1_fastqc.html")
    params:
        out_dir = join(snakemake_dir, 'fastqc_out/')
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        fastqc {input} -t {threads} -o {params.out_dir}
        """    

rule fastqc_R2:
    input:
        join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz")
    output:
        join(snakemake_dir, "fastqc_out/{sample}_trimmed_R2_fastqc.zip"),
        join(snakemake_dir, "fastqc_out/{sample}_trimmed_R2_fastqc.html")
    params:
        out_dir = join(snakemake_dir, 'fastqc_out/')
    resources:
        mem = scheduler_memory(20)
    threads:
        8
    shell:
        """
        fastqc {input} -t {threads} -o {params.out_dir}
        """    

rule star_map:
    input:
        R1 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"),
        R2 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz")
    output:
        bam_file = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam"),
        log_out = join(snakemake_dir, "bam/{sample}_Log.out"),
        final_out = join(snakemake_dir, "bam/{sample}_Log.final.out"),
        progress_out = temp(join(snakemake_dir, "bam/{sample}_Log.progress.out")),
        SJ_out = temp(join(snakemake_dir, "bam/{sample}_SJ.out.tab")),
        genome_folder = temp(directory(join(snakemake_dir, "bam/{sample}__STARgenome"))),
        pass_folder = temp(directory(join(snakemake_dir, "bam/{sample}__STARpass1")))
    threads:
        8
    resources:
        mem = scheduler_memory(41),
        time = scheduler_time_hours(8)
    message:
        "mapping with star"
    params:
        genome_index = config['star_index'],
        output_prefix = join(snakemake_dir, "bam/{sample}_"),
    shell:
        """
        STAR --genomeDir {params.genome_index} --readFilesIn {input.R1} {input.R2} --readFilesCommand zcat --runMode alignReads --outSAMattributes All --outSAMstrandField intronMotif --limitBAMsortRAM 40000000000 --outSAMtype BAM SortedByCoordinate --runThreadN {threads} --twopassMode Basic --outFileNamePrefix {params.output_prefix} 
        """

rule bam_index:
    input:
        join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam")
    output:
        temp(join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam.bai"))
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)
    message:
        "indexing bam file"
    shell:
        """    
        samtools index -@ 7 {input} 
        """

rule collate_bam:
    input:
        bam = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam"),
        bai = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam.bai")
    output:
        temp(join(snakemake_dir, "bam/{sample}_collated.bam"))
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)
    message:
        "collating bam file; first step of duplicate removal"
    shell:
        """    
        samtools collate -@ 7 -o {output} {input.bam} 
        """

rule fixmate_bam:
    input:
        join(snakemake_dir, "bam/{sample}_collated.bam")
    output:
        temp(join(snakemake_dir, "bam/{sample}_mate_fixed.bam"))
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)
    message:
        "fixing mate of bam file; second step of duplicate removal"
    shell:
        """    
        samtools fixmate -@ 7 -m {input} {output} 
        """

rule sort_fixmate_bam:
    input:
        join(snakemake_dir, "bam/{sample}_mate_fixed.bam")
    output:
        temp(join(snakemake_dir, "bam/{sample}_mate_fixed_sorted.bam"))
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)
    message:
        "sorting fixed mate bam file; third step of duplicate removal"
    shell:
        """    
        samtools sort -@ 7 -o {output} {input} 
        """

rule remove_duplicates:
    input:
        join(snakemake_dir, "bam/{sample}_mate_fixed_sorted.bam")
    output:
        bam = temp(join(snakemake_dir, "bam/{sample}_duplicate_removed.bam")),
        stats = join(snakemake_dir, "dedup_stats/{sample}_stats.txt")
    threads:
        8
    resources:
        mem = scheduler_memory(40),
        time = scheduler_time_hours(1)
    message:
        "sorting fixed mate bam file; third step of duplicate removal"
    shell:
        """    
        samtools markdup -r -f {output.stats} -@ 7 {input} {output.bam}
        """

rule index_deduped_bam:
    input:
        join(snakemake_dir, "bam/{sample}_duplicate_removed.bam")
    output:
        temp(join(snakemake_dir, "bam/{sample}_duplicate_removed.bam.bai"))
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)
    message:
        "indexing deduped bam file"
    shell:
        """    
        samtools index -@ 7 {input} 
        """

rule stats_deduped:
    input:
        bam = join(snakemake_dir, "bam/{sample}_duplicate_removed.bam"),
        bai = join(snakemake_dir, "bam/{sample}_duplicate_removed.bam.bai") 
    output:
        join(snakemake_dir, "bam_stats/deduped/{sample}_stats.txt")
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)  
    message:
        "getting bam stats for deduped bam"
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output}
        """

rule stats_normal:
    input:
        bam = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam"),
        bai = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam.bai")
    output:
        join(snakemake_dir, "bam_stats/normal/{sample}_stats.txt")
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(1)  
    message:
        "getting bam stats for normal bam"
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output}
        """

rule kallisto_quant:
    input:
        R1 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"),
        R2 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz")
    output:
        abundance_h5 = join(snakemake_dir, "kallisto_out/normal/{sample}/abundance.h5"),
        abundance_tsv = join(snakemake_dir, "kallisto_out/normal/{sample}/abundance.tsv"),
        run_info = join(snakemake_dir, "kallisto_out/normal/{sample}/run_info.json")
    params:
        kallisto_index = join(snakemake_dir, "references/kallisto.index"),
        out_dir = join(snakemake_dir, "kallisto_out/normal/{sample}")
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(6)
    message:
        "quantifying with kallisto"
    shell:
        """
        kallisto quant --bias -b 100 -i {params.kallisto_index} -t {threads} -o {params.out_dir} {input.R1} {input.R2}
        """

rule junction_count_normal:
    input:
        bam = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam"),
        bai = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam.bai")
    output:
        join(snakemake_dir, "junction_counts/normal/{sample}_junctions.out")
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(6)
    message:
        "Counting junctions"
    shell:
        """
        regtools junctions extract -a 8 -m 50 -M 500000 {input.bam} -o {output} -s XS
        """    

rule junction_count_dedup:
    input:
        bam = join(snakemake_dir, "bam/{sample}_duplicate_removed.bam"),
        bai = join(snakemake_dir, "bam/{sample}_duplicate_removed.bam.bai") 
    output:
        join(snakemake_dir, "junction_counts/dedup/{sample}_junctions.out")
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(6)
    message:
        "Counting junctions deduped"
    shell:
        """
        regtools junctions extract -a 8 -m 50 -M 500000 {input.bam} -o {output} -s XS
        """    

rule rnaseqc_normal:
    input:
        collapsed_gtf = join(snakemake_dir, "references/gencode.v43.basic.annotation_collapsed.gtf"),
        bam = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam"),
        bai = join(snakemake_dir, "bam/{sample}_Aligned.sortedByCoord.out.bam.bai")
    output:
        temp(join(snakemake_dir, "bam_stats/normal/{sample}.exon_reads.gct")),
        temp(join(snakemake_dir, "bam_stats/normal/{sample}.gene_fragments.gct")),
        temp(join(snakemake_dir, "bam_stats/normal/{sample}.gene_reads.gct")),
        temp(join(snakemake_dir, "bam_stats/normal/{sample}.gene_tpm.gct")),
        join(snakemake_dir, "bam_stats/normal/{sample}.metrics.tsv")
    params:
        out_dir = join(snakemake_dir, "bam_stats/normal/"),
        sample_id = "{sample}"
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(2)
    message:
        "running rnaseqc normal samples"
    shell:
        """
        rnaseqc -s {params.sample_id} {input.collapsed_gtf} {input.bam} {params.out_dir} 
        """    

rule rnaseqc_dedup:
    input:
        collapsed_gtf = join(snakemake_dir, "references/gencode.v43.basic.annotation_collapsed.gtf"),
        bam = join(snakemake_dir, "bam/{sample}_duplicate_removed.bam"),
        bai = join(snakemake_dir, "bam/{sample}_duplicate_removed.bam.bai") 
    output:
        temp(join(snakemake_dir, "bam_stats/deduped/{sample}.exon_reads.gct")),
        temp(join(snakemake_dir, "bam_stats/deduped/{sample}.gene_fragments.gct")),
        temp(join(snakemake_dir, "bam_stats/deduped/{sample}.gene_reads.gct")),
        temp(join(snakemake_dir, "bam_stats/deduped/{sample}.gene_tpm.gct")),
        join(snakemake_dir, "bam_stats/deduped/{sample}.metrics.tsv")
    params:
        out_dir = join(snakemake_dir, "bam_stats/deduped/"),
        sample_id = "{sample}"
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(2)
    message:
        "running rnaseqc deduped samples"
    shell:
        """
        rnaseqc -s {params.sample_id} {input.collapsed_gtf} {input.bam} {params.out_dir} 
        """   
 
rule labrat_quant:
    input:
        R1 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R1.fastq.gz"),
        R2 = join(snakemake_dir, "trimmed_fastq/{sample}_trimmed_R2.fastq.gz")
    output:
        join(snakemake_dir, "labrat_quant/{sample}/quant.sf")
    params:
        TF_ref = join(snakemake_dir, "references/TFseqs.fasta"),
        out_dir = join(snakemake_dir, "labrat_quant/{sample}")
    threads:
        8
    resources:
        mem = scheduler_memory(20),
        time = scheduler_time_hours(2)
    shell:
        """
        LABRAT.py --mode runSalmon --librarytype RNAseq --reads1 {input.R1} --reads2 {input.R2} --txfasta {params.TF_ref} --samplename {params.out_dir} --threads {threads}
        """

