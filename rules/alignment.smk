ruleorder: coverage_counts_strand_extract > coverage_counts_strand

def get_mem_star(wildcards):
    return config["params"]["star-map-"+wildcards.organism+"-ram"] / 1e6

rule star_index_gtf:
    input:
        fasta = lambda wildcards: config["ref"]["fasta"][wildcards.organism],
        gtf = lambda wildcards: config["ref"]["gtf"][wildcards.organism]
    output:
        files = expand("references/STAR_INDEX/{{organism}}/{suf}", 
                suf = ["chrNameLength.txt", "Log.out", "chrStart.txt"]),
    log:
        "logs/STAR/{organism}/star_index.log"
    threads: 20
    params:
        ram = lambda wildcards: config["params"]["star-map-" + wildcards.organism + "-ram"],
        star_usr = lambda wildcards: config["params"]["star-index-" + wildcards.organism],
        dir = directory("references/STAR_INDEX/{organism}")
    wrapper:
        "file:popart-workflow/wrapper/star_index"

rule star_pe:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = "results/trimmed/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = "results/trimmed/{sample}-{unit}/{sample}_{unit}_R2.fq.gz",
        idx = expand("references/STAR_INDEX/{{organism}}/{suf}", suf = ["chrNameLength.txt", "Log.out", "chrStart.txt"])
    output:
        # see STAR manual for additional output files
        bam = "results/star/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bam",
        log = "results/star/{organism}/{sample}-{unit}/Log.out",
        log_final = "results/star/{organism}/{sample}-{unit}/Log.final.out",
        unmapped = expand("results/star/{{organism}}/{{sample}}-{{unit}}/Unmapped.out.mate{read}", read = ["1", "2"])
    #wildcard_constraints:
    #    organism="^(homo|sars)"
    log:
        "logs/STAR_align/{organism}/{sample}-{unit}.log"
    resources:
        mem_mb = get_mem_star
    threads: 16
    params:
        # path to STAR reference genome index
        idx="references/STAR_INDEX/{organism}/",
        # optional parameters
        extra=lambda wildcards: config["params"]["star-map-" + wildcards.organism],
        tmpdir = "results/star/{organism}/{sample}-{unit}"
    wrapper:
        "file:popart-workflow/wrapper/star_align"

rule bam_index:
    input: "results/{sample}.bam"
    output: "results/{sample}.bai"
    log:
        "logs/index/{sample}.log",
    wrapper:
        "file:popart-workflow/wrapper/samtools_index"

rule split_by_strand:
    input: "results/{subdir}/{organism}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    output: fwd = "results/{subdir}/{organism}/{sample}-{unit}/{sample}-{unit}_fwd.bam", 
            rev = "results/{subdir}/{organism}/{sample}-{unit}/{sample}-{unit}_rev.bam",
            paired = "results/{subdir}/{organism}/{sample}-{unit}/{sample}-{unit}_paired.bam"
    log:
        "log/samtools/split_by_strand/{organism}/{subdir}/{sample}-{unit}.log"
    params:
        extra = "",  # optional params string
        tmpdir = temp("results/{subdir}/{organism}/{sample}-{unit}/tmp")
    threads: 1
    wrapper:
        "file:popart-workflow/wrapper/samtools_split_strand"

rule coverage_counts_strand_extract:
    input:
        bam="results/extract/only_{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
        bai="results/extract/only_{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bai"
    output:
        bigwig="results/bigwig/only_{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bw",
    log:
        "logs/bigwig/only_{organism}/{sample}-{unit}_{strand}.log",
    params:
        extra = "--binSize 1",
        strand = "" 
    wrapper:
         "file:popart-workflow/wrapper/deeptools_coverage"

rule coverage_counts_strand:
    input:
        bam="results/star/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bam",
        bai="results/star/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bai"
    output:
        bigwig="results/bigwig/{organism}/{sample}-{unit}/{sample}-{unit}_{strand}.bw",
    log:
        "logs/bigwig/{organism}/{sample}-{unit}_{strand}.log",
    params:
        extra = "--binSize 1",
        strand = "" 
    wrapper:
         "file:popart-workflow/wrapper/deeptools_coverage"

rule fastq_sort:
    input:
        "{star_dir}/{sample}-{unit}/Unmapped.out.mate{read}"
    output:
        temp("{star_dir}/{sample}-{unit}/{sample}_{unit}_R{read}.fq")
    conda: "../envs/fastq_tools.yaml"
    shell: "fastq-sort --id {input} > {output}"

rule pack_reads:
    input:
        "{star_dir}/{sample}-{unit}/{sample}_{unit}_R{read}.fq"
    output:
        "{star_dir}/{sample}-{unit}/{sample}_{unit}_R{read}.fq.gz"
    threads: 1
    conda: "../envs/fastq_tools.yaml"
    shell: "gzip {input}"

rule extract_virus:
    input: "results/star/homo_sapiens_{suffix}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    output: bam_sorted = temp("results/star/homo_sapiens_{suffix}/{sample}-{unit}/{sample}-{unit}_aligned_sorted.bam"),
            bam = "results/extract/only_{suffix}/{sample}-{unit}/{sample}-{unit}_aligned.bam"
    log:
        "logs/extract_sars_{suffix}/{sample}/{unit}.log",
    params:
        region=lambda wildcards: config["ref"]["chromosomes"][wildcards.suffix]
    threads: 1
    wrapper:
        "file:popart-workflow/wrapper/samtools_view"

rule ratio_viral_human:
    input: all = "results/star/homo_sapiens_{suffix}/{sample}-{unit}/{sample}-{unit}_{suffix}.bam",
           viral = "results/star/only_{suffix}/{sample}-{unit}/{sample}-{unit}_{suffix}.bam",
    output: txt = "results/human_viral_ratio_{suffix}/{sample}-{unit}.txt"
    log: "logs/ratio/ratio_viral_human_{suffix}/{sample}-{unit}.log"
    threads: 4
    params:
        read_type = "paired" # single for single end read
    wrapper:
        "file:popart-workflow/wrapper/ratio_viral_human"

rule bw_tracks_norm_paired:
    input: paired = "results/extract/{organism}/{sample}/{sample}_paired.bam",
           genome_fasta = lambda wildcards: config["ref"]["fasta"][wildcards.organism]
    output: raw = "results/normed_tracks/{organism}/Raw/{sample}/{sample}_RAW_paired.bw",
            rpm = "results/normed_tracks/{organism}/RPM/{sample}/{sample}_RPM_paired.bw"
    threads: 1
    params: raw = "results/normed_tracks/{organism}/Raw/{sample}/{sample}_RAW",
            rpm = "results/normed_tracks/{organism}/RPM/{sample}/{sample}_RPM"
    log: "logs/bw_tracks_norm/{organism}/{sample}_paired.log"
    wrapper: "file:popart-workflow/wrapper/norm_tracks_sample"

rule coverage_counts_paired_CPM:
    input:
        bam="results/extract/{organism}/{sample}-{unit}/{sample}-{unit}_paired.bam",
        bai="results/extract/{organism}/{sample}-{unit}/{sample}-{unit}_paired.bai"
    output:
        bigwig="results/bigwig/{organism}/{sample}-{unit}/{sample}-{unit}_CPM_paired.bw",
    log:
        "logs/bigwig/{organism}/{sample}-{unit}_paired.log",
    params:
        extra = "--binSize 1 --normalizeUsing CPM",
        strand = "" 
    wrapper:
         "file:popart-workflow/wrapper/deeptools_coverage"
