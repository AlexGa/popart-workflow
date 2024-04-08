def get_raw_fq(wildcards):
    u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    return([ f"{u.fq1}", f"{u.fq2}" ])

def get_raw_reads(sample, unit, fq):
    files = [units.loc[sample].loc[unit, fq]]
    return files

def get_trim_pipe_input(wildcards):
    return get_raw_reads(wildcards.sample, wildcards.unit, wildcards.fq)

def get_trim_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    return expand(
    "pipe/trim/{S}/{U}.{{read}}.fq{E}".format(
            S=wildcards.sample, U=wildcards.unit, E=ending
        ),
        read=["fq1", "fq2"],
    )

rule trim_pipe:
    input:
        get_trim_pipe_input,
    output:
        pipe("pipe/trim/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/trim/{sample}-{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fq|fq\.gz",
    threads: 0 
    shell:
        "cat {input} > {output} 2> {log}"

rule cutadapt_pe:
    input: get_trim_input
    output:
        fastq1="results/trimmed/{sample}-{unit}/{sample}_{unit}_R1.fq.gz",
        fastq2="results/trimmed/{sample}-{unit}/{sample}_{unit}_R2.fq.gz",
        qc="results/trimmed/{sample}-{unit}/{sample}_{unit}.paired.qc.txt"
    log:
        "logs/cutadapt/{sample}-{unit}/{sample}_{unit}.log",
    params:
        others=config["params"]["cutadapt-pe"],
        # adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 8
    wrapper:
        "file:popart-workflow/wrapper/cutadapt_v0.30.0/pe"