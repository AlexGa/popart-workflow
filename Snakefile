from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("7.8.0")


##### setup report #####
configfile: "config/config.yaml"


##### setup singularity #####


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### load rules #####


include: "popart-workflow/rules/trim.smk"
include: "popart-workflow/rules/alignment.smk"


##### target rules #####

units = pd.read_csv(config["units"], sep="\t", dtype={"sample": str, "unit": str}).set_index(["sample", "unit"], drop=False).sort_index()

rule all:
    input:
        expand("results/bigwig/only_sars_cov2/{file.sample}-{file.unit}/{file.sample}-{file.unit}_paired.bw", file = units.itertuples())