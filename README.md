# Sequence analysis for RNA purification

Nora Schmidt, Sebastian Zielinski, Alexander Gabel, Jens Aydin, and Mathias Munschauer

**Identification of proteins directly bound to SARS-CoV-2 RNA in infected cells using RNA antisense purification and mass spectrometry**.

## General description from Schmidt et al. (2023)

*To achieve a clean biochemical separation of SARS-CoV-2 genomic RNA (gRNA) and subgenomic mRNA (sgmRNA), we modified the RNA antisense purification strategy implemented in RAP-MS<sup>1,2,3</sup>.
We construct antisense capture probes that hybridize to the ORF1ab region, which is not present in sgmRNAs. This initial capture step enriches gRNAs, while sgmRNAs remain in the flowthrough. After removal of excess capture probes, we add a different pool of antisense probes that hybridize to the S-ORF10 region and capture all SARS-CoV-2 sgmRNAs.*

*Using this strategy, we performed RAP-MS in Huh-7 cells at 24 hours post-infection (hpi)8. We first extracted UV-crosslinked RNA from proteins purified by RAP-MS and sequenced the recovered RNA8. When purifying gRNAs and sgmRNAs separately, crosslinked RNA was overwhelmingly recovered from targeted RNA regions, confirming the successful separation of genomic from subgenomic RNA...*

## Short description

- adapter clipping and quality trimming with cutadapt
- map reads to host and viral genome with star
- separation of mapped reads from host and virus
- normalization of mapped reads

## Recommended requirements

It is recommended to install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and using conda it is an easy task to install [snakemake](https://snakemake.github.io/).

All tools needed for the analysis are installed during the run of Snakemake while all specifications for the different analysis are written in separate yaml environment files.

## Directory structure

Please define the following directory structure so that the snakemake workflow is able to detect all necessary files.

```bash
├── config
│   ├── cluster.json
│   ├── config.yaml
│   └── units.tsv
├── data
├── popart-workflow
├── Log_Err
├── Log_Out
├── references
│    ├── fasta
│    │   ├── GRCh38_SarsCov2.dna.toplevel.fa
│    │   ├── Homo_sapiens.GRCh38.dna_sm.toplevel.fa
│    │   ├── Sars_cov_2.ASM985889v3.dna.toplevel.fa
│    ├── gtf
│    │   ├── GRCh38_SarsCov2.gtf
│    │   ├── Homo_sapiens.GRCh38.106.gtf
│    │   └── Sars_cov_2.ASM985889v3.101.gtf
│    └── STAR_INDEX
│        ├── homo_sapiens_sars_cov2
└── Snakefile
```

# Citation:
1. Schmidt, N., Lareau, C.A., Keshishian, H., Ganskih, S., Schneider, C., Hennig, T., Melanson, R., Werner, S., Wei, Y., Zimmer, M., et al. (2020). The SARS-CoV-2 RNA-protein interactome in infected human cells. Nat Microbiol 6, 339–353. 10.1038/s41564-020-00846-z.

2. McHugh, C.A., Chen, C.-K., Chow, A., Surka, C.F., Tran, C., McDonel, P., Pandya-Jones, A., Blanco, M., Burghard, C., Moradian, A., et al. (2015). The Xist lncRNA interacts directly with SHARP to silence transcription through HDAC3. Nature 521, 232–236. 10.1038/nature14443.

3. Munschauer, M., Nguyen, C.T., Sirokman, K., Hartigan, C.R., Hogstrom, L., Engreitz, J.M., Ulirsch, J.C., Fulco, C.P., Subramanian, V., Chen, J., et al. (2018). The NORAD lncRNA assembles a topoisomerase complex critical for genome stability. Nature 561, 132–136. 10.1038/s41586-018-0453-z.

4. Schmidt, N.,  Ganskih, S., Wei, Y., Gabel, A,  Zielinski, S., Keshishian, H., Lareau, C. A., Zimmermann, L., Makroczyova, J., Pearce, C., Krey, K., Hennig, T., Stegmaier S., Moyon, L., Horlacher, M. , Werner, S., Aydin, J., Olguin-Nava, M., Potabattula, R., Kibe, A., Dölken, L., Smyth, R.P., Caliskan, N., Marsico, A., Krempl, C., Bodem, J., Pichlmair, A., Carr, S.A., Chlanda, P., Erhard, F., and Munschauer, M. (2023): SND1 binds SARS-CoV-2 negative-sense RNA and promotes viral RNA synthesis through NSP9. Cell. 2023 Oct 26;186(22):4834-4850.e23. doi: 10.1016/j.cell.2023.09.002
