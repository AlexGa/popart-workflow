import numpy as np
import pysam
import pandas as pd
import csv
import re
import itertools
from pathlib import Path
import os
import pyBigWig

class InvalidSetException(Exception):
    "Raised when two sets are not equal"
    pass

def read_genome(fasta_file):
    """
    Reading genome fasta file and return a dict with 
    chromosome names and length of the chromosomes
    """
    chr_dict = {}
    chr_length = 0
    chr_name = ""
    with open(fasta_file, 'r') as read_fa:
        for line in read_fa.readlines():
            line = line.strip()
            if line.startswith('>'):
                if chr_length > 0:
                    chr_dict[chr_name] = chr_length
                    chr_length = 0
                chr_name = line.split(' ')[0].replace('>', '')
            else:
                chr_length = chr_length + len(line)
    chr_dict[chr_name] = chr_length
    return(chr_dict)

def paired_read_number(bam_file):
    read_stats = pysam.flagstat(bam_file)
    [line.split("+") for line in read_stats.split("\n")]# 12.Element
    proper_pairs = read_stats.split("\n")
    number_prob_paired = int([line.split("+") for line in read_stats.split("\n")][11][0])
    return(number_prob_paired * 2)

# def count_bin_reads(bamFile, chr, start, end):
#     count_vec = np.zeros(end-start)
#     bam_content = pysam.Samfile(bamFile, "rb")
#     rList = []
#     for p2_rds in bam_content.fetch(chr, max(0,start-2000), end+2000):
#             if p2_rds.qname in rList:
#                 continue
#             if p2_rds.is_reverse:
#                 insPos = p2_rds.pos + p2_rds.alen
#             else:
#                 insPos = p2_rds.pos
#             # if within window
#             if insPos >= start and insPos < end:
#                 for i in range(start - 1, end-1):
#                     count_vec[i] += 1
#                 rList.append(p2_rds.qname)
#     '''
#     pile_iter = bam_content.pileup(contig = chr, start = start, stop = end, stepper = "nofilter")
#     i = 0
#     for pileupcolumn in pile_iter:
#         count_vec[pileupcolumn.pos] = pileupcolumn.n
        
#         print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
#         '''
#     return(count_vec)

def count_bin_reads(bamFile, chr, start, end):
    count_vec = np.zeros(end-start)
    bam_content = pysam.Samfile(bamFile, "rb")
    pile_iter = bam_content.pileup(contig = chr, start = start, stop = end, stepper = "nofilter", max_depth = 9e8)
    i = 0
    for pileupcolumn in pile_iter:
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                count_vec[pileupcolumn.pos] = pileupcolumn.n
        '''
        print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        '''
    return(count_vec)

def calculate_norm_tracks_stranded(bamFile_fwd, bamFile_rev, chr_dict):#, chr, start, end):
    """
    Calculate sequencing depth
    """
    libsize_fwd = paired_read_number(bamFile_fwd)
    libsize_rev = paired_read_number(bamFile_rev)
    libsize = libsize_fwd + libsize_rev
    """
    Count raw reads and normalize reads over lib size
    for each chromosome of the genome
    """
    tracks_dict = {}
    for chr_name in chr_dict.keys():
        """ 
        count base coverage for each chromosome in the genome
        """
        raw_fwd = count_bin_reads(bamFile_fwd, chr_name, 0, chr_dict[chr_name])
        raw_rev = count_bin_reads(bamFile_rev, chr_name, 0, chr_dict[chr_name])
        raw_df = pd.DataFrame({'fwd': raw_fwd.tolist(),
                               'rev': raw_rev.tolist()})
        """
        Normalize bin counts
        """
        norm_fwd = raw_fwd / libsize
        norm_rev = raw_rev / libsize
        norm_df = pd.DataFrame({'fwd': norm_fwd.tolist(),
                                'rev': norm_rev.tolist()})
        """
        Normalize bin counts
        """
        rpm_fwd = norm_df['fwd'] * 1e6
        rpm_rev = norm_df['rev'] * 1e6
        rpm_df = pd.DataFrame({'fwd': rpm_fwd.tolist(),
                               'rev': rpm_rev.tolist()})
        """
        sum into data frame
        """
        tracks_dict[chr_name] = {'raw': raw_df, 'norm': norm_df, 'rpm': rpm_df}
    return(tracks_dict)

def calculate_norm_tracks_paired(bamFile, chr_dict):#, chr, start, end):
    """
    Calculate sequencing depth
    """
    libsize = paired_read_number(bamFile)
    """
    Count raw reads and normalize reads over lib size
    for each chromosome of the genome
    """
    tracks_dict = {}
    for chr_name in chr_dict.keys():
        """ 
        count base coverage for each chromosome in the genome
        """
        raw_reads = count_bin_reads(bamFile, chr_name, 0, chr_dict[chr_name])
        raw_df = pd.DataFrame({'paired': raw_reads.tolist()})
        """
        Normalize bin counts
        """
        norm_reads = raw_reads / libsize
        norm_df = pd.DataFrame({'paired': norm_reads.tolist()})
        """
        Normalize bin counts
        """
        rpm_reads = norm_reads * 1e6
        rpm_df = pd.DataFrame({'paired': rpm_reads.tolist()})
        """
        sum into data frame
        """
        tracks_dict[chr_name] = {'raw': raw_df, 'norm': norm_df, 'rpm': rpm_df}
    return(tracks_dict)

def log2fc_tracks(ip_dict, sm_dict):
    lfc2_dict = {}
    for chr_name in ip_dict.keys():
        lfc2_dict[chr_name] = np.log2(ip_dict[chr_name]['norm'] + 1e-6) - np.log2(sm_dict[chr_name]['norm'] + 1e-6)
    return(lfc2_dict)

def subtracted_tracks(ip_dict, sm_dict):
    subtracted_dict = {}
    for chr_name in ip_dict.keys():
        subtracted_dict[chr_name] = ip_dict[chr_name]['norm'] - sm_dict[chr_name]['norm']
    return(subtracted_dict)

def kl_norm_tracks(ip_dict, lfc_dict):
    kl_dict = {}
    for chr_name in ip_dict.keys():
        kl_dict[chr_name] = ip_dict[chr_name]['norm'] * lfc_dict[chr_name]
    return(kl_dict)

def write_bed(filename, gen_dict, values):
    """
    Write bed file 
    """
    chroms = list(itertools.chain(*[[key]*value for key, value in gen_dict.items()]))
    starts = list(itertools.chain(*[range(1, value+1) for value in gen_dict.values()]))
    bed_df = pd.DataFrame({'chr': chroms,
                           'start': starts,
                           'end': starts,
                           'strand': '*',
                           'score': values})
    bed_df.to_csv(filename, sep = "\t", header = False, index = False)

def write_bw(filename, gen_dict, values):
    """
    Write bigwig file
    """
    bw = pyBigWig.open(filename, "w")
    bw.addHeader(list(gen_dict.items()))
    chroms = list(itertools.chain(*[[key]*value for key, value in gen_dict.items()]))
    starts = list(itertools.chain(*[range(1, value+1) for value in gen_dict.values()]))
    bw.addEntries(chroms, starts = [value - 1 for value in starts], ends = starts, values = values)
    bw.close()


bam_files = {'fwd': "results/extract/only_sars_cov2/SCoV2_gRNA-rep2/SCoV2_gRNA-rep2_fwd.bam", 
             'rev': "results/extract/only_sars_cov2/SCoV2_gRNA-rep2/SCoV2_gRNA-rep2_rev.bam"}

output_prefix_raw = "results/normed_tracks/only_sars_cov2/Raw/SCoV2_gRNA-rep2/SCoV2_gRNA-rep2_RAW"
output_prefix_rpm = "results/normed_tracks/only_sars_cov2/RPM/SCoV2_gRNA-rep2/SCoV2_gRNA-rep2_RPM"

# bam_files = {'paired': "results/extract/only_sars_cov2/SCoV2_sgmRNA-rep2/SCoV2_sgmRNA-rep2_paired.bam"}
genome_fa = "references/fasta/Sars_cov_2.ASM985889v3.dna.toplevel.fa"

genome_fa = snakemake.input.genome_fasta

n = len(snakemake.input) - 1

bam_files = {}

if n == 1:
    bam_files = {'paired': snakemake.input[0]}
if n == 2:
    bam_files = {'fwd': snakemake.input.fwd, 
                 'rev': snakemake.input.rev}

gen_dict = read_genome(genome_fa)

output_prefix_raw = snakemake.params.raw
output_prefix_rpm = snakemake.params.rpm

bin_dict = {}
if n == 1:
    bin_dict = calculate_norm_tracks_paired(bam_files['paired'], gen_dict)
if n == 2:
    bin_dict = calculate_norm_tracks_stranded(bam_files['fwd'], bam_files['rev'], gen_dict)

for strand in bam_files.keys():
    bin_raw_list = np.array(list(itertools.chain(*[sub_dict['raw'][strand].tolist() for sub_dict in bin_dict.values()])))
    raw_file = output_prefix_raw + "_" + strand + ".bw"
    write_bw(raw_file, gen_dict, bin_raw_list)
    rpm_raw_list = np.array(list(itertools.chain(*[sub_dict['rpm'][strand].tolist() for sub_dict in bin_dict.values()])))
    rpm_file = output_prefix_rpm + "_" + strand + ".bw"
    write_bw(rpm_file, gen_dict, rpm_raw_list)


