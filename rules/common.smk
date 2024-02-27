### COMMON FUNCTIONS AND VARIABLES ###
import errno
import pandas as pd
import numpy as np
import os
import psutil
import glob
from pathlib import Path
from os.path import splitext

def tmp_path(path=""):
    """
    if does not exists, create path and return it. If any errors, return
    default path
    :param path: path
    :return: path
    """
    default_path = os.getenv("TMPDIR", config.get("paths").get("java_tmp_dir"))
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
        return path
    return default_path

def get_mem_mb(mem, attempt):
    # Adjust memory based on attempt number
    if attempt == 1:
        return mem * 1024
    elif attempt == 2:
        return mem * 1024 * 1.25
    else:
        return mem * 1024 * 1.5

def java_params(
    tmp_dir=config['paths']['java_tmp_dir'],
    stock_cpu=2,
    custom_ncpu=4,
    custom_mem=4):

    mem_max = custom_mem
    mem_min = 4
    GC_threads = config['java_params']['gc-threads']
    ncpus = custom_ncpu
    tmpdir = tmp_path(tmp_dir)

    params_template = f"'-Xms{mem_min}G -Xmx{mem_max}G -XX:ActiveProcessorCount={ncpus} -XX:ParallelGCThreads={GC_threads} -Djava.io.tmpdir={tmpdir}'"

    return params_template

def extract_read_group(wildcards):
    return '_'.join(wildcards.sample.split('_')[-3:-1])


def get_known_variation_vcfs():
    vcf_files = glob.glob(f"{input_vcf_path}/*.vcf.gz")
    return vcf_files


def get_regions_param(regions=config["references"].get("bed"), default=""):
    if regions:
        params = f"--intervals {regions} "
        padding = config["params"].get("gatk").get("region-padding")
        if padding:
            params += f"--interval-padding {padding} "
        return params
    return default


def extract_sample_name(file_name):
    return '_'.join(file_name.split('_')[:-1])

def check_bwa_index(genome_fasta, bwai_prefix):
    required_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    index_files = [f"{bwai_prefix}{ext}" for ext in required_extensions]
    missing_files = [file for file in index_files if not os.path.exists(file)]
    if missing_files:
        return genome_fasta  # Return the fasta file to trigger indexing
    return []

# sample df preprocessing
sample_df = pd.read_csv(config['paths']['samples'], sep='\t')
sample_df = sample_df.sort_values(sample_df.columns[0])
sample_df['sample_basename'] = sample_df[sample_df.columns[0]].apply(extract_sample_name)

# sample list definitions
samples = sample_df['sample_basename'].unique().tolist()
normal_samples = sample_df[sample_df[sample_df.columns[2]] == 'normal']['sample_basename'].unique()
tumor_samples = sample_df[sample_df[sample_df.columns[2]] == 'tumor']['sample_basename'].unique()
tumor_to_normal = {tumor: normal for normal, tumor in zip(normal_samples, tumor_samples)}

# path definitions
RESULT_DIR = Path(config['paths']['results_dir'])
DATA_DIR = Path(config['paths']['data_dir'])
genome_fasta = config['paths']['genome_fasta']
genome_name = genome_fasta.split('.')[0]
input_vcf_path = Path(config['paths']['input_known_vcf'])
output_vcf_path = Path(f'{input_vcf_path}/output')
bwai_prefix = splitext(genome_fasta)[0]
