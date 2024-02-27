import pandas as pd

#### config ####
configfile: "config.yaml"

#### load rules
include: "rules/common.smk"
include: "rules/ref.smk"
include: "rules/mapping.smk"
include: "rules/somatic_calling_filtering.smk"
include: "rules/qc.smk"

def result_files():
    output = list()

    if config['mode'] == 'somatic':
        output.extend([
            RAW_FASTQC,
            TRIMMED_FASTQC,
            TUMOR_BAM,
            NORMAL_BAM,
            MUTECT_MAF,
            STRELKA_SNV_MAF,
            STRELKA_INDEL_MAF,
            LANCET_MAF,
            TUMOR_PILEUPS,
            NORMAL_PILEUPS,
            SAMPLE_SIMILARITY,
            MULTIQC
        ])
    elif config['mode'] == 'map_only':
        output.extend([
            RAW_FASTQC,
            TRIMMED_FASTQC,
            RECAL_BAM,
            MULTIQC
        ])
    else:
        print('Mode not selected correctly, define mode as germline, somatic or map_only in config.yaml')
    
    return output

### PREPROCESSING FOR ANALYSIS-READY BAMS ###
RAW_FASTQC = expand(RESULT_DIR / '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_fastqc.html', sample=samples)
TRIMMED_FASTQC = expand(RESULT_DIR / '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_val_1_fastqc.html', sample=samples)
RECAL_BAM = expand(RESULT_DIR / '02_mapping' / '{sample}' / '{sample}.recal.bam', sample=samples)
TUMOR_BAM = expand(RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam', sample_tumor=tumor_samples)
NORMAL_BAM = expand(RESULT_DIR / '02_mapping' / '{sample_normal}' / '{sample_normal}.recal.bam', sample_normal=normal_samples)

### SOMATIC TUMOR-NORMAL VARIANT CALLING ###
TUMOR_PILEUPS = expand(RESULT_DIR / '03_somatic_vc' / 'get_pileup_summaries' / 'tumor' / '{sample_tumor}' / '{sample_tumor}_pileups.table', sample_tumor=tumor_samples)
NORMAL_PILEUPS = expand(RESULT_DIR / '03_somatic_vc' / 'get_pileup_summaries' / 'normal' / '{sample_normal}' / '{sample_normal}_pileups.table', sample_normal=normal_samples)
SAMPLE_SIMILARITY = expand(RESULT_DIR /  '00_qc' / 'sample_similarity' / '{sample_tumor}' / '{sample_tumor}_SampleSimilarity.txt', sample_tumor=tumor_samples)
MUTECT_MAF = expand(RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_mutect.maf',sample_tumor=tumor_samples)
STRELKA_SNV_MAF = expand(RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_strelka_snv.maf', sample_tumor=tumor_samples)
STRELKA_INDEL_MAF = expand(RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_strelka_indel.maf', sample_tumor=tumor_samples)
LANCET_MAF = expand(RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_lancet.maf', sample_tumor=tumor_samples)
MULTIQC = RESULT_DIR / '00_qc' / 'multiqc' /'multiqc.html'

RESULT_FILES = result_files()

rule all:
    input: RESULT_FILES