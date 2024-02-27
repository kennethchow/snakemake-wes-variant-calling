import os

rule genome_dict:
    input:
        f"{genome_fasta}",
    output:
        f"{genome_name}.dict",
    log:
        "logs/samtools/create_dict.log",
    shell:
        "samtools dict {input} > {output} 2> {log} "

rule bwa_index:
    input:
        lambda wildcards: check_bwa_index(genome_fasta, bwai_prefix),
    output:
        multiext(f'{bwai_prefix}', '.amb', '.ann', '.bwt', '.pac', '.sa'),
    log:
        'logs/02_mapping/bwa_index.log',
    resources:
        mem_mb=lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['bwaindex'], attempt),
    params:
        block_size = config['params']['bwa-mem']['block_size'],
        prefix = bwai_prefix
    shell:
        'bwa index '
        '-b {params.block_size}M '
        '-p {params.prefix} '
        '-a bwtsw '
        '{input} '
        '&> {log}'

rule generate_fai:
    input:
        f"{genome_fasta}",
    output:
        f"{genome_fasta}.fai",
    log:
        'logs/samtools_faidx.log'
    shell:
        'samtools faidx {input}'

rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        vcfs = get_known_variation_vcfs(),
        fai = rules.generate_fai.output,
    output:
        vcf=f"{output_vcf_path}/variation.vcf.gz",
    log:
        "logs/get-known-variants.log",
    shell:
       'bcftools concat -Oz --naive {input.vcfs} | '
       'bcftools reheader --fai {input.fai} -o {output.vcf} '
       '&> {log}'


rule remove_iupac_codes:
    input:
        rules.get_known_variation.output['vcf'],
    output:
        f"{output_vcf_path}/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | "
        "bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        rules.remove_iupac_codes.output,
    output:
        f"{output_vcf_path}/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/tabix/variation.log",
    shell:
        'tabix -p vcf {input} &> {log}'

rule sort_bed:
    input:
        raw_bed = config['references']['bed']
    output:
        sorted_bed = f"{os.path.splitext(config['references']['bed'])[0]}.sorted.bed"
    shell:
        'sort -k1,1V -k2,2n {input.raw_bed} > {output.sorted_bed}'

rule filter_gnomad_vcf_bed:
    input:
        reference_vcf = config['references']['known_variants']['gnomad'],
        sorted_bed = rules.sort_bed.output['sorted_bed']
    output:
        filtered_reference_vcf = f"{os.path.splitext(config['references']['known_variants']['gnomad'])[0]}.filtered.vcf.gz",
    shell:
        'bcftools view -R {input.sorted_bed} {input.reference_vcf} -Oz -o {output.filtered_reference_vcf}'