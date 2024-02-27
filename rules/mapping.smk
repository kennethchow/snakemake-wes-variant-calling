### Rules for alignemnt and base quality recalibration ###
from os.path import splitext
from pathlib import Path

rule trim_reads_pe:
    input:
        r1 = DATA_DIR / '{sample}_1.fq.gz',
        r2 = DATA_DIR / '{sample}_2.fq.gz',
    output:
        r1_trim=temp(RESULT_DIR / '01_trim' / '{sample}' / '{sample}_1_val_1.fq.gz'),
        r2_trim=temp(RESULT_DIR / '01_trim' / '{sample}' / '{sample}_2_val_2.fq.gz'),
        r1_trimreport = RESULT_DIR / '01_trim' / '{sample}' / '{sample}_1.fq.gz_trimming_report.txt',
    threads: config['resources']['threads']['trim-galore']
    params:
        outprefix=f'{RESULT_DIR}/01_trim/{{sample}}',
        extra = config['params']['trim-galore']
    log: 
        'logs/01_trim/{sample}.log'
    benchmark: 
        'benchmarks/01_trim/{sample}.benchmark'
    shell:
        'mkdir -p {params.outprefix} && trim_galore {params.extra} '
        '--cores {threads} '
        '--paired '
        '--output_dir {params.outprefix} '
        '{input.r1} '
        '{input.r2} &> {log}'


rule map_reads:
    input:
        r1=rules.trim_reads_pe.output['r1_trim'],
        r2=rules.trim_reads_pe.output['r2_trim'],
        idx=rules.bwa_index.output,
    output:
        sortedbam = temp(RESULT_DIR / '02_mapping' / 'tumor' / '{sample}' / '{sample}.sorted.bam'),
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        readgroup = lambda wildcards: extract_read_group(wildcards)
    threads: config['resources']['threads']['bwa_mem']
    log: 
        'logs/02_mapping/bwa_{sample}.log'
    benchmark: 
        'benchmarks/02_mapping/bwa_{sample}.benchmark'
    shell:
        "(bwa mem "
        "-t {threads} "
        "-M "
        "-R '@RG\\tID:{params.readgroup}\\tSM:{wildcards.sample}\\tPL:illumina' "
        "{params.index} "
        "{input.r1} {input.r2} "
        "| samtools sort -@{threads} -o {output.sortedbam} -) &> {log}"



rule mark_duplicates:
    input:
        sortedbam=rules.map_reads.output,
    output:
        bam=temp(RESULT_DIR / "02_mapping" / '{sample}' / '{sample}.dedup.bam'),
        metrics=RESULT_DIR / "00_qc" / "dedup" / '{sample}.metrics.txt',
    log:
        "logs/picard/dedup/{sample}.log",
    benchmark:
        "benchmarks/picard/dedup/{sample}.benchmark",
    params:
        extra=config['params']['picard_MarkDuplicates'],
    shell:
        '(picard MarkDuplicates '
        '{params.extra} '
        'INPUT={input.sortedbam} '
        'OUTPUT={output.bam} '
        'METRICS_FILE={output.metrics}) '
        '2> {log}'


rule recalibrate_base_qualities:
    input:
        bam=rules.mark_duplicates.output['bam'],
        bai=rules.mark_duplicates.output['bam'] + '.bai',
        ref=genome_fasta,
        dictionary=genome_name + ".dict",
        known=output_vcf_path / "variation.noiupac.vcf.gz",
        known_idx=output_vcf_path / "variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=RESULT_DIR / "02_mapping" / '{sample}' / '{sample}.grp',
    log:
        "logs/gatk/bqsr/{sample}.log",
    benchmark:
        "benchmarks/gatk/bqsr/{sample}.benchmark",    
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    shell:
        "(gatk BaseRecalibrator {params.extra} "
        "-R {input.ref} -I {input.bam} "
        "-O {output.recal_table} "
        "--known-sites {input.known} "
        "--sequence-dictionary {input.dictionary}) &> {log}"


rule apply_base_quality_recalibration:
    input:
        bam=rules.mark_duplicates.output['bam'],
        bai=rules.mark_duplicates.output['bam'] + '.bai',
        ref=genome_fasta,
        dictionary=genome_name + ".dict",
        recal_table=rules.recalibrate_base_qualities.output,
    output:
        bam=protected(RESULT_DIR / '02_mapping' / '{sample}' / '{sample}.recal.bam'),
        bai=protected(RESULT_DIR / '02_mapping' / '{sample}' / '{sample}.recal.bam.bai'),
    log:
        "logs/gatk/apply-bqsr/{sample}.log",
    benchmark:
        "benchmarks/gatk/apply-bqsr/{sample}.benchmark",
    params:
        extra=get_regions_param(),
    resources:  
        mem_mb=lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['applybqr'], attempt),
    shell:
        "(gatk ApplyBQSR "
        "-R {input.ref} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.recal_table} "
        "-O {output.bam} "
        "--create-output-bam-index true) "
        "&> {log}"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    params:
        threads=config['resources']['threads']['samtools_index'],
        extra=config['params']['samtools_index']
    shell:
        'samtools index -@ {params.threads} {params.extra} '
        '{input} > {output} 2> {log}'