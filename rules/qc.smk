### QC ###

rule multiqc:
    input:
        expand(RESULT_DIR / '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_fastqc.html', sample=samples),
        expand(RESULT_DIR / '01_trim' / '{sample}' / '{sample}_1.fq.gz_trimming_report.txt', sample=samples),
        expand(RESULT_DIR / '00_qc' / 'dedup' / '{sample}.metrics.txt', sample=samples),
        expand(RESULT_DIR / '02_mapping' / '{sample}' / '{sample}.recal.bam', sample=samples)  #require recalibrated bam for baserecalibration reports
    output:
        f'{RESULT_DIR}/00_qc/multiqc/multiqc.html'
    params:
        target_dir= RESULT_DIR,
        params=config['params']['multiqc'],
        outdir=f'{RESULT_DIR}/00_qc/multiqc',
        outname="multiqc.html",
    log:
        "logs/multiqc/multiqc.log"
    shell:
        "multiqc "
        "{params.target_dir} "
        "{params.params} "
        "-o {params.outdir} "
        "-n {params.outname} "
        ">& {log}"

rule fastqc_raw:
    input:
        r1 = DATA_DIR / '{sample}_1.fq.gz',
        r2 = DATA_DIR / '{sample}_2.fq.gz',
    output:
        html = RESULT_DIR /  '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_fastqc.html',
        zip = RESULT_DIR /  '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_fastqc.zip',
    threads: 
        config['resources']['threads']['fastqc_raw']
    params:
        outprefix = f'{RESULT_DIR}/00_qc/fastqc/{{sample}}'
    log:
        "logs/00_qc/{sample}.log"
    benchmark:
        'benchmarks/00_qc/{sample}.benchmark'
    shell:
        'mkdir -p {params.outprefix} && fastqc {input.r1} {input.r2} '
        '--threads {threads} '
        '--outdir {params.outprefix} '
        '--extract '
        '--format fastq '  
        '--dir {params.outprefix} &> {log}'

rule fastqc_trimmed:
    input:
        r1=rules.trim_reads_pe.output['r1_trim'],
        r2=rules.trim_reads_pe.output['r2_trim'],
    output:
        html = RESULT_DIR / '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_val_1_fastqc.html',
        zip = RESULT_DIR / '00_qc' / 'fastqc' / '{sample}' / '{sample}_1_val_1_fastqc.zip',
    threads: 
        config['resources']['threads']['fastqc_trimmed']
    params:
        outprefix = f'{RESULT_DIR}/00_qc/fastqc/{{sample}}'
    log:
        "logs/00_qc/{sample}.log"
    benchmark:
        "benchmarks/00_qc/{sample}.benchmark"
    shell:
        'mkdir -p {params.outprefix} && fastqc {input.r1} {input.r2} '
        '--threads {threads} '
        '--outdir {params.outprefix} '
        '--extract '
        '--format fastq '  
        '--dir {params.outprefix} &> {log}'

