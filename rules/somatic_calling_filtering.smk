### VARIANT CALLING ###

rule mutect2:
    input:
        normal_bam = lambda wildcards: f"{RESULT_DIR}/02_mapping/{tumor_to_normal[wildcards.sample_tumor]}/{tumor_to_normal[wildcards.sample_tumor]}.recal.bam",
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
        ref = genome_fasta,
        intervals = config['references']['bed'],
        gnomad = rules.filter_gnomad_vcf_bed.output['filtered_reference_vcf'],
    output:
        vcf = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_mutect2.vcf',
        f1r2 = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_f1r2.tar.gz',
        bam_output = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}.bam',
    params:
        java_opts=java_params(custom_ncpu=config['java_params']['ncpu']['mutect'],
         custom_mem=config['java_params']['mem']['mutect']),
        tumor_lod = config['params']['mutect2']['lod'],
        tumor_id = lambda wildcards: wildcards.sample_tumor,  
        normal_id = lambda wildcards: tumor_to_normal.get(wildcards.sample_tumor, "default_normal"),
    resources:
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['java_params']['mem']['mutect'], attempt),
        cpus = config['java_params']['ncpu']['mutect']
    log:
        'logs/gatk/mutect2/{sample_tumor}.log'
    benchmark:
        'benchmarks/gatk/mutect2/{sample_tumor}.benchmark'
    shell:
        "(gatk --java-options {params.java_opts} Mutect2 "
        "-R {input.ref} "
        "-I {input.tumor_bam} "
        "-I {input.normal_bam} "
        "-normal {params.normal_id} "
        "-tumor {params.tumor_id} "
        "--germline-resource {input.gnomad} "
        "--intervals {input.intervals} " 
        "--f1r2-tar-gz {output.f1r2} "
        "--bam-output {output.bam_output} "
        "--tumor-lod-to-emit {params.tumor_lod} "
        "-O {output.vcf}) &> {log}"

rule manta:
    input:
        normal_bam = lambda wildcards: f"{RESULT_DIR}/02_mapping/{tumor_to_normal[wildcards.sample_tumor]}/{tumor_to_normal[wildcards.sample_tumor]}.recal.bam",
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
        ref = genome_fasta
    output:
        small_indels = RESULT_DIR / '03_somatic_vc' / 'manta' / '{sample_tumor}' / 'results' / 'variants' / 'candidateSmallIndels.vcf.gz',
        manta_run_dir = directory(RESULT_DIR / '03_somatic_vc' / 'manta' / '{sample_tumor}')
    params:
        manta_install_path = "/usr/src/app/manta",
        mem = config['resources']['mem']['manta'],
        ncpu = config['resources']['threads']['manta']  # Adjust as needed
    resources:
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['manta'], attempt),
        cpus = config['resources']['threads']['manta']
    log:
        'logs/manta/{sample_tumor}.log'
    benchmark:
        'benchmarks/manta/{sample_tumor}.benchmark'
    shell:
        "(python2 {params.manta_install_path}/bin/configManta.py "
        "--normalBam {input.normal_bam} "
        "--tumorBam {input.tumor_bam} "
        "--referenceFasta {input.ref} "
        "--runDir {output.manta_run_dir} "
        "--exome "
        "&& python2 {output.manta_run_dir}/runWorkflow.py -g {params.mem} -j {params.ncpu}) &> {log}"

rule strelka2:
    input:
        normal_bam = lambda wildcards: f"{RESULT_DIR}/02_mapping/{tumor_to_normal[wildcards.sample_tumor]}/{tumor_to_normal[wildcards.sample_tumor]}.recal.bam",
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
        ref = genome_fasta,
        bed = config['references']['bed'] + '.gz',
        indel_candidates = rules.manta.output['small_indels'],
    output:
        snv_vcf = RESULT_DIR / '03_somatic_vc' / 'strelka2' / '{sample_tumor}' / 'results' / 'variants' / 'somatic.snvs.vcf.gz',
        indel_vcf = RESULT_DIR / '03_somatic_vc' / 'strelka2' / '{sample_tumor}' / 'results' / 'variants' / 'somatic.indels.vcf.gz',
    params:
        strelka_install_path = "/usr/src/app/strelka2",
        strelka_run_dir = RESULT_DIR / '03_somatic_vc' / 'strelka2' / '{sample_tumor}',
        ncpu = config['resources']['threads']['strelka'] 
    resources:
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['strelka'], attempt),
        cpus = config['resources']['threads']['strelka']
    log:
        'logs/strelka2/{sample_tumor}.log'
    benchmark:
        'benchmarks/strelka2/{sample_tumor}.benchmark'
    shell:
        "({params.strelka_install_path}/bin/configureStrelkaSomaticWorkflow.py "
        "--normalBam {input.normal_bam} "
        "--tumorBam {input.tumor_bam} "
        "--referenceFasta {input.ref} "
        "--runDir {params.strelka_run_dir} "
        "--exome "
        "--callRegions {input.bed} "
        "--indelCandidates {input.indel_candidates} "
        "--reportEVSFeatures "
        "--outputCallableRegions "
        "&& {params.strelka_run_dir}/runWorkflow.py -m local -j {params.ncpu}) &> {log}"

rule lancet:
    input:
        normal_bam = lambda wildcards: f"{RESULT_DIR}/02_mapping/{tumor_to_normal[wildcards.sample_tumor]}/{tumor_to_normal[wildcards.sample_tumor]}.recal.bam",
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
        ref = genome_fasta,
        bed = config['references']['bed']
    output:
        vcf = RESULT_DIR / '03_somatic_vc' / 'lancet' / '{sample_tumor}' / '{sample_tumor}_lancet.vcf'
    params:
        lancet_path = "/usr/src/app/lancet-1.1.0/lancet",
        ncpu = config['resources']['threads']['lancet']
    resources:
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['lancet'], attempt),
        cpus = config['resources']['threads']['lancet']
    log:
        'logs/lancet/{sample_tumor}.log'
    benchmark:
        'benchmarks/lancet/{sample_tumor}.benchmark'
    shell:
        "({params.lancet_path} --tumor {input.tumor_bam} --normal {input.normal_bam} "
        "--ref {input.ref} --bed {input.bed} --num-threads {params.ncpu} > {output.vcf}) &> {log}"

### FILTERING ###

rule get_pileup_summaries_t:
    input:
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
        ref = genome_fasta,
        intervals = config['references']['bed'],
        variants = rules.filter_gnomad_vcf_bed.output['filtered_reference_vcf'],
    output:
        tumor_pileups = RESULT_DIR / '03_somatic_vc' / 'get_pileup_summaries' / 'tumor' / '{sample_tumor}' / '{sample_tumor}_pileups.table',
    params:
        java_opts=java_params(custom_ncpu=config['java_params']['ncpu']['getpileup'],
         custom_mem=config['java_params']['mem']['getpileup']),
    resources:
        cpus = config['java_params']['ncpu']['getpileup'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['java_params']['mem']['getpileup'], attempt)
    log:
        'logs/gatk/get_pileup_summaries/{sample_tumor}.log'
    benchmark:
        'benchmarks/gatk/get_pileup_summaries/{sample_tumor}.benchmark'
    shell:
        "(gatk --java-options {params.java_opts} GetPileupSummaries "
        "-I {input.tumor_bam} "
        "--interval-set-rule INTERSECTION -L {input.intervals} "
        "-V {input.variants} "
        "-L {input.variants} "
        "-O {output.tumor_pileups}) &> {log}"

rule get_pileup_summaries_n:
    input:
        normal_bam = RESULT_DIR / '02_mapping' / '{sample_normal}' / '{sample_normal}.recal.bam',
        ref = genome_fasta,
        intervals = config['references']['bed'],
        variants = rules.filter_gnomad_vcf_bed.output['filtered_reference_vcf'],
    output:
        normal_pileups = RESULT_DIR / '03_somatic_vc' / 'get_pileup_summaries' / 'normal' / '{sample_normal}' / '{sample_normal}_pileups.table',
    params:
        java_opts=java_params(custom_ncpu=config['java_params']['ncpu']['getpileup'],
        custom_mem=config['java_params']['mem']['getpileup']),
    resources:
        cpus = config['java_params']['ncpu']['getpileup'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['java_params']['mem']['getpileup'], attempt)
    log:
        'logs/gatk/get_pileup_summaries/{sample_normal}.log'
    benchmark:
        'benchmarks/gatk/get_pileup_summaries/{sample_normal}.benchmark'
    shell:
        "(gatk --java-options {params.java_opts} GetPileupSummaries "
        "-I {input.normal_bam} "
        "--interval-set-rule INTERSECTION -L {input.intervals} "
        "-V {input.variants} "
        "-L {input.variants} "
        "-O {output.normal_pileups}) &> {log}"

rule learn_read_orientation_model:
    input:
        f1r2 = rules.mutect2.output['f1r2']
    output:
        artifact_prior = RESULT_DIR / '03_somatic_vc' / 'mutect2' / 'artifact_prior' / '{sample_tumor}' / '{sample_tumor}_artifact_prior.tar.gz'
    params:
        java_opts=java_params(custom_ncpu=config['java_params']['ncpu']['learn_read_orientation'],
        custom_mem=config['java_params']['mem']['learn_read_orientation']),
    resources:
        cpus = config['java_params']['ncpu']['learn_read_orientation'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['java_params']['mem']['learn_read_orientation'], attempt)
    log:
        'logs/gatk/read_orientation_model/{sample_tumor}.log'
    benchmark:
        'benchmarks/gatk/learn_read_orientation_model/{sample_tumor}.benchmark'
    shell:
        "(gatk --java-options {params.java_opts} LearnReadOrientationModel "
        "-I {input.f1r2} "
        "-O {output.artifact_prior}) &> {log}"

rule calculate_contamination:
    input:
        tumor_pileups = rules.get_pileup_summaries_t.output['tumor_pileups'],
        normal_pileups = lambda wildcards: f"{RESULT_DIR}/03_somatic_vc/get_pileup_summaries/normal/{tumor_to_normal[wildcards.sample_tumor]}/{tumor_to_normal[wildcards.sample_tumor]}_pileups.table",
    output:
        contamination_table = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_contamination.table',
        tumor_segmentation = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_segments.table'
    params:
        java_opts=java_params(custom_ncpu=config['java_params']['ncpu']['calculate_contamination'],
        custom_mem=config['java_params']['mem']['calculate_contamination']),
    resources:
        cpus = config['java_params']['ncpu']['calculate_contamination'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['java_params']['mem']['calculate_contamination'], attempt)
    log:
        'logs/gatk/calculate_contamination/{sample_tumor}.log'
    benchmark:
        'benchmarks/gatk/calculate_contamination/{sample_tumor}.benchmark'
    shell:
        "(gatk --java-options {params.java_opts} CalculateContamination "
        "-I {input.tumor_pileups} "
        "-O {output.contamination_table} "
        "--tumor-segmentation {output.tumor_segmentation} "
        "--matched {input.normal_pileups}) &> {log}"

rule filter_mutect_calls:
    input:
        unfiltered_vcf = rules.mutect2.output['vcf'],
        ref = genome_fasta,  
        contamination_table = rules.calculate_contamination.output['contamination_table'],
        tumor_segmentation = rules.calculate_contamination.output['tumor_segmentation'],
        ob_priors = rules.learn_read_orientation_model.output['artifact_prior']
    output:
        filtered_vcf = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_mutect2.filtered.vcf',
        filtering_stats = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_mutect2_filtering.stats'
    params:
        java_opts=java_params(custom_ncpu=config['java_params']['ncpu']['filter_mutect_calls'],
        custom_mem=config['java_params']['mem']['filter_mutect_calls'])
    resources:
        cpus = config['java_params']['ncpu']['filter_mutect_calls'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['java_params']['mem']['filter_mutect_calls'], attempt)
    log:
        'logs/gatk/filter_mutect_calls/{sample_tumor}.log'
    benchmark:
        'benchmarks/gatk/filter_mutect_calls/{sample_tumor}.benchmark'
    shell:
        "(gatk --java-options {params.java_opts} FilterMutectCalls "
        "-V {input.unfiltered_vcf} "
        "-R {input.ref} "
        "-O {output.filtered_vcf} "
        "--contamination-table {input.contamination_table} "
        "--tumor-segmentation {input.tumor_segmentation} "
        "--ob-priors {input.ob_priors} "
        "--filtering-stats {output.filtering_stats}) &> {log}"


rule filter_vcf:
    input:
        mutect_vcf = rules.filter_mutect_calls.output['filtered_vcf'],
        strelka_snv_vcf = rules.strelka2.output['snv_vcf'],
        strelka_indel_vcf = rules.strelka2.output['indel_vcf'],
        lancet_vcf = rules.lancet.output['vcf']
    output:
        mutect_filtered_vcf = RESULT_DIR / '03_somatic_vc' / 'mutect2' / '{sample_tumor}' / '{sample_tumor}_mutect2_passfiltered.vcf',
        strelka_filtered_snv_vcf = RESULT_DIR / '03_somatic_vc' / 'strelka2' / '{sample_tumor}' / '{sample_tumor}_filtered_snv.vcf',
        strelka_filtered_indel_vcf = RESULT_DIR / '03_somatic_vc' / 'strelka2' / '{sample_tumor}' / '{sample_tumor}_filtered_indel.vcf',
        lancet_filtered_vcf = RESULT_DIR / '03_somatic_vc' / 'lancet' / '{sample_tumor}' / '{sample_tumor}_lancet_filtered.vcf'
    params:
        tumor_id = lambda wildcards: wildcards.sample_tumor,  
        normal_id = lambda wildcards: tumor_to_normal.get(wildcards.sample_tumor, "default_normal"),
    log:
        'logs/filter_vcf/{sample_tumor}.log'
    benchmark:
        'benchmarks/filter_vcf/{sample_tumor}.benchmark'
    shell:
        """
        (set -euo pipefail;
        grep -E -- '^#|PASS' {input.mutect_vcf} > {output.mutect_filtered_vcf} &&
        grep -E -- '^#|PASS' {input.lancet_vcf} > {output.lancet_filtered_vcf} &&
        zcat {input.strelka_snv_vcf} | grep -E -- '^#|PASS' | sed 's/NORMAL/{params.normal_id}/' | sed 's/TUMOR/{params.tumor_id}/' > {output.strelka_filtered_snv_vcf} &&
        zcat {input.strelka_indel_vcf} | grep -E -- '^#|PASS' > {output.strelka_filtered_indel_vcf}) &> {log}
        """

rule sample_similarity:
    input:
        normal_bam = lambda wildcards: f"{RESULT_DIR}/02_mapping/{tumor_to_normal[wildcards.sample_tumor]}/{tumor_to_normal[wildcards.sample_tumor]}.recal.bam",
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
    output:
        similarity_txt = RESULT_DIR /  '00_qc' / 'sample_similarity' / '{sample_tumor}' / '{sample_tumor}_SampleSimilarity.txt'
    params:
        sample_similarity = "/usr/src/app/ngs-bits/bin/SampleSimilarity"
    log:
        'logs/sample_similarity/{sample_tumor}.log'
    benchmark:
        'benchmarks/sample_similarity/{sample_tumor}.benchmark'
    shell:
        "({params.sample_similarity} -in {input.normal_bam} {input.tumor_bam} "
        "-mode bam -build hg38 "
        "-out {output.similarity_txt}) &> {log}"

rule SOBdetector_mutect:
    input:
        mutect_vcf = rules.filter_vcf.output["mutect_filtered_vcf"],
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
    output:
        mutect_variant_file = RESULT_DIR / '03_somatic_vc' / 'SOBdetector' / '{sample_tumor}' / '{sample_tumor}_SOBdetector_mutect.vcf',
    params:
        sobdetector_jar = "/usr/src/app/SOBdetector/SOBDetector_v1.0.4.jar",
    log:
        'logs/SOBdetector/{sample_tumor}_mutect.log',
    benchmark:
        'benchmarks/SOBdetector/{sample_tumor}_mutect.benchmark',
    shell:
        """
        java -jar {params.sobdetector_jar} \
            --input-type VCF \
            --input-variants {input.mutect_vcf} \
            --input-bam {input.tumor_bam} \
            --output-variants {output.mutect_variant_file} \
            --only-passed true &> {log}
        """

rule SOBdetector_strelka_snv:
    input:
        strelka_snv_vcf = rules.filter_vcf.output["strelka_filtered_snv_vcf"],
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
    output:
        strelka_snv_variant_file = RESULT_DIR / '03_somatic_vc' / 'SOBdetector' / '{sample_tumor}' / '{sample_tumor}_SOBdetector_strelka_snv.vcf',
    params:
        sobdetector_jar = "/usr/src/app/SOBdetector/SOBDetector_v1.0.4.jar",
    log:
        'logs/SOBdetector/{sample_tumor}_strelka_snv.log',
    benchmark:
        'benchmarks/SOBdetector/{sample_tumor}_strelka_snv.benchmark',
    shell:
        """
        java -jar {params.sobdetector_jar} \
            --input-type VCF \
            --input-variants {input.strelka_snv_vcf} \
            --input-bam {input.tumor_bam} \
            --output-variants {output.strelka_snv_variant_file} \
            --only-passed true &> {log}
        """

rule SOBdetector_lancet:
    input:
        lancet_vcf = rules.filter_vcf.output["lancet_filtered_vcf"],
        tumor_bam = RESULT_DIR / '02_mapping' / '{sample_tumor}' / '{sample_tumor}.recal.bam',
    output:
        lancet_variant_file = RESULT_DIR / '03_somatic_vc' / 'SOBdetector' / '{sample_tumor}' / '{sample_tumor}_SOBdetector_lancet.vcf',
    params:
        sobdetector_jar = "/usr/src/app/SOBdetector/SOBDetector_v1.0.4.jar",
    log:
        'logs/SOBdetector/{sample_tumor}_lancet.log',
    benchmark:
        'benchmarks/SOBdetector/{sample_tumor}_lancet.benchmark',
    shell:
        """
        java -jar {params.sobdetector_jar} \
            --input-type VCF \
            --input-variants {input.lancet_vcf} \
            --input-bam {input.tumor_bam} \
            --output-variants {output.lancet_variant_file} \
            --only-passed true &> {log}
        """

rule vcf2maf_mutect:
    input:
        mutect_vcf = rules.SOBdetector_mutect.output['mutect_variant_file'],
    output:
        mutect_maf = RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_mutect.maf',
    params:
        vcf2maf_path = "/data/pipelines/variantcalling/pkgs/vcf2maf-1.6.21.pl",  # path to the vcf2maf script
        ref_fasta = genome_fasta,  # path to the reference genome
        vep_forks = config['resources']['threads']['vcf2maf'],  # number of CPUs to use for VEP
        vep_data = config['references']['known_variants']['vep'],  # path to VEP data
        vep_path = "/usr/src/app/ensembl-vep",  # path to the VEP script
        tumor_id = lambda wildcards: wildcards.sample_tumor,  # use wildcard for tumor ID
        normal_id = lambda wildcards: tumor_to_normal.get(wildcards.sample_tumor, "default_normal"),  # replace with normal sample ID
        ncbi_build = "GRCh38",  # NCBI build
        retain_info = "numF1R2Alt,numF2R1Alt,numF1R2Ref,numF2R1Ref,numF1R2Other,numF2R1Other,SOB,pArtifact,artiStatus",  # Info fields to retain
    resources:
        cpus = config['resources']['threads']['vcf2maf'],
        mem_mb = config['resources']['mem']['vcf2maf'] * 1024
    log:
        'logs/vcf2maf/{sample_tumor}_mutect.log',
    benchmark:
        'benchmarks/vcf2maf/{sample_tumor}_mutect.benchmark',
    shell:
        """
        perl {params.vcf2maf_path} --input-vcf {input.mutect_vcf} --output-maf {output.mutect_maf} \
        --ref-fasta {params.ref_fasta} --vep-overwrite --vep-forks {params.vep_forks} --vep-data {params.vep_data} \
        --vep-path {params.vep_path} --tumor-id {params.tumor_id} --normal-id {params.normal_id} \
        --ncbi-build {params.ncbi_build} --retain-info {params.retain_info} &> {log}
        """

rule vcf2maf_strelka_snv:
    input:
        strelka_snv_vcf = rules.SOBdetector_strelka_snv.output['strelka_snv_variant_file'],
    output:
        strelka_snv_maf = RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_strelka_snv.maf',
    params:
        vcf2maf_path = "/data/pipelines/variantcalling/pkgs/vcf2maf-1.6.21.pl",  # path to the vcf2maf script
        ref_fasta = genome_fasta,  # path to the reference genome
        vep_forks = config['resources']['threads']['vcf2maf'],  # number of CPUs to use for VEP
        vep_data = config['references']['known_variants']['vep'],  # path to VEP data
        vep_path = "/usr/src/app/ensembl-vep",  # path to the VEP script
        tumor_id = lambda wildcards: wildcards.sample_tumor,  # use wildcard for tumor ID
        normal_id = lambda wildcards: tumor_to_normal.get(wildcards.sample_tumor, "default_normal"),  # replace with normal sample ID
        ncbi_build = "GRCh38",  # NCBI build
        retain_info = "numF1R2Alt,numF2R1Alt,numF1R2Ref,numF2R1Ref,numF1R2Other,numF2R1Other,SOB,pArtifact,artiStatus",  # Info fields to retain
    resources:
        cpus = config['resources']['threads']['vcf2maf'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['vcf2maf'], attempt)
    log:
        'logs/vcf2maf/{sample_tumor}_strelka_snv.log',
    benchmark:
        'benchmarks/vcf2maf/{sample_tumor}_strelka_snv.benchmark',
    shell:
        """
        perl {params.vcf2maf_path} --input-vcf {input.strelka_snv_vcf} --output-maf {output.strelka_snv_maf} \
        --ref-fasta {params.ref_fasta} --vep-overwrite --vep-forks {params.vep_forks} --vep-data {params.vep_data} \
        --vep-path {params.vep_path} --tumor-id {params.tumor_id} --normal-id {params.normal_id} \
        --ncbi-build {params.ncbi_build} --retain-info {params.retain_info} &> {log}
        """

rule vcf2maf_strelka_indel:
    input:
        strelka_indel_vcf = rules.filter_vcf.output['strelka_filtered_indel_vcf'],
    output:
        strelka_indel_maf = RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_strelka_indel.maf',
    params:
        vcf2maf_path = "/data/pipelines/variantcalling/pkgs/vcf2maf-1.6.21.pl",  # path to the vcf2maf script
        ref_fasta = genome_fasta,  # path to the reference genome
        vep_forks = config['resources']['threads']['vcf2maf'],  # number of CPUs to use for VEP
        vep_data = config['references']['known_variants']['vep'],  # path to VEP data
        vep_path = "/usr/src/app/ensembl-vep",  # path to the VEP script
        tumor_id = lambda wildcards: wildcards.sample_tumor,  # use wildcard for tumor ID
        normal_id = lambda wildcards: tumor_to_normal.get(wildcards.sample_tumor, "default_normal"),  # replace with normal sample ID
        ncbi_build = "GRCh38",  # NCBI build
        retain_info = "numF1R2Alt,numF2R1Alt,numF1R2Ref,numF2R1Ref,numF1R2Other,numF2R1Other,SOB,pArtifact,artiStatus",  # Info fields to retain
    resources:
        cpus = config['resources']['threads']['vcf2maf'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['vcf2maf'], attempt)
    log:
        'logs/vcf2maf/{sample_tumor}_strelka_indel.log',
    benchmark:
        'benchmarks/vcf2maf/{sample_tumor}_strelka_indel.benchmark',
    shell:
        """
        perl {params.vcf2maf_path} --input-vcf {input.strelka_indel_vcf} --output-maf {output.strelka_indel_maf} \
        --ref-fasta {params.ref_fasta} --vep-overwrite --vep-forks {params.vep_forks} --vep-data {params.vep_data} \
        --vep-path {params.vep_path} --tumor-id {params.tumor_id} --normal-id {params.normal_id} \
        --ncbi-build {params.ncbi_build} --retain-info {params.retain_info} &> {log}
        """

rule vcf2maf_lancet:
    input:
        lancet_vcf = rules.SOBdetector_lancet.output['lancet_variant_file'],
    output:
        lancet_maf = RESULT_DIR / '03_somatic_vc' / 'vcf2maf' / '{sample_tumor}' / '{sample_tumor}_lancet.maf',
    params:
        vcf2maf_path = "/data/pipelines/variantcalling/pkgs/vcf2maf-1.6.21.pl",  # path to the vcf2maf script
        ref_fasta = genome_fasta,  # path to the reference genome
        vep_forks = config['resources']['threads']['vcf2maf'],  # number of CPUs to use for VEP
        vep_data = config['references']['known_variants']['vep'],  # path to VEP data
        vep_path = "/usr/src/app/ensembl-vep",  # path to the VEP script
        tumor_id = lambda wildcards: wildcards.sample_tumor,  # use wildcard for tumor ID
        normal_id = lambda wildcards: tumor_to_normal.get(wildcards.sample_tumor, "default_normal"),  # replace with normal sample ID
        ncbi_build = "GRCh38",  # NCBI build
        retain_info = "numF1R2Alt,numF2R1Alt,numF1R2Ref,numF2R1Ref,numF1R2Other,numF2R1Other,SOB,pArtifact,artiStatus",  # Info fields to retain
    resources:
        cpus = config['resources']['threads']['vcf2maf'],
        mem_mb = lambda wildcards, attempt: get_mem_mb(config['resources']['mem']['vcf2maf'], attempt)
    log:
        'logs/vcf2maf/{sample_tumor}_lancet.log',
    benchmark:
        'benchmarks/vcf2maf/{sample_tumor}_lancet.benchmark',
    shell:
        """
        perl {params.vcf2maf_path} --input-vcf {input.lancet_vcf} --output-maf {output.lancet_maf} \
        --ref-fasta {params.ref_fasta} --vep-overwrite --vep-forks {params.vep_forks} --vep-data {params.vep_data} \
        --vep-path {params.vep_path} --tumor-id {params.tumor_id} --normal-id {params.normal_id} \
        --ncbi-build {params.ncbi_build} --retain-info {params.retain_info} &> {log}
        """

