from utils import ConfigManager

config = ConfigManager(config)


rule all:
    input:
        config.get_pipeline_files()

rule build_map_wig:
    input:
        config.map_file
    output:
        temp(config.map_template)
    params:
        c=",".join(config.chromosomes),
        w=lambda wc: config.parse_region_size(wc.bin_size)
    conda:
        "envs/hmmcopy-utils.yaml"
    log:
        config.get_log_file(config.map_template)
    resources:
        mem="8G"
    shell:
        "mapCounter -c {params.c} -w {params.w} -s {input} "
        "> {output} "
        "2> {log}"

rule build_gc_wig:
    input:
        config.ref_genome_file
    output:
        temp(config.gc_template)
    params:
        c=",".join(config.chromosomes),
        w=lambda wc: config.parse_region_size(wc.bin_size)
    conda:
        "envs/hmmcopy-utils.yaml"
    log:
        config.get_log_file(config.gc_template)
    resources:
        mem="8G"
    shell:
        "gcCounter -c {params.c} -w {params.w} -s {input} "
        "> {output} "
        "2>{log}"

rule build_reads_chrom:
    input:
        config.get_bam_file
    output:
        temp(config.reads_chrom_template)
    params:
        q=config.min_mqual,
        s=lambda wc: config.parse_region_size(wc.bin_size)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.reads_chrom_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/get_coverage.py "
        "-b {input} "
        "-o {output} "
        "-c {wildcards.chrom} "
        "-q {params.q} "
        "-s {params.s}) >{log} 2>&1"

rule build_reads:
    input:
        lambda wc: [
            str(config.reads_chrom_template).format(chrom=c,bin_size=wc.bin_size,sample=wc.sample)
            for c in config.chromosomes
        ]
    output:
        temp(config.reads_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.reads_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/merge_tables.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"

rule build_rdr:
    input:
        i=config.reads_template,
        n=lambda wc: str(config.reads_template).format(sample="normal",bin_size=wc.bin_size)
    output:
        temp(config.rdr_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/get_rdr.py "
        "-i {input.i} "
        "-n {input.n} "
        "-o {output}) >{log} 2>&1"

rule build_rdr_corrected:
    input:
        i=config.rdr_template,
        g=config.gc_template,
        m=config.map_template
    output:
        config.rdr_corrected_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_corrected_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/get_corrected_coverage.py "
        "-i {input.i} "
        "-g {input.g} "
        "-m {input.m} "
        "-o {output} "
        "--sample-id {wildcards.sample}) >{log} 2>&1"

rule build_allele_counts_chrom:
    input:
        b=config.get_bam_file,
        s=config.get_snp_file
    output:
        temp(config.allele_counts_chrom_template)
    params:
        config.min_bqual
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.allele_counts_chrom_template)
    resources:
        mem="8G"
    shell:
        "(python scripts/get_allele_counts.py "
        "-b {input.b} "
        "-s {input.s} "
        "-o {output} "
        "-c {wildcards.chrom} "
        "-q {params}) >{log} 2>&1"

rule build_allele_counts:
    input:
        lambda wc: [
            str(config.allele_counts_chrom_template).format(chrom=c,prog=wc.prog,sample=wc.sample)
            for c in config.chromosomes
        ]
    output:
        config.allele_counts_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.allele_counts_template)
    shell:
        "(python scripts/merge_tables.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"

rule build_hap_bin_counts:
    input:
        b=config.get_bam_file,
        i=config.allele_counts_template
    output:
        config.allele_bin_counts_template
    params:
        c=" ".join(config.chromosomes),
        s=lambda wc: config.parse_region_size(wc.bin_size)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.allele_bin_counts_template)
    shell:
        "(python scripts/get_allele_bin_counts.py "
        "-b {input.b} "
        "-i {input.i} "
        "-o {output} "
        "-c {params.c} "
        "-s {params.s} "
        "--sample-id {wildcards.sample}) >{log} 2>&1"

rule plot_baf:
    input:
        config.allele_bin_counts_template
    output:
        config.baf_plot_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.baf_plot_template)
    shell:
        "(python scripts/plot_baf.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"

rule plot_rdr:
    input:
        config.rdr_corrected_template
    output:
        config.rdr_plot_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_plot_template)
    shell:
        "(python scripts/plot_rdr.py "
        "-i {input} "
        "-o {output}) >{log} 2>&1"
