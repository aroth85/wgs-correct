# from utils import ConfigManager
include: "utils.py"

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
        w=config.parse_region_size(config.bin_size_rdr)
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
        w=config.parse_region_size(config.bin_size_rdr)
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

if config.chromosome_parallel:
    rule build_reads_chrom:
        input:
            b=config.get_bam_file,
            script=workflow.source_path("scripts/get_coverage.py")
        output:
            temp(config.reads_chrom_template)
        params:
            q=config.min_mqual,
            s=config.parse_region_size(config.bin_size_rdr)
        conda:
            "envs/python.yaml"
        log:
            config.get_log_file(config.reads_chrom_template)
        resources:
            mem="8G"
        shell:
            "(python {input.script} "
            "-b {input.b} "
            "-o {output} "
            "-c {wildcards.chrom} "
            "-q {params.q} "
            "-s {params.s}) >{log} 2>&1"

    rule build_reads:
        input:
            i=lambda wc: [
                str(config.reads_chrom_template).format(chrom=c,sample=wc.sample)
                for c in config.chromosomes
            ],
            script=workflow.source_path("scripts/merge_tables.py")
        output:
            temp(config.reads_template)
        conda:
            "envs/python.yaml"
        log:
            config.get_log_file(config.reads_template)
        resources:
            mem="8G"
        shell:
            "(python {input.script} "
            "-i {input.i} "
            "-o {output}) >{log} 2>&1"

else:
    rule build_reads:
        input:
            bam_file=config.get_bam_file,
            script=workflow.source_path("scripts/get_coverage.py")
        output:
            temp(config.reads_template)
        params:
            c=" ".join(config.chromosomes),
            q=config.min_mqual,
            s=config.parse_region_size(config.bin_size_rdr)
        conda:
            "envs/python.yaml"
        log:
            config.get_log_file(config.reads_template)
        resources:
            mem="8G"
        shell:
            "(python {input.script} "
            "-b {input.bam_file} "
            "-o {output} "
            "-c {params.c} "
            "-q {params.q} "
            "-s {params.s}) >{log} 2>&1"


rule build_rdr:
    input:
        i=config.reads_template,
        n=lambda wc: str(config.reads_template).format(sample="normal"),
        script=workflow.source_path("scripts/get_rdr.py")
    output:
        temp(config.rdr_template)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_template)
    resources:
        mem="8G"
    shell:
        "(python {input.script} "
        "-i {input.i} "
        "-n {input.n} "
        "-o {output}) >{log} 2>&1"

rule build_rdr_corrected:
    input:
        i=config.rdr_template,
        g=config.gc_template,
        m=config.map_template,
        script=workflow.source_path("scripts/get_corrected_rdr.py")
    output:
        config.rdr_corrected_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_corrected_template)
    resources:
        mem="8G"
    shell:
        "(python {input.script} "
        "-i {input.i} "
        "-g {input.g} "
        "-m {input.m} "
        "-o {output} "
        "--sample-id {wildcards.sample}) >{log} 2>&1"

if config.chromosome_parallel:
    rule build_allele_counts_chrom:
        input:
            b=config.get_bam_file,
            s=config.snp_file,
            script=workflow.source_path("scripts/get_allele_counts.py")
        output:
            temp(config.allele_counts_chrom_template)
        params:
            q=config.min_bqual
        conda:
            "envs/python.yaml"
        log:
            config.get_log_file(config.allele_counts_chrom_template)
        resources:
            mem="8G"
        shell:
            "(python {input.script} "
            "-b {input.b} "
            "-s {input.s} "
            "-o {output} "
            "-c {wildcards.chrom} "
            "-q {params.q}) >{log} 2>&1"

    rule build_allele_counts:
        input:
            i=lambda wc: [
                str(config.allele_counts_chrom_template).format(chrom=c,sample=wc.sample)
                for c in config.chromosomes
            ],
            script=workflow.source_path("scripts/merge_tables.py")
        output:
            config.allele_counts_template
        conda:
            "envs/python.yaml"
        log:
            config.get_log_file(config.allele_counts_template)
        shell:
            "(python {input.script} "
            "-i {input.i} "
            "-o {output}) >{log} 2>&1"

else:
    rule build_allele_counts:
        input:
            b=config.get_bam_file,
            s=config.snp_file,
            script=workflow.source_path("scripts/get_allele_counts.py")
        output:
            config.allele_counts_template
        params:
            c=" ".join(config.chromosomes),
            q=config.min_bqual
        conda:
            "envs/python.yaml"
        log:
            config.get_log_file(config.allele_counts_template)
        resources:
            mem="8G"
        shell:
            "(python {input.script} "
            "-b {input.b} "
            "-s {input.s} "
            "-o {output} "
            "-c {params.c} "
            "-q {params.q}) >{log} 2>&1"

rule build_hap_bin_counts:
    input:
        b=config.get_bam_file,
        i=config.allele_counts_template,
        script=workflow.source_path("scripts/get_allele_bin_counts.py")
    output:
        config.allele_bin_counts_template
    params:
        c=" ".join(config.chromosomes),
        s=config.parse_region_size(config.bin_size_baf)
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.allele_bin_counts_template)
    shell:
        "(python {input.script} "
        "-b {input.b} "
        "-i {input.i} "
        "-o {output} "
        "-c {params.c} "
        "-s {params.s} "
        "--sample-id {wildcards.sample}) >{log} 2>&1"

rule plot_baf:
    input:
        i=config.allele_bin_counts_template,
        script=workflow.source_path("scripts/plot_baf.py")
    output:
        config.baf_plot_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.baf_plot_template)
    shell:
        "(python {input.script} "
        "-i {input.i} "
        "-o {output}) >{log} 2>&1"

rule plot_rdr:
    input:
        i=config.rdr_corrected_template,
        script=workflow.source_path("scripts/plot_rdr.py")
    output:
        config.rdr_plot_template
    conda:
        "envs/python.yaml"
    log:
        config.get_log_file(config.rdr_plot_template)
    shell:
        "(python {input.script} "
        "-i {input.i} "
        "-o {output}) >{log} 2>&1"
