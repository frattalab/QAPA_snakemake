'''
Steps to convert QAPA/salmon outputs to count matrices and perform differential polyA site usage analysis
'''

rule get_tx2id_tbls:
    input:
        pau = rules.qapa_quant_combined.output,
        sf = expand(os.path.join(config["out_dir"], "salmon", "quant", "{seq_type}", "{sample}", "quant.sf"),
                    zip,
                    seq_type=samples.seq_type.tolist(),
                    sample=samples.index.tolist())[0], # only need one quant.sf file as all use the same reference
        utr_bed = input_utr_bed(config["use_precomputed_bed"], config["use_custom_polya_bed"], config["utr_bed"])

    output:
        tx2apa = os.path.join(config["out_dir"], "annotation", "qapa_annotation.tx2apa.tsv"),
        tx2gene = os.path.join(config["out_dir"], "annotation", "qapa_annotation.tx2gene.tsv")

    params:
        script = os.path.join(config["scripts_dir"], "generate_tx2apaid_tbls.py"),
        species = config["species"],
        output_prefix = os.path.join(config["out_dir"], "annotation", "qapa_annotation")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "differential_apa", "get_tx2id_tbls.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "differential_apa", "get_tx2id_tbls.stderr.log")

    container:
        "docker://quay.io/biocontainers/pyranges:0.0.120--pyh7cba7a3_0" # pandas is dependency so will re-use container

    shell:
        """
        python {params.script} \
        --pau {input.pau} \
        --quant {input.sf} \
        --utr {input.utr_bed} \
        --species {params.species} \
        -o {params.output_prefix} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule tximport:  
    '''
    '''
    input:
        sfs = expand(os.path.join(config["out_dir"], "salmon", "quant", "{seq_type}", "{sample}", "quant.sf"),
               zip,
               seq_type=samples.seq_type.tolist(),
               sample=samples.index.tolist()),
        sample_tbl = config["sample_sheet"],
        tx2gene = rules.get_tx2id_tbls.output.tx2gene if config["tx2gene"] == "" else config["tx2gene"]

    output:
        tx_counts = os.path.join(config["out_dir"], "differential_apa", "all_samples.tx_counts.tsv"),
        gene_counts = os.path.join(config["out_dir"], "differential_apa", "all_samples.gene_counts.tsv")

    params:
        script = os.path.join(config["scripts_dir"], "qapa_tximport.R"),
        salmon_dir = os.path.join(config["out_dir"], "salmon", "quant"),
        output_prefix = os.path.join(config["out_dir"], "differential_apa", "all_samples"),
        cfa = config["counts_from_abundance"]
    
    log:
        stdout = os.path.join(config["out_dir"], "logs", "differential_apa", "tximport.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "differential_apa", "tximport.stderr.log")

    container:
        "docker://sambrycesmith/qapa_snakemake_r:7d789ca83"

    shell:
        """
        Rscript {params.script} \
        -s {input.sample_tbl} \
        -d {params.salmon_dir} \
        -g {input.tx2gene} \
        --countsFromAbundance {params.cfa} \
        -o {params.output_prefix} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule saturn_apa:
    input:
        counts = rules.tximport.output.tx_counts,
        le2gene = rules.get_tx2id_tbls.output.tx2gene if config["tx2gene"] == "" else config["tx2gene"],
        sample_tbl = config["sample_sheet"]

    output:
        saturn_tbl = os.path.join(config["out_dir"], "differential_apa", "saturn_apa.results.tsv"),
        rda = os.path.join(config["out_dir"], "differential_apa", "saturn_apa.image.RData"),
        normed_counts = os.path.join(config["out_dir"], "differential_apa", "saturn_apa.filtered_normed_counts.tsv")

    params:
        script = os.path.join(config["scripts_dir"], "run_differential_usage.R"),
        output_prefix = os.path.join(config["out_dir"], "differential_apa", "saturn_apa"),
        min_mean_count = config["min_mean_count"],
        base_condition = config["base_condition"]
    
    threads:
        config["saturn_threads"]
    
    log:
        stdout = os.path.join(config["out_dir"], "logs", "differential_apa", "saturn_apa.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "differential_apa", "saturn_apa.stderr.log")

    container:
        "docker://sambrycesmith/qapa_snakemake_r:7d789ca83"

    shell:
        """
        Rscript {params.script} \
        -i {input.counts} \
        -g {input.le2gene} \
        -s {input.sample_tbl} \
        -b {params.base_condition} \
        -c {threads} \
        --min-mean-count {params.min_mean_count} \
        -o {params.output_prefix} \
        > {log.stdout} \
        2> {log.stderr}
        """