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
        "docker://sambrycesmith/qapa_snakemake_r_dexseq:334b834"

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


rule make_formulas_txt:
    output:
        os.path.join(config["out_dir"], "differential_apa", "formulas.txt")
    
    params:
        full = config["dexseq_formula_full"],
        reduced = config["dexseq_formula_reduced"]
    
    run:
        with open(output[0], "w") as outfile:
            outfile.write(params.full + "\n")
            outfile.write(params.reduced + "\n")


rule dexseq_apa:
    input:
        counts = rules.tximport.output.tx_counts,
        tx2gene = rules.get_tx2id_tbls.output.tx2gene if config["tx2gene"] == "" else config["tx2gene"],
        formulas = rules.make_formulas_txt.output,
        sample_tbl = config["sample_sheet"]

    output:
        os.path.join(config["out_dir"], "differential_apa", "dexseq_apa.{contrast}.results.tsv")

    params:
        script = os.path.join(config["scripts_dir"], "run_dexseq.R"),
        output_prefix = os.path.join(config["out_dir"], "differential_apa", "dexseq_apa.{contrast}"),
        min_mean_count = config["min_mean_count"],
        min_rel_usage = config["min_relative_usage"],
        contrast_name = "{contrast}",
        base_key = lambda wildcards: contrasts.loc[wildcards.contrast, "base_key"],
        contrast_key = lambda wildcards: contrasts.loc[wildcards.contrast, "contrast_key"],
        condition_col = lambda wildcards: contrasts.loc[wildcards.contrast, "column_name"],

    log:
        stdout = os.path.join(config["out_dir"], "logs", "differential_apa", "dexseq_apa.{contrast}.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "differential_apa", "dexseq_apa.{contrast}.stderr.log")

    threads:
        config["dexseq_threads"]
    
    container:
        "docker://sambrycesmith/qapa_snakemake_r_dexseq:334b834"

    shell:
        """
        Rscript {params.script} \
        -i {input.counts} \
        -g {input.tx2gene} \
        -s {input.sample_tbl} \
        --formulas {input.formulas} \
        --base-key {params.base_key} \
        --contrast-key {params.contrast_key} \
        -n {params.contrast_name} \
        --condition-col {params.condition_col} \
        -m {params.min_mean_count} \
        -r {params.min_rel_usage} \
        -c {threads} \
        -o {params.output_prefix} \
        > {log.stdout} \
        2> {log.stderr}
        """


rule process_dexseq_tbl:
    input:
        dexseq_tbl = rules.dexseq_apa.output,
        tx2apa = rules.get_tx2id_tbls.output.tx2apa if config["tx2apa"] == "" else config["tx2apa"],
        qapa_quant = rules.qapa_quant_combined.output

    output:
        os.path.join(config["out_dir"], "differential_apa", "dexseq_apa.{contrast}.results.processed.tsv")


    params:
        script = os.path.join(config["scripts_dir"], "process_dexseq_tbl.R"),
        output_prefix = os.path.join(config["out_dir"], "differential_apa", "dexseq_apa.{contrast}.results")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "differential_apa", "process_dexseq_tbl.{contrast}.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "differential_apa", "process_dexseq_tbl.{contrast}.stderr.log")

    container:
        "docker://sambrycesmith/qapa_snakemake_r_dexseq:334b834"

    shell:
        """
        Rscript {params.script} \
        -i {input.dexseq_tbl} \
        -a {input.tx2apa} \
        -p {input.qapa_quant} \
        -o {params.output_prefix} \
        > {log.stdout} \
        2> {log.stderr}
        """