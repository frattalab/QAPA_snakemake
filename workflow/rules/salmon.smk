# Place any additional SUBWORKFLOWS here.

# was getting weird AttributeError: 'numpy.flagsobj' object has no attribute 'update' with single-end data - looked like it was expecting paired end output (was selecting that rule)
# if single_end:
#     ruleorder: salmon_quant_se > salmon_quant_pe

# else:
#     ruleorder: salmon_quant_pe > salmon_quant_se


rule salmon_index:
    input:
        utr_fa = rules.qapa_fasta_decoy.output.utr_fa if config["add_genome_decoys"] else rules.qapa_fasta.output.utr_fa,
        decoys = rules.qapa_fasta_decoy.output.decoys if config["add_genome_decoys"] else rules.qapa_fasta.output.utr_fa # dummy output to keep snakemake happy

    output:
        seq = os.path.join(config["out_dir"], "salmon", "index", "seq.bin"),
        pos = os.path.join(config["out_dir"], "salmon", "index", "pos.bin")

    params:
        k = config["salmon_kmer_size"],
        outdir = os.path.join(config["out_dir"], "salmon", "index", ""),
        decoys = lambda wildcards, input: " ".join(["--decoys", input.decoys]) if config["add_genome_decoys"] else ""

    threads:
        config["salmon_index_threads"]

    container:
        "docker://quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

    log:
        stdout = os.path.join(config["out_dir"], "logs",
                    "salmon",
                    "salmon_index.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs",
                    "salmon",
                    "salmon_index.stderr.log")


    shell:
        """
        salmon index \
        -t {input.utr_fa} \
        -i {params.outdir} \
        {params.decoys} \
        -k {params.k} \
        -p {threads} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule salmon_quant_pe:
    input:
        fast1 = lambda wildcards: samples.loc[wildcards.sample, "fastq1"],
        fast2 = lambda wildcards: samples.loc[wildcards.sample, "fastq2"],
        index = os.path.join(config["salmon_index_dir"], "seq.bin") if config["use_precomputed_salmon_index"] else rules.salmon_index.output.seq

    output:
        os.path.join(config["out_dir"], "salmon", "quant", "pe", "{sample}", "quant.sf")

    params:
        index_dir = config["salmon_index_dir"] if config["use_precomputed_salmon_index"] else rules.salmon_index.params.outdir,
        output_dir = os.path.join(config["out_dir"], "salmon", "quant", "pe", "{sample}"),
        salmon_extra_parameters = " ".join(config["salmon_extra_parameters"]),
        libtype = "A"

    threads:
        config["salmon_quant_threads"]

    container:
        "docker://quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

    log:
        stdout = os.path.join(config["out_dir"], "logs",
                    "salmon",
                    "salmon_quant_pe.{sample}.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs",
                    "salmon",
                    "salmon_quant_pe.{sample}.stderr.log")



    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        --mates1 {input.fast1} \
        --mates2 {input.fast2} \
        --threads {threads} \
        -o {params.output_dir} \
        {params.salmon_extra_parameters} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule salmon_quant_se:
    input:
        fast1 = lambda wildcards: samples.loc[wildcards.sample, "fastq1"],
        index = os.path.join(config["salmon_index_dir"], "seq.bin") if config["use_precomputed_salmon_index"] else rules.salmon_index.output.seq

    output:
        os.path.join(config["out_dir"], "salmon", "quant", "se", "{sample}", "quant.sf")

    params:
        index_dir = config["salmon_index_dir"] if config["use_precomputed_salmon_index"] else rules.salmon_index.params.outdir,
        output_dir = os.path.join(config["out_dir"], "salmon", "quant", "se", "{sample}"),
        salmon_extra_parameters = " ".join(config["salmon_extra_parameters"]),
        libtype = "A"

    container:
        "docker://quay.io/biocontainers/salmon:1.4.0--hf69c8f4_0"

    log:
        stdout = os.path.join(config["out_dir"], "logs",
                    "salmon",
                    "salmon_quant_se.{sample}.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs",
                    "salmon",
                    "salmon_quant_se.{sample}.stderr.log")

    threads:
        config["salmon_quant_threads"]

    shell:
        """
        salmon quant \
        --index {params.index_dir} \
        --libType {params.libtype} \
        -r {input.fast1} \
        --threads {threads} \
        -o {params.output_dir} \
        {params.salmon_extra_parameters} \
        1> {log.stdout} \
        2> {log.stderr}
        """
