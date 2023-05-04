rule qapa_build_standard:
    '''
    Perform qapa build step using PolyASite database and Gencode manual PolyA track
    '''
    input:
        genepred = rules.gtf_to_genepred.output,
        tx_db = rules.make_identifiers_tbl.output,
        gencode_pas = config["gencode_polya"] if config["gencode_polya"].endswith(".bed") else rules.make_gencode_polya_bed.output,
        polyasite = config["polyasite_bed"]

    output:
        os.path.join(config["out_dir"], "qapa_build", "utrs.standard_polya.bed")

    params:
        species = config["species"],
        num_extends = config["num_exons_extend_5p"]

    log:
        stdout = os.path.join(config["out_dir"], "logs", "qapa_build", "qapa_build_standard.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "qapa_build", "qapa_build_standard.stderr.log"),
    
    container:
        "docker://sambrycesmith/qapa_fork:231dd0c"

    shell:
        """
        qapa build \
        --db {input.tx_db} \
        --gencode_polya {input.gencode_pas} \
        --polyasite {input.polyasite} \
        --species {params.species} \
        -e {params.num_extends} \
        {input.genepred} \
        > {output} \
        2> {log.stderr}
        """
    
rule qapa_build_custom:
    '''
    Run qapa build using a custom BED file of polyA sites
    '''
    input:
        genepred = rules.gtf_to_genepred.output,
        tx_db = rules.make_identifiers_tbl.output,
        custom_pas = config["custom_polya_bed"]

    output:
        os.path.join(config["out_dir"], "qapa_build", "utrs.custom_polya.bed")
    
    params:
        species = config["species"],
        num_extends = config["num_exons_extend_5p"]


    log:
        stdout = os.path.join(config["out_dir"], "logs", "qapa_build", "qapa_build_custom.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "qapa_build", "qapa_build_custom.stderr.log"),
    
    container:
        "docker://sambrycesmith/qapa_fork:231dd0c"

    shell:
        """
        qapa build \
        --db {input.tx_db} \
        --other {input.custom_pas} \
        --species {params.species} \
        -e {params.num_extends} \
        {input.genepred} \
        > {output} \
        2> {log.stderr}
        """


rule extend_qapa_bed:
    input:
        bed = rules.qapa_build_custom.output if config["use_custom_polya_bed"] else rules.qapa_build_standard.output,
        gtf = config["gtf"]

    output: 
        os.path.join(config["out_dir"], "qapa_build", "utrs.custom_polya.extend_upstream.bed12") if config["use_custom_polya_bed"] else os.path.join(config["out_dir"], "qapa_build", "utrs.standard_polya.extend_upstream.bed12")

    params:
        script = os.path.join(config["scripts_dir"], "extend_qapa_bed.py"),
        species = config["species"]

    log:
        stdout = os.path.join(config["out_dir"], "logs", "qapa_build", "extend_qapa_bed.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "qapa_build", "extend_qapa_bed.stderr.log")

    container:
        "docker://quay.io/biocontainers/pyranges:0.0.120--pyh7cba7a3_0"

    shell:
        '''
        python {params.script} \
        {input.bed} \
        {input.gtf} \
        {params.species} \
        {output} \
        > {log.stdout} \
        2> {log.stderr}
        '''