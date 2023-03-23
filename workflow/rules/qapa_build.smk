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