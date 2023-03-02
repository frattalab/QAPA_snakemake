rule make_identifiers_tbl:
    input:
        config["gtf"]

    output:
        os.path.join(config["out_dir"], "annotation", "annotation_identifers.txt")

    params:
        script = os.path.join(config["scripts_dir"], "GTFtoMartExport.py")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "preprocess_annotation", "make_identifiers_tbl.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "preprocess_annotation", "make_identifiers_tbl.stderr.log"),
        
    shell:
        """
        python {params.script} \
        {input} {output} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule gtf_to_genepred:
    input:
        config["gtf"] 
    output:
        os.path.join(config["out_dir"], "annotation", "genes.genePred")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "preprocess_annotation", "make_gencode_polya_bed.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "preprocess_annotation", "make_gencode_polya_bed.stderr.log")

    container:
        "docker://quay.io/biocontainers/ucsc-gtftogenepred:377--ha8a8165_5"

    shell:
        """
        gtfToGenePred -genePredExt {input} {output}
        """ 


rule make_gencode_polya_bed:
    input:
        config["gencode_polya"]
    
    output:
        os.path.join(config["out_dir"], "annotation", "gencode.polyA_sites.bed")

    params:
        script = os.path.join(config["scripts_dir"], "gencode_polya_gtf2bed.py")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "preprocess_annotation", "make_gencode_polya_bed.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "preprocess_annotation", "make_gencode_polya_bed.stderr.log")

    container:
        "docker://quay.io/biocontainers/pyranges:0.0.120--pyh7cba7a3_0"
        
    shell:
        """
        python {params.script} \
        {input} {output} \
        1> {log.stdout} \
        2> {log.stderr}
        """