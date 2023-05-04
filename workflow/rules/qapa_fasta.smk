def input_utr_bed(use_precomputed_bed: bool,
                  use_custom_polya_bed: bool,
                  extend_upstream: bool,
                  utr_bed: str):
    '''
    Decide target 3'UTR annotation BED file based on workflow parameterisation
    '''

    if use_precomputed_bed:
        return utr_bed
    
    elif extend_upstream:
        return rules.extend_qapa_bed.output

    elif use_custom_polya_bed:
        return rules.qapa_build_custom.output

    else:
        return rules.qapa_build_standard.output
    


rule qapa_fasta_decoy:
    input: 
        utr_bed = input_utr_bed(config["use_precomputed_bed"], config["use_custom_polya_bed"], config["extend_to_upstream_exon"], config["utr_bed"]),
        genome_fa = config["genome_fasta"]

    output:
        utr_fa = os.path.join(config["out_dir"], "qapa_fasta", "utr_sequences.fa"),
        decoys = os.path.join(config["out_dir"], "qapa_fasta", "decoys.txt")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "qapa_fasta", "qapa_fasta_decoy.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "qapa_fasta", "qapa_fasta_decoy.stderr.log"),
    
    container:
        "docker://sambrycesmith/qapa_fork:231dd0c"

    shell:
        """
        qapa fasta \
        -f {input.genome_fa} \
        --decoys \
        -d {output.decoys} \
        {input.utr_bed} {output.utr_fa} \
        2> {log.stderr}
        """ 


rule qapa_fasta:
    input: 
        utr_bed = input_utr_bed(config["use_precomputed_bed"], config["use_custom_polya_bed"], config["extend_to_upstream_exon"], config["utr_bed"]),
        genome_fa = config["genome_fasta"]

    output:
        utr_fa = os.path.join(config["out_dir"], "qapa_fasta", "utr_sequences.fa")

    log:
        stdout = os.path.join(config["out_dir"], "logs", "qapa_fasta", "qapa_fasta_decoy.stdout.log"),
        stderr = os.path.join(config["out_dir"], "logs", "qapa_fasta", "qapa_fasta_decoy.stderr.log"),
    
    container:
        "docker://sambrycesmith/qapa_fork:231dd0c"

    shell:
        """
        qapa fasta \
        -f {input.genome_fa} \
        {input.utr_bed} {output.utr_fa} \
        2> {log.stderr}
        """ 