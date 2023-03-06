'''
Run qapa quant per-sample to generate individual qapa results (for Yuk-Kei's challenge outputs script)
Also generate a combined table (run qapa quant on all samples)
'''


# rule qapa_quant_single:
#     input:
#         sf = rules.salmon_quant_se.output if single_end else rules.salmon_quant_pe.output,
#         ids = rules.make_identifiers_tbl.output


#     output:
#         os.path.join(config["out_dir"], "qapa_quant", "{sample}.pau_results.txt")

#     log:
#         stderr = os.path.join(config["out_dir"], "logs", "qapa_fasta", "qapa_quant_single.{sample}.stderr.log"),

#     container:
#         "docker://sambrycesmith/qapa_fork:231dd0c"

#     shell:
#         """
#         qapa quant \
#         --db {input.ids} \
#         {input.sf} > {output} \
#         2> {log.stderr}
#         """


rule qapa_quant_combined:
    input:
        sfs = expand(os.path.join(config["out_dir"], "salmon", "quant", "{seq_type}", "{sample}", "quant.sf"), zip, seq_type=samples.seq_type.tolist(), sample=samples.index.tolist()),
        ids = config["metadata_txt"] if config["use_precomputed_bed"] and config["metadata_txt"] != "" else rules.make_identifiers_tbl.output

    output:
        os.path.join(config["out_dir"], "qapa_quant", "all_samples.pau_results.txt")

    log:
        stderr = os.path.join(config["out_dir"], "logs", "qapa_quant", "qapa_quant_combined.stderr.log"),

    container:
        "docker://sambrycesmith/qapa_fork:408d4ae"

    shell:
        """
        qapa quant \
        --db {input.ids} \
        {input.sfs} > {output} \
        2> {log.stderr}
        """
