#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import argparse

'''
Go from a line in compute_pau.R output like:
ENST00000378424,ENST00000378425,ENST00000378427,ENST00000419916	ENSG00000157870	PRXL2B	chr1	2589409	2589589	+	2589427	2589589

to a string like below (transcript_id label in salmon index/quant.sf output file)
ENST00000378424_PRXL2B,ENST00000378425_PRXL2B,ENST00000378427_PRXL2B,ENST00000419916_PRXL2B_hsa_chr1_2589409_2589589_+_utr_2589427_2589589::chr1:2589409-2589589(+)

Notes:
- Order of transcript IDs in compute_pau.R output (i.e. first field) is preserved with order in quant.sf transcript name
- _hsa_ is 'species specific' i.e. will need to provide this label to the script
- utr label seems to be general, but will need to add it manually as not present in compute_pau.R output
'''


# columns from output of compute_pau.R that contain minimal info required to reconstruct tx names in quant.sf
pau_annot_cols = "APA_ID	Transcript	Gene	Gene_Name	Chr	LastExon.Start	LastExon.End	Strand	UTR3.Start	UTR3.End".split("\t")


def eprint(*args, **kwargs):
    '''
    Nice lightweight function to print to STDERR (saves typing, I'm lazy)
    Credit: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python (MarcH)
    '''
    print(*args, file=sys.stderr, **kwargs)


def main(pau_tsv_path,
         quant_sf_path,
         species_str,
         genepred_path,
         output_prefix):
    '''
    '''

    eprint(f"Reading in PAU results table (output of compute_pau.R) - {pau_tsv_path}")
    pau = pd.read_csv(pau_tsv_path, sep="\t", usecols=pau_annot_cols, dtype=np.string_)

    # eprint(pau)
    # eprint(pau.dtypes)
    
    if genepred_path is None:
        eprint("-g/--genePred or path to GTF-derived genePred file not provided - will use 'gene_name' to construct transcript ID")
    
    else:
        eprint(f"Reading in genePred file to extract transcript-gene ID mappings - {genepred_path}")
        # only need gene_id with version number
        genepred = pd.read_csv(genepred_path, sep="\t", usecols=[11], names=["Gene_Name"])
        
        # gene ID versions not present in pau tbl - need a column with numbers removed before join
        genepred.loc[:, "Gene"] = genepred["Gene_Name"].str.split(".", expand=True)[0]
        genepred.drop_duplicates(inplace=True)
        # eprint(genepred)
        
        # join versioned gene_id to pau tbl
        # eprint(pau)
        pau.drop(columns=["Gene_Name"], inplace=True)
        pau = pau.merge(genepred, on="Gene", how="left")
        # eprint(pau)

    eprint("Reconstructing transcript IDs/names in 3'UTR library for each APA_ID...")

    #1. Split 'Transcript' by commma to get list of transcript IDs for each APA_ID
    pau.Transcript = pau.Transcript.str.split(",", expand=False)

    # eprint(pau)

    # for every tx, combine with gene_name
    # e.g. [ENST00000666110, ENST00000217961] --> ENST00000666110_STS,ENST00000217961_STS
    pau.loc[:, "tx_name"] = pau.apply(lambda row: ",".join([tx + "_" + str(row["Gene_Name"]) for tx in row["Transcript"]]), axis="columns")
    pau.drop(columns=["Transcript", "Gene_Name"], inplace=True)

    # eprint(pau)

    #2. create 'sequence info' string
    # e.g. hsa_chr1_2589409_2589589_+_utr_2589427_2589589
    # <species>_<chromosome>_<exon_start>_<exon_end>_<strand>_utr_<utr_start>_<utr_end>

    # First add dummy columns for species and utr
    pau.loc[:, "species"] = species_str
    pau.loc[:, "utr"] = "utr"

    # eprint(pau)

    # Chr	LastExon.Start	LastExon.End	Strand	UTR3.Start	UTR3.End
    pau.loc[:, "sequence"] = pau["species"].str.cat(pau[["Chr", "LastExon.Start", "LastExon.End", "Strand", "utr", "UTR3.Start", "UTR3.End"]], sep="_")

    pau.drop(columns=["species", "utr", "UTR3.Start", "UTR3.End"], inplace=True)

    # eprint(pau)

    # Finally need to generate 'bedtools getfasta' string
    # e.g. ::chr1:2589409-2589589(+)
    # ::<chrom>:<exon_start>-<exon_end>(<strand>)

    # Generate 'strand' string
    pau.Strand = "(" + pau.Strand + ")"

    # eprint(pau.Strand)

    # generate coordinate string
    pau.loc[:, "exon_coords"] = pau["LastExon.Start"].str.cat(pau["LastExon.End"], sep="-")

    # eprint(pau.exon_coords)

    # last exon start & end no longer needed
    pau.drop(columns=["LastExon.Start", "LastExon.End"], inplace=True)

    # generate 'chromosome' string
    pau.Chr = "::" + pau.Chr + ":"

    # eprint(pau.Chr)

    # finally generate 'getfasta string'
    pau.loc[:, "getfasta_coords"] = pau["Chr"].str.cat(pau[["exon_coords", "Strand"]].astype(str), sep="")

    # drop cols no longer needed
    pau.drop(columns=["Chr", "Strand", "exon_coords"], inplace=True)

    # eprint(pau)

    # finally create transcript ID - combo of tx_name, sequence & getfasta_coords
    # ENST00000378424_PRXL2B,ENST00000378425_PRXL2B,ENST00000378427_PRXL2B,ENST00000419916_PRXL2B_hsa_chr1_2589409_2589589_+_utr_2589427_2589589::chr1:2589409-2589589(+)

    # first join sequence & getfasta_coords with no delimiter
    pau = pau.assign(seq_getfasta=lambda df: df["sequence"].str.cat(df["getfasta_coords"], sep="")).drop(columns=["sequence", "getfasta_coords"])

    # eprint(pau)

    pau = pau.assign(Transcript_ID=lambda df: df["tx_name"].str.cat(df["seq_getfasta"], sep="_")).drop(columns=["tx_name", "seq_getfasta"])

    # eprint(pau)

    ### now check that transcript IDs are found in quant.sf files
    # Since use same reference for all quant.sf files, a single file is sufficient
    eprint(f"Reading in quant.sf file - {quant_sf_path}")

    quant = pd.read_csv(quant_sf_path, sep="\t", usecols=["Name"])

    eprint(f"Checking that constructed transcript IDs/Names are present in quant.sf file...")

    quant_tx_ids = set(quant.Name)
    pau_tx_ids = set(pau.Transcript_ID)

    eprint(f"Number of generated tx IDs from PAU results table - {len(pau_tx_ids)}")
    eprint(f"Number of tx IDs in quant.sf table - {len(quant_tx_ids)}")

    pau_missing_ids = pau_tx_ids - quant_tx_ids
    quant_missing_ids = quant_tx_ids - pau_tx_ids

    eprint(f"number of tx IDs from PAU results that are not represented in quant.sf - {len(pau_missing_ids)}")
    eprint(f"number of tx IDs from quant.sf that are not represented in PAU results - {len(quant_missing_ids)}")

    pau_frac_missing = len(pau_missing_ids) / len(pau_tx_ids)

    eprint(f"Fraction of tx IDs generated from pau_results that are missing in quant.sf - {pau_frac_missing}")
    eprint(f"If not 0, consider whether you used inconsistent gene annotation (to make reference 3'UTR library) & identifier versions or provided the wrong species label")

    if pau_frac_missing > 0.5:
        raise Exception(f"At least 50 % of constructed transcript IDs do not match IDs in quant.sf")

    # Output missing IDs to separate TSV files of Transcript_ID | APA_ID
    if len(pau_missing_ids) > 0 or len(quant_missing_ids) > 0:
        eprint("Writing missing transcript IDs to file...")
        pau_missing_df = pau.loc[pau["Transcript_ID"].isin(pau_missing_ids), ["Transcript_ID", "APA_ID"]]

        # eprint(pau_missing_df)

        quant_missing_df = quant.copy()
        quant_missing_df = quant_missing_df.loc[quant_missing_df["Name"].isin(quant_missing_ids), :]
        quant_missing_df.loc[:, "APA_ID"] = np.nan
        quant_missing_df.rename(columns={"Name": "Transcript_ID"}, inplace=True)

        # eprint(quant_missing_df)

        eprint(f"Writing missing IDs from PAU results table - {output_prefix + '.pau_results_missing_ids.tsv'}")
        pau_missing_df.sort_values(by="Transcript_ID").to_csv(output_prefix + ".pau_results_missing_ids.tsv", sep="\t", index=False, header=True)

        eprint(f"Writing missing IDs from quant.sf table - {output_prefix + '.quant_sf_missing_ids.tsv'}")
        quant_missing_df.sort_values(by="Transcript_ID").to_csv(output_prefix + ".quant_sf_missing_ids.tsv", sep="\t", index=False, header=True, na_rep="NA")


    # reorder PAU columns and output to TSV
    pau = pau[["Transcript_ID", "APA_ID", "Gene"]]

    eprint(f"Writing Transcript_ID | APA_ID assignment table to file - {output_prefix + '.tx2apa.tsv'}")
    pau[["Transcript_ID", "APA_ID"]].to_csv(output_prefix + ".tx2apa.tsv", sep="\t", header=True, index=False)
    eprint(f"Writing Transcript_ID | Gene assignment table to file - {output_prefix + '.tx2gene.tsv'}")
    pau[["Transcript_ID", "Gene"]].to_csv(output_prefix + ".tx2gene.tsv", sep="\t", header=True, index=False)


if __name__ == '__main__':

    descrpn = """Generate a table mapping QAPA's APA_IDs to transcript IDs in 3'UTRome fasta/Salmon quantification output"""

    parser = argparse.ArgumentParser(description=descrpn,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, # Add defaults to end of help strings
                                     )

    parser.add_argument("-p",
                        "--pau",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to 'PAU results' table output by qapa quant (compute_pau.R)"
            )

    parser.add_argument("-q",
                        "--quant",
                        required=True,
                        default=argparse.SUPPRESS,
                        type=str,
                        help="Path to quant.sf file output by salmon quant for 3'UTRome/reference library")

    parser.add_argument("--species",
                        default="hsa",
                        type=str,
                        help="Species identifier used to construct reference 3'UTRome library. Provided human annotations are labelled with 'hsa', mouse 'mmu' (qapa build default is 'unk')")
    
    parser.add_argument("-g",
                        "--genePred",
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to extended GenePred format BED file used as input to qapa build. Applicable if BED was generated using gtfToGenePred -genePredExt as done by QAPA_snakemake pipeline, because this will lead to gene IDs being included in the transcript names as opposed to gene names (if e.g. follow QAPA's mysql query). If unsure, double check the 'Name' field in the quant.sf file. If using pre-provided annotations ignore this flag")

    parser.add_argument("-o",
                        "--output-prefix",
                        default="qapa_annotation",
                        type=str,
                        help="Output prefix for output tables. Suffixed with '.tx2apa.tsv' for table mapping transcript IDs to APA_IDs & '.tx2gene.tsv' for mapping transcript IDs to gene IDs. Any generated transcript IDs absent from PAU results or quant.sf these are output to tables suffixed with '.pau_results_missing_ids.tsv' & '.quant_sf_missing_ids.tsv' respectively")


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()


    main(args.pau, args.quant, args.species, args.genePred, args.output_prefix)
