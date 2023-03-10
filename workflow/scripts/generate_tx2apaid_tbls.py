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
         utr_bed_path,
         species_str,
         output_prefix):
    '''
    '''

    eprint(f"Reading in BED file used to construct 3'UTR library - {utr_bed_path}")
    bed = pd.read_csv(utr_bed_path,
                      sep="\t",
                      names=["chrom", "start", "end", "name", "score", "Strand", "gene_id"],
                      usecols=["chrom", "start", "end", "name", "Strand"], # score & gene_id not used to construct transcript name or downstream matching
                      dtype=np.string_
                      )

    # eprint(bed)

    # eprint("Reconstructing transcript IDs/names in 3'UTR library for each APA_ID...")


    ## Construct the generated transcript ID from info in BED file
    # ENST00000378424_PRXL2B,ENST00000378425_PRXL2B,ENST00000378427_PRXL2B,ENST00000419916_PRXL2B_hsa_chr1_2589409_2589589_+_utr_2589427_2589589::chr1:2589409-2589589(+)
    # <name>::<chrom>:<start>-<end>(<strand>)

    # 1. Generate 'bedtools getfasta' string
    # e.g. ::chr1:2589409-2589589(+)
    # ::<chrom>:<start>-<end>(<strand>)

    # generate coordinate string
    bed.loc[:, "exon_coords"] = bed["start"].str.cat(bed["end"], sep="-")

    # eprint(bed.exon_coords)

    # last exon start & end coords no longer needed
    bed.drop(columns=["start", "end"], inplace=True)

    # Generate 'strand' string
    bed.loc[:, "strand"] = "(" + bed.Strand + ")"

    # combine start & end with strand
    bed.loc[:, "exon_coords"] = bed["exon_coords"].str.cat(bed["strand"], sep="")

    # now combine start-end(strand) with chrom
    bed.loc[:, "exon_coords"] = bed["chrom"].str.cat(bed["exon_coords"], sep=":")

    # now combine name with exon coords
    # name::chr:start-end(strand)
    bed.loc[:, "getfasta_name"] = bed["name"].str.cat(bed["exon_coords"], sep="::")

    # eprint(bed["getfasta_name"])

    # Extract 3'UTR coordinates from name - allow to join up to correct APA_ID
    # ENSMUST00000000985_ENSMUSG00000000959.7_mmu_chr14_54368317_54368934_+_utr_54368451_54368934
    # 54368451_54368934
    utr_coords = bed["name"].str.split("_utr_", expand=True)[1]
    utr_coords = utr_coords.str.split("_", expand=True)
    utr_coords.rename(columns={0: "UTR3.Start", 1: "UTR3.End"}, inplace=True)

    # eprint(utr_coords)

    bed = bed.merge(utr_coords, left_index=True, right_index=True)

    # eprint(bed)

    bed = bed[["UTR3.Start", "UTR3.End", "Strand", "getfasta_name"]].rename(columns={"getfasta_name": "Transcript_ID"})

    # eprint(bed)

    eprint(f"Reading in PAU results table (output of compute_pau.R) - {pau_tsv_path}")
    pau = pd.read_csv(pau_tsv_path, sep="\t", usecols=pau_annot_cols, dtype=np.string_)

    # eprint(pau)
    # # eprint(pau.dtypes)

    # Join constructed names to pau results to associate APA_IDs with transcript names
    pau = pau.merge(bed, on = ["UTR3.Start", "UTR3.End", "Strand"], how="left")

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
    if len(pau_missing_ids) > 0:
        
        pau_missing_df = pau.loc[pau["Transcript_ID"].isin(pau_missing_ids), ["Transcript_ID", "APA_ID"]]
        # eprint(pau_missing_df)
        
        eprint(f"Writing APA_IDs from PAU results table with missing IDs in quant.sf to file - {output_prefix + '.pau_results_missing_ids.tsv'}")
        pau_missing_df.sort_values(by="Transcript_ID").to_csv(output_prefix + ".pau_results_missing_ids.tsv", sep="\t", index=False, header=True)

    if len(quant_missing_ids) > 0:

        quant_missing_df = quant.copy()
        quant_missing_df = quant_missing_df.loc[quant_missing_df["Name"].isin(quant_missing_ids), :]
        quant_missing_df.loc[:, "APA_ID"] = np.nan
        quant_missing_df.rename(columns={"Name": "Transcript_ID"}, inplace=True)

        # eprint(quant_missing_df)
        
        eprint(f"Writing transcript IDs from quant.sf without associated APA_IDs from PAU results table to file - {output_prefix + '.quant_sf_missing_ids.tsv'}")
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

    parser.add_argument("-b",
                        "--utr",
                        required=True,
                        type=str,
                        default=argparse.SUPPRESS,
                        help="Path to 3'UTR BED file generated by qapa build used to construct transcripts sequences (with qapa fasta)")

    parser.add_argument("--species",
                        default="hsa",
                        type=str,
                        help="Species identifier used to construct reference 3'UTRome library. Provided human annotations are labelled with 'hsa', mouse 'mmu' (qapa build default is 'unk')")

    parser.add_argument("-o",
                        "--output-prefix",
                        default="qapa_annotation",
                        type=str,
                        help="Output prefix for output tables. Suffixed with '.tx2apa.tsv' for table mapping transcript IDs to APA_IDs & '.tx2gene.tsv' for mapping transcript IDs to gene IDs. Any generated transcript IDs absent from PAU results or quant.sf these are output to tables suffixed with '.pau_results_missing_ids.tsv' & '.quant_sf_missing_ids.tsv' respectively")


    if len(sys.argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()


    main(args.pau, args.quant, args.utr, args.species, args.output_prefix)
