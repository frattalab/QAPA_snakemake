import pyranges as pr
import pandas as pd
import numpy as np
from functools import reduce
import sys

'''
Script adds penultimate exons to the output BED file of QAPA build, returning a BED12 output file compatible with bedtools getfasta (and downstream QAPA/Salmon steps)
This results in all 3'UTR/APA transcripts receiving the same penultimate exons
tried qapa build with -e 1 (extend by 1 exon), but didn't seem to add to region (as I assumed)
Purpose of doing this is to rescue additional fragments where the left-most read in pair originates from the penultimate exon (e.g. junction spanning reads)
This should increase sensitivity for shorter last exons (more reads assigned) and potentially accuracy of expression estimates too
It will also allow to quantify proximal polyA sites that occur very close to 5'end of terminal exon (e.g. within 100bp - QAPA filters for minimum length of region otherwise Salmon has difficulties mapping reads to the regions)

1. Find introns directly bookending the last exons
'''

def _df_match_3p_adj(df: pd.DataFrame,
                     suffix: str = "_b") -> pd.Series:
    '''Match overlapping regions that are directly adjacent/bookended features at the 3'end of the target feature

    Intended to be applied after a gr.join(gr2, slack=1) call (e.g. <target>.join(<query>, slack=1, strandedness="same")
    So the 3'end of the target is exactly matching the 5'end of the overlapping query (True), otherwise False
    
    The output of this function is compatible with pr.subset()

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    suffix : str, optional
        _description_, by default "_b"

    Returns
    -------
    pd.Series
        like-indexed boolean series
    '''


    if (df.Strand == "+").all():
        # Compare End (3'end of feature) to Start (5' end of overlapping feature)
        decisions = np.where(df['End'] == df['Start' + suffix],
                             True,
                             False)

    elif (df.Strand == "-").all():
        # Compare Start (3'end of feature) to End (5' end of overlapping feature)
        decisions = np.where(df['Start'] == df['End' + suffix],
                             True,
                             False)

    return pd.Series(decisions, index=df.index)


def _df_to_bed12(df: pd.DataFrame,
                 id_col: str = "Name",
                 length_col: str = "exon_length",
                 item_rgb: str = "0,0,255",
                 _score: float = 0.0) -> pd.DataFrame:
    '''Collapse a group of intervals (e.g. exons) to a single row in BED12 format

    intended to be applied pd.Dataframes internal to a pr.PyRanges object (e.g. pr.apply(lambda df: _df_to_bed12(df)))

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    id_col : str, optional
        _description_, by default "Name"
    length_col : str, optional
        _description_, by default "exon_length"
    item_rgb : str, optional
        _description_, by default "0,0,255"
    _score : float, optional
        _description_, by default 0.0

    Returns
    -------
    pd.DataFrame
        pd.DataFrame compatible with pr.PyRanges
    '''
    
    grouped = df.groupby(id_col)
    
    # First define start and end coords for region ()
    chroms = grouped["Chromosome"].first().reset_index()
    strands = grouped["Strand"].first().reset_index()
    starts = grouped["Start"].min().reset_index()
    ends = grouped["End"].max().reset_index()
    
    # blockCount - number of exons
    blockcount = grouped.size().rename("blockCount").reset_index()
    
    # blocksizes collapse exon lengths to comma separated string
    blocksizes = grouped[length_col].agg(lambda col: ",".join(col.astype(str))).rename("blockSizes").reset_index()
    
    # calculate blockStarts: exon Start - starts & collapse to string
    blockstarts = grouped["Start"].apply(lambda x: ",".join((x - x.min()).astype(str))).rename("blockStarts").reset_index()
    
    # define thickStart and thickEnd - default set to start and end
    thickstarts = starts.rename(columns={"Start": "thickStart"})
    thickends = ends.rename(columns={"End": "thickEnd"})
    
    # set a default value for itemRgb & score
    itemrgbs = pd.DataFrame({"Name": df[id_col], "itemRgb": [item_rgb]*len(df)})
    scores = pd.DataFrame({"Name": df[id_col], "Score": [_score]*len(df)})
    
    l = [chroms, starts, ends, scores, strands, thickstarts, thickends, itemrgbs, blockcount, blocksizes, blockstarts]
    out = reduce(lambda s1, s2: pd.merge(s1, s2, on="Name"), l)
    
    
    return out
    

def gr_to_bed12(gr: pr.PyRanges, id_col="Name") -> pr.PyRanges:
    '''_summary_

    Parameters
    ----------
    gr : pr.PyRanges
        _description_
    id_col : str, optional
        _description_, by default "Name"

    Returns
    -------
    pr.PyRanges
        _description_
    '''
    
    # Calculate exon lengths
    gr.exon_length = gr.lengths()
    
    # make sure sorted by start end
    gr = gr.sort()
    
    out_gr = gr.apply(lambda df: _df_to_bed12(df, id_col=id_col))
    
    return out_gr
    
    
    


def main(qapa_bed_path: str,
         gtf_path: str,         
         species_str: str,
         output_bed: str):
    
    '''_summary_

    _extended_summary_
    '''
    
    # read in QAPA BED file of last exons
    qapa_bed = pr.read_bed(qapa_bed_path)
    
    # Extract out transcript IDs from the name field
    # e.g. ENSMUST00000044369_ENSMUSG00000033793.12_mmu_chr1_5162104_5162529_+_utr_5162165_5162529
    # or ENSMUST00000160777_ENSMUSG00000025905.14,ENSMUST00000239100_ENSMUSG00000025905.14_mmu_chr1_5602251_5606152_+_utr_5602871_5606152
    
    # split on species value to extract tx-gene pairs
    tx_gene = qapa_bed.Name.str.split("_" + species_str + "_", expand=True)[0]
    
    # Generate columns of each tx_gene pair for field
    tx_gene = tx_gene.str.split(",", expand=True)
    
    # for each column, split on underscore and retain the transcript_id
    tx_df = tx_gene.apply(lambda col: col.str.split("_", expand=True)[0], axis = "index")
    
    # print(tx_df)
    
    # collapse columns into a single col
    txs = set(tx_df.stack().dropna())
    
    print(f"Number of transcript IDs extracted from qapa build BED file - {len(txs)}")
    
    
    # read in input GTF, retaining transcript and exon entries only
    print("Reading input GTF, can take a while for large references...")
    gtf = pr.read_gtf(gtf_path).subset(lambda df: df["Feature"].isin(["transcript", "exon"]))
    
    # drop to minimal cols to save lugging around un-necessary info
    gtf = gtf[["Feature", "gene_id", "transcript_id", "gene_name"]]

    print("Subsetting GTF for transcripts in qapa build BED file")
    # QAPA strips version ID suffix from transcript IDs - remove to allow matching
    gtf = gtf.assign("transcript_id", lambda df: df["transcript_id"].str.split(".", expand=True)[0])
    
    # print(gtf[["gene_id", "transcript_id"]])
    
    # subset for transcripts in qapa build BED file
    gtf = gtf.subset(lambda df: df.transcript_id.isin(txs))
    
    # print(gtf[["gene_id", "transcript_id"]])
    
    # Double check if transcripts from BED are missing in GTF
    bed_missing_tx = txs - set(gtf.transcript_id)
    print(f"Number of transcript IDs in QAPA BED file missing from GTF - {len(bed_missing_tx)} ")

    # extract introns
    introns = gtf.features.introns(by="transcript")
    
    # print(introns[["Feature", "transcript_id", "gene_id"]])
    # print(introns.columns)
    
    
    print("Finding terminal introns for each last exon in QAPA BED file")
    # Find introns that are directly bookended with the 5'end of last exons 
    # 1. Overlap introns with last exons from BED
    # 2. Subset for introns where 3'end coordinate equals 5'end coord of BED
    introns_m = (introns.join(qapa_bed,
                   how=None, # only keep introns that overlap last exon
                   strandedness="same", 
                   slack=1 # bookended intervals can overlap
                   )
     .subset(lambda df: _df_match_3p_adj(df))
     )
    
    # print(introns_m)
    
    # now extract exons that are directly adjacent to the matching introns (i.e. are the upstream exon)
    print("Finding penultimate exons for each transcript")
    exons = gtf.subset(lambda df: df.Feature == "exon")
    
    exons_m = (exons.join(introns_m,
               how=None, # only keep introns that overlap last exon
               strandedness="same", 
               slack=1, # bookended intervals can overlap
               suffix="_intron"
               )
               .subset(lambda df: _df_match_3p_adj(df, suffix="_intron"))
               )
    
    # print(exons_m)
    # print(exons_m.columns)
    
    # Sometimes upstream last exon will be identical between multiple tx ids - just keep one for merging
    exons_m = exons_m[["Name"]].apply(lambda df: df.drop_duplicates(subset=["Name", "Start", "End"]))
    all_exons = pr.concat([exons_m, qapa_bed])
    
    # print(all_exons)
    
    print("Reformatting transcript regions to BED12 region specifications...")
    bed12 = gr_to_bed12(all_exons)
    
    # Get duplicate rows (maybe because multiple matching penultimate exons?) - remove
    bed12 = bed12.apply(lambda df: df.drop_duplicates(subset=["Name", "Start", "End"]))
    
    
    out_col_order = ["Chromosome", "Start", "End", "Name", "Score", "Strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]

    print(f"Writing BED12 to - {output_bed}")
    bed12.as_df()[out_col_order].to_csv(output_bed, sep="\t",header=False,index=False)
    
    # all_exons.to_bed(output_bed + ".all_exons.bed")
    
    # Output BED12
    

if __name__ == '__main__':
    
    help="""python extend_qapa_bed.py QAPA_BED GTF SPECIES OUTPUT_BED

    Add penultimate exons to transcript models produced by qapa build for PAS quantification with Salmon

    QAPA_BED - BED file produced by qapa build command
    GTF - transcript annotation file in GTF format used as input to qapa build command
    SPECIES - 'species' string used in qapa build command (e.g. 'hsa', 'mmu')
    OUTPUT_BED - name of output BED12 file
    """
    
    if any(i in sys.argv for i in ["-h","--help"]) or len(sys.argv) == 1:
        print(help)
        sys.exit()
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])