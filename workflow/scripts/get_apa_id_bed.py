#!/usr/bin/env python3

import pyranges as pr
import pandas as pd
import numpy as np
import sys


def main(utr_bed_path,
         tx2apa_path,
         out_path):
    '''
    '''

    utr_bed = pd.read_csv(utr_bed_path, sep="\t",
    names=["chrom", "start", "end", "name", "score", "Strand", "gene_id"], dtype=np.string_)

    tx2apa = pd.read_csv(tx2apa_path, sep="\t")

    # strip bedtools suffix
    tx2apa.loc[:, "name"] = tx2apa["Transcript_ID"].str.split("::", expand=True)[0]

    # print(tx2apa)

    utr_bed = utr_bed.merge(tx2apa, how="left", on="name")
    utr_bed = utr_bed[["chrom", "start", "end", "APA_ID", "score", "Strand", "gene_id"]]
    utr_bed.to_csv(out_path, sep="\t", index=False, header=False)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
