#!/usr/bin/env python3

import pyranges as pr
import sys

def main(gencode_polya_gtf,
         out_bed):
    
    
    gtf = pr.read_gtf(gencode_polya_gtf)
    
    # Filter for polyA sites only ('polyA_site' in 3rd field/'Feature' column)
    gtf = gtf.subset(lambda df: df["Feature"] == "polyA_site")
    
    print(f"Number of polyA sites - {len(gtf)}")
    
    # Subset to minimal columns for BED - need score and name field
    # QAPA fills with 'polyA_site' so will do the same
    bed = gtf[["Score", "Feature"]].apply(lambda df: df.rename(columns={"Feature": "Name"}))
    
    bed.to_bed(out_bed)
    

if __name__ == '__main__':
    
    help="""usage: python gencode_polya_gtf2bed.py POLYA_GTF OUTPUT_BED
    
    Extracts 'polyA_site' features from Gencode manual polyA track annotation GTF file and converts to a BED file
    
    Arguments:
    POLYA_GTF - path to Gencde polyA features GTF file
    OUTPUT_BED - name of output BED
    """
    
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print(help)
        sys.exit()
        
    main(sys.argv[1], sys.argv[2])
