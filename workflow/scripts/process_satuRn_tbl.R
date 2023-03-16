suppressPackageStartupMessages(library(optparse))

option_list <- list(make_option(c("-i", "--input-tsv"),
                                type="character",
                                dest = "results",
                                help="Path to <prefix>.results.tsv file storing results table from SatuRn differential usage analysis (output by run_differential_usage.R)"),
                    make_option(c("-a","--tx2apa"),
                                type="character",
                                help = "Path to <prefix>.tx2apa.tsv file storing transcript ID to QAPA APA_ID assignment (output by generate_tx2apaid_tbls.py)"),
                    make_option(c("-p","--pau"),
                                type="character",
                                dest="qapa_quant",
                                help = "Path to PAU results file storing APA_ID and annotation information (output by qapa quant)"),
                    make_option(c("-o", "--output-prefix"),
                                type = "character",
                                dest = "output_prefix",
                                default="differential_usage.results",
                                help = "Prefix to name of augmented SatuRn results output table <prefix>.processed.tsv (default = %default)")
)

opt_parser <- OptionParser(option_list = option_list)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)




suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


### Functions

#' Add a strand-aware proximal-distal order of poly(A) sites along a region (along with total count)
#' Most proximal/5' site is - 1, most distal/3'site is - n (where n = number of PAS in region)
#' region_id_col - can be a single character of vector of column names (recommended to use multiple to account for _PAR_Y genes with same stripped gene ID)
add_polya_site_numbers = function(df,
                                  site_outcol = "Pas_Number",
                                  total_outcol = "Total_Pas",
                                  region_id_col = c("Gene", "Chr"),
                                  event_id_col = "APA_ID",
                                  strand_col = "Strand",
                                  start_col = "LastExon.Start",
                                  end_col = "LastExon.End") {
  
  df %>%
    dplyr::group_by(across(all_of(region_id_col))) %>%
    # group_vars()
    #if + strand then End coord = 3'end - smallest value = left-most = most proximal
    #if - strand then Start coord = 5'end - largest value = left-most = most proximal
    dplyr::mutate("{site_outcol}" := dplyr::if_else(.data[[strand_col]] == '+',
                                              true = row_number(.data[[end_col]]), # smallest value assigned 1
                                              false = row_number(desc(.data[[start_col]])) # largest value assigned 1
                                              ),
    "{total_outcol}" := n_distinct(.data[[event_id_col]])
    ) %>%
    dplyr::ungroup()
    
}

###


# output of run_differential_usage.R
saturn_df <- read_tsv(opt$results)

# Output of qapa quant ...
qapa_df <- read_tsv(opt$qapa_quant, col_select = c(1:12))

# transcript_ID to APA_ID assignment table - output of generate_tx2apaid_tbls.py
tx2apa <- read_tsv(opt$tx2apa,
                   col_names = c("isoform_id", "APA_ID"), # satuRn renames tx_id to isoform_id
                   skip = 1 # skip the header present in file
                   )


# add APA_ID to saturn_df
saturn_df <- left_join(saturn_df, tx2apa, by = "isoform_id")

# Join saturn df with reference information from qapa results df
saturn_df <- left_join(saturn_df, qapa_df, by = "APA_ID") 


# Add a strand-aware polyA site number and number of sites
# May differ to QAPA df if e.g. lowly expressed isoforms have been filtered out
saturn_df <- add_polya_site_numbers(saturn_df)

# Rearrange columns
saturn_df <- saturn_df %>% relocate(APA_ID, Transcript, Gene, Gene_Name, .before = isoform_id) %>% relocate(isoform_id, .after = last_col())

write_tsv(saturn_df, file = paste(opt$output_prefix, ".processed.tsv", sep = ""), col_names = T)