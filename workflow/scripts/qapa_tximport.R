suppressPackageStartupMessages(library(optparse))

# suppressPackageStartupMessages(library(tidyr))
# suppressPackageStartupMessages(library(glue))

option_list <- list(make_option(c("-s", "--sample-table"),
                                type="character",
                                dest="sample_table",
                                help="Path to sample table CSV file used as input to PAPA pipeline. By default the first value in the 'condition' column is taken as the 'base key'"),
                    make_option(c("-d", "--salmon-dir"),
                                dest = "salmon_dir",
                                type = "character",
                                help = "Path to top-level directory under which per-sample Salmon quantification outputs are stored"),
                    make_option(c("-g","--tx2gene"),
                                type="character",
                                help = "Path to <prefix>'.tx2gene.tsv file storing transcript ID to gene ID assignment (output of generate_tx2apaid_tbls.py)"),
                    make_option(c("-c", "--countsFromAbundance"),
                                dest="countsFromAbundance",
                                type="character",
                                default = "dtuScaledTPM",
                                help = "Option to pass to tximport to scale count estimates when generating a combined matrix across samples. Can be one of 'dtuScaledTPM', 'lengthScaledTPM' or 'scaledTPM' ([default= %default])"),
                    make_option(c("-o", "--output-prefix"),
                                dest = "output_prefix",
                                default = "quant_combined",
                                help = "Prefix to names of output files storing per-transcript counts (<output_prefix>.tx_counts.tsv), TPMs (<output_prefix.tx_tpms.tsv) & per-gene summarised counts (<output_prefix>.gene_counts.tsv) & TPMs (<output_prefix>.gene_tpms.tsv) ([default= %default])")
                    )

opt_parser <- OptionParser(option_list = option_list)

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  print_help(opt_parser)
  stop()
}

opt <- parse_args(opt_parser)

suppressPackageStartupMessages(library(tximport))
# suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

### Functions

#' check whether all transcript IDs in tx2gene table are found in the quant.sf file (return True/False)
#' Uses the first column of each df to pull out transcripts
check_tx_overlap <- function(tx2gene, quant_sf) {

    tbl_tx <- tx2gene[, 1] # transcript ID column
    sf_tx <- quant_sf[, 1] # transcript ID column - usually with 'Name' label

    all(sf_tx %in% tbl_tx)

}

### Main script_logic

#1. Read in sample table & construct paths to all Salmon quantification files

# opt_parser$sample_table
sample_tbl <- read.table(opt$sample_table, sep=",", header =TRUE, stringsAsFactors=FALSE)

# identify pe/se datasets based on whether fastq2 is empty (so can construct path to quant.sf)
sample_tbl <- mutate(sample_tbl, seq_type = if_else(is.na(fastq2), "se", "pe"))

quant_paths <- file.path(opt$salmon_dir, sample_tbl$seq_type, sample_tbl$sample_name, "quant.sf")
names(quant_paths) <- sample_tbl$sample_name

if (isFALSE(all(file.exists(quant_paths)))) {
  stop(paste("Not all provided paths to quant.sf files exist, has -d/--salmon-dir been provided correctly?",
             paste(quant_paths, file.exists(quant_paths), sep = " - "), sep = "\n")
       )
}

#2 read in tx2gene - columns should be in that order (names don't matter) but won't be checked by pipeline
tx2gene <- read.table(opt$tx2gene, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# First generate transcript-level matrices - combine across samples, 
# dtuScaledTPM - counts are scaled using the median transcript length among isoforms of a gene, and then the library size
# dtuScaledTPM works such that within a gene, values from all samples and all transcripts get scaled by the same fixed median transcript length.

message("Checking whether all transcripts in quant.sf file are represented in tx2gene")

eg_quant_sf <- read.table(quant_paths[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
all_tx_match <- check_tx_overlap(tx2gene, eg_quant_sf)

if (!all_tx_match) {
    warning("Not all transcripts in quant.sf file are present in tx2gene table. These transcripts will be silently dropped...")
}

if (!all_tx_match & opt$countsFromAbundance == "dtuScaledTPM") {
    stop("dtuScaledTPM requires that all transcripts in quant.sf are present in tx2gene table")
}


message("Generating transcript-level count matrix...")
txi.tx <- tximport(quant_paths,
                   type = "salmon",
                   txOut = TRUE,
                   tx2gene = tx2gene,
                   countsFromAbundance = opt$countsFromAbundance)

# Also generate gene abundance count matrices by summing expression across 3'UTR isoforms
# lengthScaledTPM - scaled using the average transcript length over samples and then the library size ()
# Counts are no longer correlated across samples with transcript length, and so the length offset matrix should not be used for ds analysis

message("Generating gene-level count matrix...")
txi.gene <- tximport(quant_paths,
                     type = "salmon",
                     txOut = FALSE,
                     tx2gene = tx2gene,
                     countsFromAbundance = "lengthScaledTPM")

isoform_id_col <- colnames(tx2gene)[1]
gene_id_col <- colnames(tx2gene)[1]
message(paste("isoform ID column:", isoform_id_col, sep=" "))
message(paste("gene ID column:", gene_id_col, sep=" "))

# Extract count, TPM matrices & convert rownames (tramscript_id/gene_id) to a column
pas_counts <- cbind(rownames(txi.tx$counts), as.data.frame(txi.tx$counts))
pas_tpm <- cbind(rownames(txi.tx$abundance), as.data.frame(txi.tx$abundance))

colnames(pas_counts) <- c(isoform_id_col, colnames(as.data.frame(txi.tx$counts)))
colnames(pas_tpm) <- c(isoform_id_col, colnames(as.data.frame(txi.tx$abundance)))

gene_counts <- cbind(rownames(txi.gene$counts), as.data.frame(txi.gene$counts))
gene_tpm <- cbind(rownames(txi.gene$abundance), as.data.frame(txi.gene$abundance))

colnames(gene_counts) <- c(gene_id_col, colnames(as.data.frame(txi.gene$counts)))
colnames(gene_tpm) <- c(gene_id_col, colnames(as.data.frame(txi.gene$abundance)))


message(paste("Writing transcript counts matrix to - ", opt$output_prefix, ".tx_counts.tsv", sep = ""))

write.table(pas_counts,
          file = paste(opt$output_prefix, ".tx_counts.tsv", sep=""),
          sep="\t",
          col.names = TRUE,row.names=FALSE, quote=FALSE)

message(paste("Writing transcript tpm matrix to - ", opt$output_prefix, ".tx_tpms.tsv", sep = ""))

write.table(pas_tpm,
          file = paste(opt$output_prefix, ".tx_tpms.tsv", sep = ""),
          sep="\t",
          col.names = TRUE,row.names=FALSE, quote=FALSE)

message(paste("Writing gene counts matrix to - ", opt$output_prefix, ".gene_counts.tsv", sep = ""))

write.table(gene_counts,
          file = paste(opt$output_prefix, ".gene_counts.tsv", sep=""),
          sep="\t",
          col.names = TRUE,row.names=FALSE, quote=FALSE)

message(paste("Writing gene tpm matrix to - ", opt$output_prefix, ".gene_tpms.tsv", sep = ""))

write.table(gene_tpm,
          file = paste(opt$output_prefix, ".gene_tpms.tsv", sep = ""),
          sep="\t",
          col.names = TRUE,row.names=FALSE, quote=FALSE)

