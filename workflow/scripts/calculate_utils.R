suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

#' Convert QAPA results from wide to long format
#'
#' @param qapa_df A data frame containing QAPA results in wide format with .TPM and .PAU columns
#' @param metric Character string specifying which metric to pivot: "TPM", "PAU", or "both"
#' @param sample_name_col Name for the output column containing sample names (default = "sample_name")
#'
#' @return For metric="TPM" or "PAU", returns a data frame in long format (1 row per APA_ID and sample) with columns 
#'         sample_name_col and the specified metric. For metric="both", returns a data frame 
#'         with both TPM and PAU columns for each APA_ID and sample_name_col combination
#' @export
#'
#' @examples
#' pivot_qapa_results(qapa_df, metric = "both", sample_name_col = "sample_id")
pivot_qapa_results <- function(qapa_df, metric = c("TPM", "PAU", "both"), 
                              sample_name_col = "sample_name") {
 metric <- match.arg(metric)
 
 if (metric %in% c("PAU", "both")) {
   pau_long <- qapa_df %>%
     select(APA_ID, ends_with(".PAU")) %>%
     pivot_longer(
       cols = ends_with(".PAU"),
       names_to = sample_name_col,
       values_to = "PAU",
       names_pattern = "(.*).PAU"
     )
 }
 
 if (metric %in% c("TPM", "both")) {
   tpm_long <- qapa_df %>%
     select(APA_ID, ends_with(".TPM")) %>%
     pivot_longer(
       cols = ends_with(".TPM"),
       names_to = sample_name_col,
       values_to = "TPM",
       names_pattern = "(.*).TPM"
     )
 }
 
  if (metric == "TPM") return(tpm_long)
  if (metric == "PAU") return(pau_long)
 
  inner_join(pau_long, tpm_long, by = c("APA_ID", sample_name_col))
}