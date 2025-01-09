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


#' Calculate summary TPM or PAU values (mean, median) from a long-format QAPA results table
#' 
#' @param df A data frame containing QAPA quant results in long format (e.g. output of pivot_qapa_results)
#' @param group_cols Character vector of columns to group by. Default: c("condition", "APA_ID")
#' @param metric_cols Character vector of metric columns to summarize (can be one argument / both). Default: c("TPM", "PAU")
#' @param names_sep String to separate names in pivoted columns. Default: "."
#'
#' @return A data frame in wide format with mean/median values for each metric and condition (one row per APA_ID)
#'      
#' 
#' @examples
#' @export
summarise_qapa_results <- function(df, 
                                     group_cols = c("condition", "APA_ID"),
                                     metric_cols = c("TPM", "PAU"),
                                     names_sep = ".") {
  
  # Validate inputs
  stopifnot(
    is.data.frame(df),
    all(group_cols %in% colnames(df)),
    all(metric_cols %in% colnames(df))
  )
  
  # Get non-APA_ID grouping columns
  pivot_name_cols <- setdiff(group_cols, "APA_ID")
  
  df %>%
    # first calculate summary statistics for each group and metric
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(
        all_of(metric_cols),
        list(
          mean = mean,
          median = median
        ),
        .names = "{.col}{.fn}"
      ),
      .groups = "drop"
    ) %>%
    # put summary values into columns (one row per APA_ID)
    pivot_wider(
      id_cols = APA_ID,
      names_from = all_of(pivot_name_cols),
      values_from = matches("mean|median"),
      names_sep = names_sep
    )
}

test_df <- tibble::tibble(
  APA_ID = rep(c("gene1_1", "gene1_2", "gene2_1", "gene2_2"), each = 2),
  condition = "a",
  replicate = rep(1:2, 4),
  TPM = c(10, 15, 20, 25, 30, 35, 5, 8),
  PAU = c(80, 75, 20, 25, 60, 65, 40, 35)
)

# default (calculate TPM and PAU summaries for each condition and APA_ID)
summarise_qapa_results(test_df)
