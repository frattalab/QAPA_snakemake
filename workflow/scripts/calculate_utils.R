suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
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

# the TPMs are nonsense, but priority here is just getting PAUs and their deltas correct
test_dfa <- tibble::tibble(
  APA_ID = rep(c("gene1_1", "gene1_2", "gene2_1", "gene2_2"), each = 2),
  condition = "a",
  replicate = rep(1:2, 4),
  TPM = c(10, 15, 20, 25, 30, 35, 5, 8),
  PAU = c(80, 75, 20, 25, 60, 65, 40, 35)
)

# 2nd condition, swap PAUs around between genes
test_dfb <- tibble::tibble(
  APA_ID = rep(c("gene1_1", "gene1_2", "gene2_1", "gene2_2"), each = 2),
  condition = "b",
  replicate = rep(1:2, 4),
  TPM = c(10, 15, 20, 25, 30, 35, 5, 8),
  PAU = c(20, 25, 80, 75, 40, 35, 60, 65)
)

# combine two conditions
test_df <- bind_rows(test_dfa, test_dfb)

# default (calculate TPM and PAU summaries for each condition and APA_ID)
test_df_summ <- summarise_qapa_results(test_df)

#' Generate column name pairs for delta calculations based on contrasts column
#' 
#' @param contrasts_df A dataframe containing contrast_key and base_key columns
#' @param prefix Character vector of prefixes to prepend to column names (to identify columns containing quantification values)
#' @param names_sep Separator between prefix and condition names (default = ".")
#' @param check_df Optional dataframe to verify column existence
#' @return A list of character vectors, each containing pairs of column names
#' @export
get_delta_cols <- function(contrasts_df, prefix, names_sep = ".", check_df = NULL) {
  # contrast df validation
  required_cols <- c("contrast_key", "base_key")
  if (!all(required_cols %in% colnames(contrasts_df))) {
    stop("contrasts_df must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # list of vectors of column names corresponding to comparisons defined in the contrasts table
  # c(contrast_key, base_key)
  # Generate column names for each prefix (i.e. quantification metric/summary)
  col_pairs <- map(prefix, function(pfx) {
    apply(contrasts_df,
          MARGIN = 1,
          function(x) paste(pfx, c(x["contrast_key"], x["base_key"]), sep = names_sep),
          simplify = FALSE
    )
  })
  
  # TODO: would be nice to name the top level of the list (simple with prefix %>% set_names() %>% map(...))

  # Check column existence if check_df provided
  if (!is.null(check_df)) {
    expected_cols <- unique(unlist(col_pairs))
    missing_cols <- setdiff(expected_cols, colnames(check_df))

    if (length(missing_cols) > 0) {
      stop("The following columns are missing from check_df: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  col_pairs
}


#' Calculate differences/deltas between pairs of conditions extracted from contrasts table
#'
#' @param df Input dataframe containing columns to compare
#' @param col_pairs List of list(s) of column name pairs (i.e. output of get_delta_cols). Must be same length as col_prefix
#' @param col_prefix Prefix(es) used in col_pairs (must match what was used in get_delta_cols, used to extract condition from column name)
#' @param output_prefix Prefix for output column names (default = "delta", prepended to col_prefix & conditions)
#' @param output_sep Separator between delta description (<output_prefix><col_prefix>) and the condition names in the output column (default = ".")
#' @param output_key_sep Separator between conditions being compared in the output column (default = "__")
#' @return Input dataframe with additional delta columns
#' @export
calculate_deltas <- function(df, col_pairs, col_prefix, output_prefix = "delta", 
                           output_sep = ".", output_key_sep = "__") {
  
  # Input column name validation
  if (!all(unique(unlist(col_pairs)) %in% colnames(df))) {
    stop("Not all required columns present in input dataframe")
  }
  
  # Calculate deltas for each prefix and its corresponding column pairs
  delta_dfs <- map2(col_prefix, col_pairs, function(pfx, pairs) {
    # Process each contrast pair for this prefix
    pairs %>%
      map(~ {
        # Extract contrast and base condition names from full column names
        contrast_condition <- str_remove_all(.x[1], paste0(pfx, '.'))
        base_condition <- str_remove_all(.x[2], paste0(pfx, '.'))
        
        # Construct output column name from components
        output_colname <- paste0(
          output_prefix,
          pfx,
          output_sep,
          contrast_condition,
          output_key_sep,
          base_condition
        )

        # Calculate delta and return new column
        reframe(df, 
                !!sym(output_colname) := !!sym(.x[1]) - !!sym(.x[2]))
                # "{output_colname}" := !!sym(.x[1]) - !!sym(.x[2]))
      }) %>%
      # combine all comparison columns into a single df
      bind_cols()
  }) 
  
  # now that comparisons computed for each 'prefix' (quantification value summary), combine prefixes into a single df
  delta_dfs <- bind_cols(delta_dfs)
  
  # Add delta columns to original dataframe
  bind_cols(df, delta_dfs)
  
}

test_contrasts_df <- tibble::tibble(comparison_name = "b__vs__a",
                                    contrast_key = "b",
                                    base_key = "a")


# test list of vectors of column names to compare - try just one prefix
get_delta_cols(test_contrasts_df,prefix = "PAUmean", names_sep = ".", check_df = test_df_summ)

# more std usage would be both summary values
test_prefixes <- c("PAUmean", "PAUmedian")
test_col_pairs <- get_delta_cols(test_contrasts_df, prefix = test_prefixes, names_sep = ".", check_df = test_df_summ)
test_col_pairs

# test both pairs
calculate_deltas(test_df_summ, test_col_pairs, test_prefixes)

# try calculating just for PAUmean
calculate_deltas(test_df_summ, 
                 list(test_col_pairs[[1]]), # need to list of lists as input (even if just 1 nested list)
                 test_prefixes[[1]])


