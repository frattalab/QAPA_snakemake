library(testthat)
library(here)
library(tidyverse)
source(here("workflow/scripts/calculate_utils.R"))

# Test data setup
test_data_pivot <- tibble(
  APA_ID = c("gene1_1", "gene1_2", "gene2_1", "gene2_2"),
  `sample1.TPM` = c(10, 20, 30, 5),
  `sample2.TPM` = c(15, 25, 35, 8),
  `sample1.PAU` = c(80, 20, 60, 40),
  `sample2.PAU` = c(75, 25, 65, 35)
)


test_that("pivot_qapa_results handles different metrics correctly", {
 # Test TPM conversion
 tpm_result <- pivot_qapa_results(test_data_pivot, "TPM")
 expect_equal(nrow(tpm_result), 8)
 expect_equal(colnames(tpm_result), c("APA_ID", "sample_name", "TPM"))
 expect_true(all(tpm_result$TPM == c(10, 15, 20, 25, 30, 35, 5, 8)))
 
 # Test PAU conversion
 pau_result <- pivot_qapa_results(test_data_pivot, "PAU")
 expect_equal(nrow(pau_result), 8)
 expect_equal(colnames(pau_result), c("APA_ID", "sample_name", "PAU"))
 expect_true(all(pau_result$PAU == c(80, 75, 20, 25, 60, 65, 40, 35)))
 
 # Test both metrics
 both_result <- pivot_qapa_results(test_data_pivot, "both")
 expect_equal(nrow(both_result), 8)
 expect_equal(colnames(both_result), c("APA_ID", "sample_name", "PAU", "TPM"))
 
 # Test PAU sums to 100 for each gene
 pau_sums <- pau_result %>%
   mutate(gene = str_remove(APA_ID, "_[0-9]+$")) %>%
   group_by(gene, sample_name) %>%
   summarise(pau_sum = sum(PAU), .groups = "drop")
 expect_true(all(pau_sums$pau_sum == 100))
})

test_that("custom sample name column works", {
 result <- pivot_qapa_results(test_data_pivot, "both", sample_name_col = "sample_id")
 expect_true("sample_id" %in% colnames(result))
 expect_false("sample_name" %in% colnames(result))
})

test_that("function handles empty data appropriately", {
 empty_data <- test_data_pivot[0,]
 expect_equal(nrow(pivot_qapa_results(empty_data, "TPM")), 0)
 expect_equal(nrow(pivot_qapa_results(empty_data, "PAU")), 0)
 expect_equal(nrow(pivot_qapa_results(empty_data, "both")), 0)
})

test_that("function validates metric argument", {
 expect_error(pivot_qapa_results(test_data_pivot, "invalid"))
})