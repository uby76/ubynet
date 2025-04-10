#' Classify MolForm into product, resistant, and disappearance
#'
#' This function takes two CSV files or data frames representing "before" and "after" datasets,
#' then classifies MolForm into "product", "resistant", or "disappearance".
#'
#' @param before Data frame or path to CSV file containing the "before" dataset.
#' @param after Data frame or path to CSV file containing the "after" dataset.
#' @param output_file Path to save the classified results (default: NULL, no file saved).
#'
#' @return A data frame with MolForm classification.
#' @export
#'
#' @examples
#' # Using data frames
#' before_data <- data.frame(MolForm = c("C6H12O6", "C2H5OH"))
#' after_data <- data.frame(MolForm = c("C2H5OH", "CH4"))
#' classify_MolForm(before_data, after_data)
#'
#' # Using CSV files
#' classify_MolForm("before.csv", "after.csv", "classified_results.csv")
classify_MolForm <- function(before, after, output_file = NULL) {
  library(dplyr)
  library(tidyr)

  # Read data if input is a file path
  if (is.character(before)) {
    before <- read.csv(before)
  }
  if (is.character(after)) {
    after <- read.csv(after)
  }

  # Check for MolForm column
  if (!"MolForm" %in% colnames(before) | !"MolForm" %in% colnames(after)) {
    stop("Error: The input data must contain a 'MolForm' column.")
  }

  # Extract MolForm and mark presence
  before_formulas <- before %>% select(MolForm) %>% mutate(in_before = TRUE)
  after_formulas  <- after %>% select(MolForm) %>% mutate(in_after = TRUE)

  # Merge and classify
  classified_data <- full_join(before_formulas, after_formulas, by = "MolForm") %>%
    mutate(
      status = case_when(
        !is.na(in_before) & is.na(in_after)  ~ "disappearance",# Only in before: resistant
        is.na(in_before) & !is.na(in_after)  ~ "product",
        !is.na(in_before) & !is.na(in_after) ~ "resistant"
      ),
      value = case_when(
        status == "disappearance"   ~ 0,
        status == "product" ~ 2,
        status == "resistant"  ~ 1,
        TRUE ~ NA_real_
      )
    ) %>%
    select(MolForm, status, value)  # Keep relevant columns

  # Save to CSV if output path is provided
  if (!is.null(output_file)) {
    write.csv(classified_data, output_file, row.names = FALSE)
    message("Classification completed. Results saved to: ", output_file)
  }

  return(classified_data)
}

#' Classify Mass into product, resistant, and resistant
#'
#' This function takes two CSV files or data frames representing "before" and "after" datasets,
#' then classifies Mass into "product", "resistant", or "resistant".
#'
#' @param before Data frame or path to CSV file containing the "before" dataset.
#' @param after Data frame or path to CSV file containing the "after" dataset.
#' @param output_file Path to save the classified results (default: NULL, no file saved).
#'
#' @return A data frame with Mass classification.
#' @export
#'
#' @examples
#' # Using data frames
#' before_data <- data.frame(Mass = c(180.16, 46.07))
#' after_data <- data.frame(Mass = c(46.07, 16.04))
#' classify_Mass(before_data, after_data)
#'
#' # Using CSV files
#' classify_Mass("before.csv", "after.csv", "classified_results.csv")
classify_Mass <- function(before, after, output_file = NULL) {
  library(dplyr)
  library(tidyr)

  # Read data if input is a file path
  if (is.character(before)) {
    before <- read.csv(before)
  }
  if (is.character(after)) {
    after <- read.csv(after)
  }

  # Check for Mass column
  if (!"Mass" %in% colnames(before) | !"Mass" %in% colnames(after)) {
    stop("Error: The input data must contain a 'Mass' column.")
  }

  # Extract Mass and mark presence
  before_masses <- before %>% select(Mass) %>% mutate(in_before = TRUE)
  after_masses  <- after %>% select(Mass) %>% mutate(in_after = TRUE)

  # Merge and classify
  classified_data <- full_join(before_masses, after_masses, by = "Mass") %>%
    mutate(
      status = case_when(
        !is.na(in_before) & is.na(in_after)  ~ "disappearance",
        is.na(in_before) & !is.na(in_after)  ~ "product",
        !is.na(in_before) & !is.na(in_after) ~ "resistant"
      ),
      value = case_when(
        status == "disappearance"   ~ 0,
        status == "product" ~ 2,
        status == "resistant"  ~ 1,
        TRUE ~ NA_real_
      )
    ) %>%
    select(Mass, status, value)  # Keep relevant columns

  # Save to CSV if output path is provided
  if (!is.null(output_file)) {
    write.csv(classified_data, output_file, row.names = FALSE)
    message("Classification completed. Results saved to: ", output_file)
  }

  return(classified_data)
}
    