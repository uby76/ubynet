#' Merge Mass and Intensity from Multiple CSV Files
#'
#' @param dir_path The directory containing CSV files.
#' @param output_mass_intensity Path to save the merged Mass-Intensity file.
#' @param output_mass_elements Path to save the Mass-Elements file.
#' @param output_mass_only Path to save the Mass-only file.
#'
#' @return NULL (files are saved to disk).
#' @export
merge_mass_intensity <- function(dir_path, output_mass_intensity, output_mass_elements) {
  library(dplyr)
  library(readr)
  library(tidyr)

  file_list <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(file_list) == 0) stop("No CSV files found in the directory.")
  
  mass_intensity_list <- list()
  mass_elements_list <- list()
  all_elements <- c()

  for (file in file_list) {
    df <- read_csv(file, show_col_types = FALSE)
    file_name <- tools::file_path_sans_ext(basename(file))
    
    if (!all(c("Mass", "intensity") %in% colnames(df))) next

    df_mass_intensity <- df %>% select(Mass, intensity) %>% rename(!!file_name := intensity)
    mass_intensity_list[[file_name]] <- df_mass_intensity

    element_cols <- setdiff(colnames(df), c("Mass", "intensity"))
    all_elements <- union(all_elements, element_cols)

    if (length(element_cols) > 0) {  
      df_elements <- df %>% select(Mass, all_of(element_cols))
      mass_elements_list[[file_name]] <- df_elements
    }
  }

  if (length(mass_intensity_list) > 0) {
    merged_mass_intensity <- Reduce(function(x, y) full_join(x, y, by = "Mass"), mass_intensity_list) %>%
      mutate(across(everything(), ~ replace_na(.x, 0))) %>%
      distinct(Mass, .keep_all = TRUE)
    write_csv(merged_mass_intensity, output_mass_intensity)
  }

  if (length(mass_elements_list) > 0) {
    merged_mass_elements <- bind_rows(mass_elements_list) %>%
      mutate(across(all_of(all_elements), ~ replace_na(.x, 0))) %>%
      arrange(Mass) %>%
      distinct(Mass, .keep_all = TRUE)
    write_csv(merged_mass_elements, output_mass_elements)
  }
}
#' Merge Data Based on Molecular Formula (MolForm) with Filtering
#'
#' This function processes CSV files in a directory, filters MolForm duplicates based on the minimum absolute Error,
#' then merges intensity and element data for further analysis. It also saves filtered data for each sample.
#'
#' @param dir_path The directory containing CSV files.
#' @param output_molform_intensity Path to save the merged MolForm-Intensity file.
#' @param output_molform_elements Path to save the MolForm-Elements file.
#' @param output_filtered_molform Path to save the filtered MolForm file.
#' @param output_filtered_samples_dir Directory to save individual filtered sample files.
#'
#' @return NULL (files are saved to disk).
#' @export
merge_molform_intensity <- function(dir_path, output_molform_intensity, output_molform_elements, output_filtered_molform, output_filtered_samples_dir) {
  library(dplyr)
  library(readr)
  library(tidyr)

  file_list <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(file_list) == 0) stop("No CSV files found in the directory.")
  
  molform_intensity_list <- list()
  molform_elements_list <- list()
  all_elements <- c()
  filtered_data_list <- list()

  if (!dir.exists(output_filtered_samples_dir)) {
    dir.create(output_filtered_samples_dir, recursive = TRUE)
  }

  for (file in file_list) {
    df <- read_csv(file, show_col_types = FALSE)
    file_name <- tools::file_path_sans_ext(basename(file))
    
    if (!all(c("MolForm", "intensity") %in% colnames(df))) next
    
    if (any(duplicated(df$MolForm))) {
      df <- df %>%
        group_by(MolForm) %>%
        arrange(abs(Error), .by_group = TRUE) %>%
        slice(1) %>%
        ungroup()
    }
    
    filtered_data_list[[file_name]] <- df
    write_csv(df, file.path(output_filtered_samples_dir, paste0(file_name, "_filtered.csv")))
    
    df_molform_intensity <- df %>% select(MolForm, intensity) %>% rename(!!file_name := intensity)
    molform_intensity_list[[file_name]] <- df_molform_intensity

    element_cols <- setdiff(colnames(df), c("MolForm", "intensity"))
    all_elements <- union(all_elements, element_cols)

    if (length(element_cols) > 0) {  
      df_elements <- df %>% select(MolForm, all_of(element_cols))
      molform_elements_list[[file_name]] <- df_elements
    }
  }
  
  if (length(filtered_data_list) > 0) {
    filtered_data <- bind_rows(filtered_data_list) %>% distinct(MolForm, .keep_all = TRUE)
    write_csv(filtered_data, output_filtered_molform)
  }

  if (length(molform_intensity_list) > 0) {
    merged_molform_intensity <- Reduce(function(x, y) full_join(x, y, by = "MolForm"), molform_intensity_list) %>%
      mutate(across(everything(), ~ replace_na(.x, 0))) %>%
      distinct(MolForm, .keep_all = TRUE)
    write_csv(merged_molform_intensity, output_molform_intensity)
  }

  if (length(molform_elements_list) > 0) {
    merged_molform_elements <- bind_rows(molform_elements_list) %>%
      mutate(across(all_of(all_elements), ~ replace_na(.x, 0))) %>%
      arrange(MolForm) %>%
      distinct(MolForm, .keep_all = TRUE)
    write_csv(merged_molform_elements, output_molform_elements)
  }
}
