#' Merge Mass and Intensity from Multiple CSV Files
#'
#' @param dir_path The directory containing CSV files.
#' @param output_mass_intensity Path to save the merged Mass-Intensity file.
#' @param output_mass_elements Path to save the Mass-Elements file.
#' @param output_meta_file Path to save the auto-generated metadata file (optional).
#'
#' @return A peakData object with the merged data.
#' @export
merge_mass_intensity <- function(dir_path, output_mass_intensity, output_mass_elements, output_meta_file = NULL) {
    install_dependencies <- function() {
        # Install devtools
        if (!requireNamespace("devtools", quietly = TRUE)) {
            install.packages("devtools")
        }
        
        # Install required packages
        required_packages <- c("readr", "dplyr", "tools", "ftmsRanalysis", "tidyr")
        for (pkg in required_packages) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                install.packages(pkg)
            }
        }
        
        # Install specific version of ftmsRanalysis
        if (!"ftmsRanalysis" %in% installed.packages()[, "Package"]) {
            devtools::install_github("EMSL-Computing/ftmsRanalysis@1.0.0")
        }
    }
    
    # Call dependency installation function
    install_dependencies()
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ftmsRanalysis)
    
    file_list <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
    if (length(file_list) == 0) stop("No CSV files found in the directory.")
    
    mass_intensity_list <- list()
    mass_elements_list <- list()
    all_elements <- c()
    
    # Create a list to store sample names for metadata generation
    sample_names <- character()
    
    for (file in file_list) {
        df <- read_csv(file, show_col_types = FALSE)
        file_name <- tools::file_path_sans_ext(basename(file))
        
        if (!all(c("Mass", "intensity") %in% colnames(df))) next
        
        # Add sample name to our list for metadata
        sample_names <- c(sample_names, file_name)
        
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
    } else {
        stop("No valid files with Mass and intensity columns found.")
    }
    
    if (length(mass_elements_list) > 0) {
        merged_mass_elements <- bind_rows(mass_elements_list) %>%
            mutate(across(all_of(all_elements), ~ replace_na(.x, 0))) %>%
            arrange(Mass) %>%
            distinct(Mass, .keep_all = TRUE)
        write_csv(merged_mass_elements, output_mass_elements)
    } else {
        stop("No element information found in the files.")
    }
    
    # Import data
    data_fticrms <- read.csv(output_mass_intensity) # Intensity data for each sample
    e_fticrms <- read.csv(output_mass_elements)    # Information on Mass Spectrometry
    
    # Critical fix: Extract actual sample IDs from the data frame columns
    # Get all column names except "Mass" - these are the sample IDs
    actual_sample_ids <- setdiff(colnames(data_fticrms), "Mass")
    
    # Generate metadata dataframe with correct sample IDs
    emeta <- data.frame(
        SampleID = actual_sample_ids,
        FileName = paste0(actual_sample_ids, ".csv"),
        SampleType = "Sample"  # Default sample type
    )
    
    # Save metadata file if path is provided
    if (!is.null(output_meta_file)) {
        write_csv(emeta, output_meta_file)
    }
    
    # Build e_fticrmsdata data frame
    columns_to_extract <- c("Mass", "NeutralMass", "Error", "C", "H", "O", "N", "C13", "S", "P", "Na")
    
    # Create base data frame with Mass from e_fticrms
    e_fticrmsdata <- e_fticrms %>% select(Mass)
    
    # Add existing columns from e_fticrms
    existing_columns <- intersect(columns_to_extract[-1], names(e_fticrms))
    if (length(existing_columns) > 0) {
        e_fticrmsdata <- e_fticrmsdata %>%
            left_join(e_fticrms %>% select(Mass, all_of(existing_columns)), by = "Mass")
    }
    
    # Add missing columns with zeros
    missing_columns <- setdiff(columns_to_extract[-1], names(e_fticrmsdata))
    for (col in missing_columns) {
        e_fticrmsdata[[col]] <- 0
    }
    
    # Ensure correct column order
    e_fticrmsdata <- e_fticrmsdata %>%
        select(all_of(columns_to_extract[columns_to_extract %in% names(e_fticrmsdata)]))
    
    # Create peakData object
    peakObj <- as.peakData(
        e_data = data_fticrms,
        f_data = emeta,
        e_meta = e_fticrmsdata,
        edata_cname = "Mass",
        mass_cname = "Mass",
        fdata_cname = "SampleID",
        c_cname = "C",
        h_cname = "H",
        o_cname = "O",
        n_cname = "N",
        s_cname = "S",
        p_cname = "P"
    )
    
    return(peakObj)
}

#' Merge Data Based on Molecular Formula (MolForm) with Filtering
#'
#' @param dir_path The directory containing CSV files.
#' @param output_molform_intensity Path to save the merged MolForm-Intensity file.
#' @param output_molform_elements Path to save the MolForm-Elements file.
#' @param output_filtered_samples_dir Directory to save individual filtered sample files.
#' @param output_meta_file Path to save the auto-generated metadata file (optional).
#'
#' @return A peakData object with the merged data.
#' @export
merge_molform_intensity <- function(dir_path, output_molform_intensity, output_molform_elements, 
                                    output_filtered_samples_dir, output_meta_file = NULL) {
    install_dependencies <- function() {
        # Install devtools
        if (!requireNamespace("devtools", quietly = TRUE)) {
            install.packages("devtools")
        }
        
        # Install required packages
        required_packages <- c("readr", "dplyr", "tools", "ftmsRanalysis", "tidyr")
        for (pkg in required_packages) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                install.packages(pkg)
            }
        }
        
        # Install specific version of ftmsRanalysis
        if (!"ftmsRanalysis" %in% installed.packages()[, "Package"]) {
            devtools::install_github("EMSL-Computing/ftmsRanalysis@1.0.0")
        }
    }
    
    # Call dependency installation function
    install_dependencies()
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ftmsRanalysis)
    
    file_list <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
    if (length(file_list) == 0) stop("No CSV files found in the directory.")
    
    molform_intensity_list <- list()
    molform_elements_list <- list()
    all_elements <- c()
    filtered_data_list <- list()
    
    # Create a list to store sample names for metadata generation
    sample_names <- character()
    
    if (!dir.exists(output_filtered_samples_dir)) {
        dir.create(output_filtered_samples_dir, recursive = TRUE)
    }
    
    for (file in file_list) {
        df <- read_csv(file, show_col_types = FALSE)
        file_name <- tools::file_path_sans_ext(basename(file))
        
        if (!all(c("MolForm", "intensity") %in% colnames(df))) next
        
        # Add sample name to our list for metadata
        sample_names <- c(sample_names, file_name)
        
        # Filter out duplicate MolForms, keeping the one with the lowest absolute Error
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
    
    if (length(molform_intensity_list) > 0) {
        merged_molform_intensity <- Reduce(function(x, y) full_join(x, y, by = "MolForm"), molform_intensity_list) %>%
            mutate(across(everything(), ~ replace_na(.x, 0))) %>%
            distinct(MolForm, .keep_all = TRUE)
        write_csv(merged_molform_intensity, output_molform_intensity)
    } else {
        stop("No valid files with MolForm and intensity columns found.")
    }
    
    if (length(molform_elements_list) > 0) {
        merged_molform_elements <- bind_rows(molform_elements_list) %>%
            mutate(across(all_of(all_elements), ~ replace_na(.x, 0))) %>%
            arrange(MolForm) %>%
            distinct(MolForm, .keep_all = TRUE)
        write_csv(merged_molform_elements, output_molform_elements)
    } else {
        stop("No element information found in the files.")
    }
    
    # Import data
    data_fticrms <- read.csv(output_molform_intensity) # Intensity data for each sample
    e_fticrms <- read.csv(output_molform_elements)    # Information on Mass Spectrometry
    
    # Critical fix: Extract actual sample IDs from the data frame columns
    # Get all column names except "MolForm" - these are the sample IDs
    actual_sample_ids <- setdiff(colnames(data_fticrms), "MolForm")
    
    # Generate metadata dataframe with correct sample IDs
    emeta <- data.frame(
        SampleID = actual_sample_ids,
        FileName = paste0(actual_sample_ids, ".csv"),
        SampleType = "Sample"  # Default sample type
    )
    
    # Save metadata file if path is provided
    if (!is.null(output_meta_file)) {
        write_csv(emeta, output_meta_file)
    }
    
    # Build e_fticrmsdata data frame
    columns_to_extract <- c("Mass", "NeutralMass", "Error", "C", "H", "O", "N", "C13", "S", "P", "Na")
    
    # Create base data frame with MolForm from e_fticrms
    e_fticrmsdata <- e_fticrms %>% select(MolForm)
    
    # Add existing columns from e_fticrms
    existing_columns <- intersect(columns_to_extract, names(e_fticrms))
    if (length(existing_columns) > 0) {
        e_fticrmsdata <- e_fticrmsdata %>%
            left_join(e_fticrms %>% select(MolForm, all_of(existing_columns)), by = "MolForm")
    }
    
    # Add missing columns with zeros
    missing_columns <- setdiff(columns_to_extract, names(e_fticrmsdata))
    for (col in missing_columns) {
        e_fticrmsdata[[col]] <- 0
    }
    
    # Ensure MolForm is the first column
    e_fticrmsdata <- e_fticrmsdata %>%
        select(MolForm, everything())
    
    # Create peakData object
    peakObj <- as.peakData(
        e_data = data_fticrms,
        f_data = emeta,
        e_meta = e_fticrmsdata,
        edata_cname = "MolForm",  # Use MolForm as identifier
        mass_cname = "Mass",
        fdata_cname = "SampleID",
        c_cname = "C",
        h_cname = "H",
        o_cname = "O",
        n_cname = "N",
        s_cname = "S",
        p_cname = "P"
    )
    
    return(peakObj)
}