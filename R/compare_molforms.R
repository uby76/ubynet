# R/compare_csv.R

#' Compare MolForm columns between two CSV files
#'
#' This function compares the MolForm columns between two CSV files and identifies
#' the intersection and differences.
#'
#' @param file1 Character. Path to the first CSV file.
#' @param file2 Character. Path to the second CSV file.
#' @param output_prefix Character. Prefix for output files. If NULL, uses the combined file names.
#' @param output_dir Character. Directory to save output files. If NULL, uses current directory.
#'
#' @return A list containing the comparison results with counts and sets of MolForm values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compare_molforms("MS_MolInfor1.csv", "MS_MolInfor2.csv", output_dir = "results")
#' }
compare_molforms <- function(file1, file2, output_prefix = NULL, output_dir = NULL) {
    # If no output prefix is specified, combine file names
    if(is.null(output_prefix)) {
        f1_name <- gsub("\\.csv$", "", basename(file1))
        f2_name <- gsub("\\.csv$", "", basename(file2))
        output_prefix <- paste0(f1_name, "_vs_", f2_name)
    }
    
    # Create output directory if provided
    if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        # Combine output prefix with directory
        output_prefix <- file.path(output_dir, output_prefix)
    }
    
    # Read CSV files
    data1 <- read.csv(file1, stringsAsFactors = FALSE)
    data2 <- read.csv(file2, stringsAsFactors = FALSE)
    
    # Extract MolForm columns for comparison
    molforms1 <- data1$MolForm
    molforms2 <- data2$MolForm
    
    # Calculate intersection and unique parts
    intersection_molforms <- intersect(molforms1, molforms2)
    only_in_data1_molforms <- setdiff(molforms1, molforms2)
    only_in_data2_molforms <- setdiff(molforms2, molforms1)
    
    # Prepare results
    result <- list(
        intersection_count = length(intersection_molforms),
        only_in_file1_count = length(only_in_data1_molforms),
        only_in_file2_count = length(only_in_data2_molforms),
        intersection = intersection_molforms,
        only_in_file1 = only_in_data1_molforms,
        only_in_file2 = only_in_data2_molforms
    )
    
    # Get complete row data based on MolForm values
    only_in_data1_full <- data1[data1$MolForm %in% only_in_data1_molforms, ]
    only_in_data2_full <- data2[data2$MolForm %in% only_in_data2_molforms, ]
    intersection_full <- merge(data1, data2, by = "MolForm", suffixes = c("_file1", "_file2"))
    
    # Save results to files
    write.csv(data.frame(MolForm = intersection_molforms), 
              paste0(output_prefix, "_molform_intersection.csv"), row.names = FALSE)
    write.csv(only_in_data1_full, 
              paste0(output_prefix, "_molform_only_in_", basename(file1)), row.names = FALSE)
    write.csv(only_in_data2_full, 
              paste0(output_prefix, "_molform_only_in_", basename(file2)), row.names = FALSE)
    
    return(result)
}

#' Compare Mass columns between two CSV files
#'
#' This function compares the Mass columns between two CSV files with a specified
#' tolerance value to find matching masses.
#'
#' @param file1 Character. Path to the first CSV file.
#' @param file2 Character. Path to the second CSV file.
#' @param output_prefix Character. Prefix for output files. If NULL, uses the combined file names.
#' @param output_dir Character. Directory to save output files. If NULL, uses current directory.
#' @param mass_tolerance Numeric. The tolerance for matching masses. Default is 0.01.
#'
#' @return A list containing the comparison results with counts and sets of Mass values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compare_mass("MS_MolInfor1.csv", "MS_MolInfor2.csv", 
#'                         output_dir = "results", mass_tolerance = 0.01)
#' }
compare_mass <- function(file1, file2, output_prefix = NULL, output_dir = NULL, mass_tolerance = 0.01) {
    # If no output prefix is specified, combine file names
    if(is.null(output_prefix)) {
        f1_name <- gsub("\\.csv$", "", basename(file1))
        f2_name <- gsub("\\.csv$", "", basename(file2))
        output_prefix <- paste0(f1_name, "_vs_", f2_name)
    }
    
    # Create output directory if provided
    if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        # Combine output prefix with directory
        output_prefix <- file.path(output_dir, output_prefix)
    }
    
    # Read CSV files
    data1 <- read.csv(file1, stringsAsFactors = FALSE)
    data2 <- read.csv(file2, stringsAsFactors = FALSE)
    
    # Ensure Mass column exists
    if(!"Mass" %in% names(data1) || !"Mass" %in% names(data2)) {
        stop("Both files must contain a Mass column")
    }
    
    # Extract Mass columns
    mass1 <- data1$Mass
    mass2 <- data2$Mass
    
    # Function to find matches within mass_tolerance range
    find_matches <- function(mass_value, mass_list, tolerance) {
        abs(mass_value - mass_list) <= tolerance
    }
    
    # Find intersection (considering mass tolerance)
    matches_matrix <- sapply(mass1, function(m1) find_matches(m1, mass2, mass_tolerance))
    
    # Find matched mass values
    matched_mass1_indices <- which(rowSums(matches_matrix) > 0)
    matched_mass2_indices <- which(colSums(matches_matrix) > 0)
    
    # Get matched and unmatched Mass values
    matched_mass1 <- mass1[matched_mass1_indices]
    matched_mass2 <- mass2[matched_mass2_indices]
    only_in_data1_mass <- mass1[!mass1 %in% matched_mass1]
    only_in_data2_mass <- mass2[!mass2 %in% matched_mass2]
    
    # Get complete row data
    only_in_data1_full <- data1[!mass1 %in% matched_mass1, ]
    only_in_data2_full <- data2[!mass2 %in% matched_mass2, ]
    
    # Create mapping relationships for matched masses
    match_pairs <- list()
    for(i in 1:length(mass1)) {
        matches <- find_matches(mass1[i], mass2, mass_tolerance)
        if(any(matches)) {
            for(j in which(matches)) {
                match_pairs[[length(match_pairs) + 1]] <- list(
                    mass1 = mass1[i],
                    mass2 = mass2[j],
                    index1 = i,
                    index2 = j,
                    diff = abs(mass1[i] - mass2[j])
                )
            }
        }
    }
    
    # Create matching dataframe - only including Mass information
    if(length(match_pairs) > 0) {
        matches_df <- do.call(rbind, lapply(match_pairs, function(pair) {
            data.frame(
                Mass_file1 = mass1[pair$index1],
                Mass_file2 = mass2[pair$index2],
                Mass_diff = pair$diff
            )
        }))
    } else {
        matches_df <- data.frame()
    }
    
    # Prepare results
    result <- list(
        matching_count = length(unique(matched_mass1)),  # Use unique to prevent duplicate counting
        only_in_file1_count = length(only_in_data1_mass),
        only_in_file2_count = length(only_in_data2_mass),
        matching_pairs = match_pairs,
        only_in_file1 = only_in_data1_mass,
        only_in_file2 = only_in_data2_mass
    )
    
    # Save results to files
    write.csv(matches_df, paste0(output_prefix, "_mass_matches.csv"), row.names = FALSE)
    write.csv(only_in_data1_full, paste0(output_prefix, "_mass_only_in_", basename(file1)), row.names = FALSE)
    write.csv(only_in_data2_full, paste0(output_prefix, "_mass_only_in_", basename(file2)), row.names = FALSE)
    
    return(result)
}

#' Compare multiple CSV datasets
#'
#' This function compares multiple CSV files either by MolForm or Mass columns.
#'
#' @param file_list Character vector. Paths to the CSV files to compare.
#' @param output_dir Character. Directory to save output files. If NULL, uses current directory.
#' @param comparison_type Character. Type of comparison, either "molform" or "mass". Default is "molform".
#' @param mass_tolerance Numeric. The tolerance for matching masses when comparison_type is "mass". Default is 0.01.
#'
#' @return A list containing the comparison results for all file pairs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' files <- c("MS_MolInfor1.csv", "MS_MolInfor2.csv", "MS_MolInfor3.csv")
#' results <- compare_multiple_datasets(files, output_dir = "results", comparison_type = "molform")
#' }
compare_multiple_datasets <- function(file_list, output_dir = NULL, comparison_type = "molform", mass_tolerance = 0.01) {
    n_files <- length(file_list)
    if(n_files < 2) {
        stop("At least two files are required for comparison")
    }
    
    # Create output directory if provided
    if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
    }
    
    # Create storage list for results
    results <- list()
    
    # First read all files
    all_data <- list()
    
    for(i in 1:n_files) {
        file_name <- file_list[i]
        data <- read.csv(file_name, stringsAsFactors = FALSE)
        all_data[[i]] <- data
    }
    
    if(comparison_type == "molform") {
        # Extract MolForm columns
        all_molforms <- lapply(all_data, function(data) data$MolForm)
        
        # Find MolForms common to all samples
        common_to_all <- Reduce(intersect, all_molforms)
        message("\nMolForms common to all ", n_files, " files: ", length(common_to_all))
        
        # Save the list of MolForms common to all samples
        if(length(common_to_all) > 0) {
            common_file <- if(!is.null(output_dir)) file.path(output_dir, "common_molform_to_all_files.csv") else "common_molform_to_all_files.csv"
            write.csv(data.frame(MolForm = common_to_all), common_file, row.names = FALSE)
        }
        
        # Compare each pair of files
        for(i in 1:(n_files-1)) {
            for(j in (i+1):n_files) {
                file1 <- file_list[i]
                file2 <- file_list[j]
                
                # Create output prefix
                f1_name <- gsub("\\.csv$", "", basename(file1))
                f2_name <- gsub("\\.csv$", "", basename(file2))
                output_prefix <- paste0(f1_name, "_vs_", f2_name)
                
                message("\nComparing MolForm between ", basename(file1), " and ", basename(file2), ":")
                result <- compare_molforms(file1, file2, output_prefix, output_dir)
                
                # Print basic results
                message("  Intersection count: ", result$intersection_count)
                message("  Only in ", basename(file1), ": ", result$only_in_file1_count)
                message("  Only in ", basename(file2), ": ", result$only_in_file2_count)
                
                # Store results
                results[[paste0(output_prefix, "_molform")]] <- result
            }
        }
        
        # Create summary for all files
        summary_data <- data.frame(
            Comparison = names(results),
            Intersection = sapply(results, function(x) x$intersection_count),
            File1_Only = sapply(results, function(x) x$only_in_file1_count),
            File2_Only = sapply(results, function(x) x$only_in_file2_count)
        )
        
        summary_file <- if(!is.null(output_dir)) file.path(output_dir, "molform_comparison_summary.csv") else "molform_comparison_summary.csv"
        write.csv(summary_data, summary_file, row.names = FALSE)
        
        # Add MolForms common to all files to the results
        results$common_to_all <- common_to_all
        
    } else if(comparison_type == "mass") {
        # Verify that all files have a Mass column
        for(i in 1:n_files) {
            if(!"Mass" %in% names(all_data[[i]])) {
                stop(paste("File", basename(file_list[i]), "is missing the Mass column"))
            }
        }
        
        # Extract Mass columns
        all_masses <- lapply(all_data, function(data) data$Mass)
        
        # Compare each pair of files
        for(i in 1:(n_files-1)) {
            for(j in (i+1):n_files) {
                file1 <- file_list[i]
                file2 <- file_list[j]
                
                # Create output prefix
                f1_name <- gsub("\\.csv$", "", basename(file1))
                f2_name <- gsub("\\.csv$", "", basename(file2))
                output_prefix <- paste0(f1_name, "_vs_", f2_name)
                
                message("\nComparing Mass between ", basename(file1), " and ", basename(file2), " (tolerance: ", mass_tolerance, "):")
                result <- compare_mass(file1, file2, output_prefix, output_dir, mass_tolerance)
                
                # Print basic results
                message("  Matching masses count: ", result$matching_count)
                message("  Only in ", basename(file1), ": ", result$only_in_file1_count)
                message("  Only in ", basename(file2), ": ", result$only_in_file2_count)
                
                # Store results
                results[[paste0(output_prefix, "_mass")]] <- result
            }
        }
        
        # Create summary for all files
        summary_data <- data.frame(
            Comparison = names(results),
            Matching = sapply(results, function(x) x$matching_count),
            File1_Only = sapply(results, function(x) x$only_in_file1_count),
            File2_Only = sapply(results, function(x) x$only_in_file2_count)
        )
        
        summary_file <- if(!is.null(output_dir)) file.path(output_dir, "mass_comparison_summary.csv") else "mass_comparison_summary.csv"
        write.csv(summary_data, summary_file, row.names = FALSE)
    } else {
        stop("comparison_type must be either 'molform' or 'mass'")
    }
    
    return(results)

# R/compare_csv.R

#' Compare MolForm columns between two CSV files
#'
#' This function compares the MolForm columns between two CSV files and identifies
#' the intersection and differences.
#'
#' @param file1 Character. Path to the first CSV file.
#' @param file2 Character. Path to the second CSV file.
#' @param output_prefix Character. Prefix for output files. If NULL, uses the combined file names.
#' @param output_dir Character. Directory to save output files. If NULL, uses current directory.
#'
#' @return A list containing the comparison results with counts and sets of MolForm values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compare_molforms("MS_MolInfor1.csv", "MS_MolInfor2.csv", output_dir = "results")
#' }
compare_molforms <- function(file1, file2, output_prefix = NULL, output_dir = NULL) {
    # If no output prefix is specified, combine file names
    if(is.null(output_prefix)) {
        f1_name <- gsub("\\.csv$", "", basename(file1))
        f2_name <- gsub("\\.csv$", "", basename(file2))
        output_prefix <- paste0(f1_name, "_vs_", f2_name)
    }
    
    # Create output directory if provided
    if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        # Combine output prefix with directory
        output_prefix <- file.path(output_dir, output_prefix)
    }
    
    # Read CSV files
    data1 <- read.csv(file1, stringsAsFactors = FALSE)
    data2 <- read.csv(file2, stringsAsFactors = FALSE)
    
    # Extract MolForm columns for comparison
    molforms1 <- data1$MolForm
    molforms2 <- data2$MolForm
    
    # Calculate intersection and unique parts
    intersection_molforms <- intersect(molforms1, molforms2)
    only_in_data1_molforms <- setdiff(molforms1, molforms2)
    only_in_data2_molforms <- setdiff(molforms2, molforms1)
    
    # Prepare results
    result <- list(
        intersection_count = length(intersection_molforms),
        only_in_file1_count = length(only_in_data1_molforms),
        only_in_file2_count = length(only_in_data2_molforms),
        intersection = intersection_molforms,
        only_in_file1 = only_in_data1_molforms,
        only_in_file2 = only_in_data2_molforms
    )
    
    # Get complete row data based on MolForm values
    only_in_data1_full <- data1[data1$MolForm %in% only_in_data1_molforms, ]
    only_in_data2_full <- data2[data2$MolForm %in% only_in_data2_molforms, ]
    intersection_full <- merge(data1, data2, by = "MolForm", suffixes = c("_file1", "_file2"))
    
    # Save results to files
    write.csv(data.frame(MolForm = intersection_molforms), 
              paste0(output_prefix, "_molform_intersection.csv"), row.names = FALSE)
    write.csv(only_in_data1_full, 
              paste0(output_prefix, "_molform_only_in_", basename(file1)), row.names = FALSE)
    write.csv(only_in_data2_full, 
              paste0(output_prefix, "_molform_only_in_", basename(file2)), row.names = FALSE)
    
    return(result)
}

#' Compare Mass columns between two CSV files
#'
#' This function compares the Mass columns between two CSV files with a specified
#' tolerance value to find matching masses.
#'
#' @param file1 Character. Path to the first CSV file.
#' @param file2 Character. Path to the second CSV file.
#' @param output_prefix Character. Prefix for output files. If NULL, uses the combined file names.
#' @param output_dir Character. Directory to save output files. If NULL, uses current directory.
#' @param mass_tolerance Numeric. The tolerance for matching masses. Default is 0.01.
#'
#' @return A list containing the comparison results with counts and sets of Mass values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compare_mass("MS_MolInfor1.csv", "MS_MolInfor2.csv", 
#'                         output_dir = "results", mass_tolerance = 0.01)
#' }
compare_mass <- function(file1, file2, output_prefix = NULL, output_dir = NULL, mass_tolerance = 0.01) {
    # If no output prefix is specified, combine file names
    if(is.null(output_prefix)) {
        f1_name <- gsub("\\.csv$", "", basename(file1))
        f2_name <- gsub("\\.csv$", "", basename(file2))
        output_prefix <- paste0(f1_name, "_vs_", f2_name)
    }
    
    # Create output directory if provided
    if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        # Combine output prefix with directory
        output_prefix <- file.path(output_dir, output_prefix)
    }
    
    # Read CSV files
    data1 <- read.csv(file1, stringsAsFactors = FALSE)
    data2 <- read.csv(file2, stringsAsFactors = FALSE)
    
    # Ensure Mass column exists
    if(!"Mass" %in% names(data1) || !"Mass" %in% names(data2)) {
        stop("Both files must contain a Mass column")
    }
    
    # Extract Mass columns
    mass1 <- data1$Mass
    mass2 <- data2$Mass
    
    # Function to find matches within mass_tolerance range
    find_matches <- function(mass_value, mass_list, tolerance) {
        abs(mass_value - mass_list) <= tolerance
    }
    
    # Find intersection (considering mass tolerance)
    matches_matrix <- sapply(mass1, function(m1) find_matches(m1, mass2, mass_tolerance))
    
    # Find matched mass values
    matched_mass1_indices <- which(rowSums(matches_matrix) > 0)
    matched_mass2_indices <- which(colSums(matches_matrix) > 0)
    
    # Get matched and unmatched Mass values
    matched_mass1 <- mass1[matched_mass1_indices]
    matched_mass2 <- mass2[matched_mass2_indices]
    only_in_data1_mass <- mass1[!mass1 %in% matched_mass1]
    only_in_data2_mass <- mass2[!mass2 %in% matched_mass2]
    
    # Get complete row data
    only_in_data1_full <- data1[!mass1 %in% matched_mass1, ]
    only_in_data2_full <- data2[!mass2 %in% matched_mass2, ]
    
    # Create mapping relationships for matched masses
    match_pairs <- list()
    for(i in 1:length(mass1)) {
        matches <- find_matches(mass1[i], mass2, mass_tolerance)
        if(any(matches)) {
            for(j in which(matches)) {
                match_pairs[[length(match_pairs) + 1]] <- list(
                    mass1 = mass1[i],
                    mass2 = mass2[j],
                    index1 = i,
                    index2 = j,
                    diff = abs(mass1[i] - mass2[j])
                )
            }
        }
    }
    
    # Create matching dataframe - only including Mass information
    if(length(match_pairs) > 0) {
        matches_df <- do.call(rbind, lapply(match_pairs, function(pair) {
            data.frame(
                Mass_file1 = mass1[pair$index1],
                Mass_file2 = mass2[pair$index2],
                Mass_diff = pair$diff
            )
        }))
    } else {
        matches_df <- data.frame()
    }
    
    # Prepare results
    result <- list(
        matching_count = length(unique(matched_mass1)),  # Use unique to prevent duplicate counting
        only_in_file1_count = length(only_in_data1_mass),
        only_in_file2_count = length(only_in_data2_mass),
        matching_pairs = match_pairs,
        only_in_file1 = only_in_data1_mass,
        only_in_file2 = only_in_data2_mass
    )
    
    # Save results to files
    write.csv(matches_df, paste0(output_prefix, "_mass_matches.csv"), row.names = FALSE)
    write.csv(only_in_data1_full, paste0(output_prefix, "_mass_only_in_", basename(file1)), row.names = FALSE)
    write.csv(only_in_data2_full, paste0(output_prefix, "_mass_only_in_", basename(file2)), row.names = FALSE)
    
    return(result)
}

#' Compare multiple CSV datasets
#'
#' This function compares multiple CSV files either by MolForm or Mass columns.
#'
#' @param file_list Character vector. Paths to the CSV files to compare.
#' @param output_dir Character. Directory to save output files. If NULL, uses current directory.
#' @param comparison_type Character. Type of comparison, either "molform" or "mass". Default is "molform".
#' @param mass_tolerance Numeric. The tolerance for matching masses when comparison_type is "mass". Default is 0.01.
#'
#' @return A list containing the comparison results for all file pairs.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' files <- c("MS_MolInfor1.csv", "MS_MolInfor2.csv", "MS_MolInfor3.csv")
#' results <- compare_multiple_datasets(files, output_dir = "results", comparison_type = "molform")
#' }
compare_multiple_datasets <- function(file_list, output_dir = NULL, comparison_type = "molform", mass_tolerance = 0.01) {
    n_files <- length(file_list)
    if(n_files < 2) {
        stop("At least two files are required for comparison")
    }
    
    # Create output directory if provided
    if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
    }
    
    # Create storage list for results
    results <- list()
    
    # First read all files
    all_data <- list()
    
    for(i in 1:n_files) {
        file_name <- file_list[i]
        data <- read.csv(file_name, stringsAsFactors = FALSE)
        all_data[[i]] <- data
    }
    
    if(comparison_type == "molform") {
        # Extract MolForm columns
        all_molforms <- lapply(all_data, function(data) data$MolForm)
        
        # Find MolForms common to all samples
        common_to_all <- Reduce(intersect, all_molforms)
        message("\nMolForms common to all ", n_files, " files: ", length(common_to_all))
        
        # Save the list of MolForms common to all samples
        if(length(common_to_all) > 0) {
            common_file <- if(!is.null(output_dir)) file.path(output_dir, "common_molform_to_all_files.csv") else "common_molform_to_all_files.csv"
            write.csv(data.frame(MolForm = common_to_all), common_file, row.names = FALSE)
        }
        
        # Compare each pair of files
        for(i in 1:(n_files-1)) {
            for(j in (i+1):n_files) {
                file1 <- file_list[i]
                file2 <- file_list[j]
                
                # Create output prefix
                f1_name <- gsub("\\.csv$", "", basename(file1))
                f2_name <- gsub("\\.csv$", "", basename(file2))
                output_prefix <- paste0(f1_name, "_vs_", f2_name)
                
                message("\nComparing MolForm between ", basename(file1), " and ", basename(file2), ":")
                result <- compare_molforms(file1, file2, output_prefix, output_dir)
                
                # Print basic results
                message("  Intersection count: ", result$intersection_count)
                message("  Only in ", basename(file1), ": ", result$only_in_file1_count)
                message("  Only in ", basename(file2), ": ", result$only_in_file2_count)
                
                # Store results
                results[[paste0(output_prefix, "_molform")]] <- result
            }
        }
        
        # Create summary for all files
        summary_data <- data.frame(
            Comparison = names(results),
            Intersection = sapply(results, function(x) x$intersection_count),
            File1_Only = sapply(results, function(x) x$only_in_file1_count),
            File2_Only = sapply(results, function(x) x$only_in_file2_count)
        )
        
        summary_file <- if(!is.null(output_dir)) file.path(output_dir, "molform_comparison_summary.csv") else "molform_comparison_summary.csv"
        write.csv(summary_data, summary_file, row.names = FALSE)
        
        # Add MolForms common to all files to the results
        results$common_to_all <- common_to_all
        
    } else if(comparison_type == "mass") {
        # Verify that all files have a Mass column
        for(i in 1:n_files) {
            if(!"Mass" %in% names(all_data[[i]])) {
                stop(paste("File", basename(file_list[i]), "is missing the Mass column"))
            }
        }
        
        # Extract Mass columns
        all_masses <- lapply(all_data, function(data) data$Mass)
        
        # Compare each pair of files
        for(i in 1:(n_files-1)) {
            for(j in (i+1):n_files) {
                file1 <- file_list[i]
                file2 <- file_list[j]
                
                # Create output prefix
                f1_name <- gsub("\\.csv$", "", basename(file1))
                f2_name <- gsub("\\.csv$", "", basename(file2))
                output_prefix <- paste0(f1_name, "_vs_", f2_name)
                
                message("\nComparing Mass between ", basename(file1), " and ", basename(file2), " (tolerance: ", mass_tolerance, "):")
                result <- compare_mass(file1, file2, output_prefix, output_dir, mass_tolerance)
                
                # Print basic results
                message("  Matching masses count: ", result$matching_count)
                message("  Only in ", basename(file1), ": ", result$only_in_file1_count)
                message("  Only in ", basename(file2), ": ", result$only_in_file2_count)
                
                # Store results
                results[[paste0(output_prefix, "_mass")]] <- result
            }
        }
        
        # Create summary for all files
        summary_data <- data.frame(
            Comparison = names(results),
            Matching = sapply(results, function(x) x$matching_count),
            File1_Only = sapply(results, function(x) x$only_in_file1_count),
            File2_Only = sapply(results, function(x) x$only_in_file2_count)
        )
        
        summary_file <- if(!is.null(output_dir)) file.path(output_dir, "mass_comparison_summary.csv") else "mass_comparison_summary.csv"
        write.csv(summary_data, summary_file, row.names = FALSE)
    } else {
        stop("comparison_type must be either 'molform' or 'mass'")
    }
    
    return(results)
    }
}