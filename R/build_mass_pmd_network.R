#' Build Mass PMD Network
#'
#' Creates a network based on paired mass differences between molecules.
#'
#' @param mol_file Character. Path to CSV file containing molecular information with 'Mass' column.
#' @param trans_file Character. Path to CSV file containing transformation database with 'Name' and 'Mass' columns.
#' @param error_term Numeric. Mass tolerance for matching transformations (default: 0.00001).
#' @param output_dir Character. Directory to save results (default: "mass_pmd_network").
#' @param mass_digits Integer. Number of decimal places for mass rounding (default: 6).
#'
#' @return Invisible NULL. Results are saved to files in the output directory.
#'
#' @details
#' This function reads molecular mass data and a transformation database, then identifies
#' potential biochemical transformations by matching mass differences. The results include:
#' \itemize{
#'   \item edges.csv - Network edges showing transformations between molecules
#'   \item network_stats.txt - Network topology statistics
#'   \item node_attributes.csv - Degree information for each mass
#'   \item transformation_summary.csv - Count of each transformation type
#'   \item sample_summary.csv - Overall dataset summary
#'   \item mass_verification.txt - Verification report for mass matching
#' }
#'
#' @examples
#' \dontrun{
#' build_mass_pmd_network(
#'     mol_file = "MS_MolInfor1.csv",
#'     trans_file = "Transformation_Database_07-2020.csv",
#'     error_term = 0.00001,
#'     output_dir = "MS_MolInfor2",
#'     mass_digits = 6
#' )
#' }
#'
#' @export
#' @importFrom dplyr filter mutate select count arrange desc bind_rows
#' @importFrom readr read_csv write_csv
#' @importFrom igraph graph_from_data_frame gorder gsize degree average.path.length diameter edge_density is.connected
#' @importFrom tools file_path_sans_ext
build_mass_pmd_network <- function(
        mol_file,
        trans_file,
        error_term = 0.00001,
        output_dir = "mass_pmd_network",
        mass_digits = 6
) {
    library(dplyr)
    library(readr)
    library(igraph)
    
    # Set precision for mass calculations
    options(digits = 15)
    
    # Create output directory
    if (!dir.exists(output_dir)) dir.create(output_dir)
    
    # Read input files
    mol_data <- read_csv(mol_file, show_col_types = FALSE)
    trans_db <- read_csv(trans_file, show_col_types = FALSE)
    trans_db$Name <- as.character(trans_db$Name)
    
    # Validate input columns
    if (!"Mass" %in% colnames(mol_data)) stop("mol_file must contain 'Mass' column.")
    if (!all(c("Name", "Mass") %in% colnames(trans_db))) stop("trans_file must contain 'Name' and 'Mass' columns.")
    
    # Ensure Mass is numeric and standardize precision
    mol_data$Mass <- as.numeric(mol_data$Mass)
    mol_data$Mass <- round(mol_data$Mass, mass_digits)
    
    # Store original masses for verification
    original_masses <- sort(unique(mol_data$Mass))
    
    # Generate all possible mass pairs and calculate PMD
    masses <- original_masses
    mass_pairs <- expand.grid(Mass.x = masses, Mass.y = masses, stringsAsFactors = FALSE) %>%
        filter(Mass.x > Mass.y) %>%
        mutate(
            Mass.x = round(Mass.x, mass_digits),
            Mass.y = round(Mass.y, mass_digits),
            PMD = round(Mass.x - Mass.y, mass_digits)
        )
    
    # Progress bar setup
    if (requireNamespace("progress", quietly = TRUE)) {
        pb <- progress::progress_bar$new(
            format = "Processing transformations [:bar] :percent (:current/:total) ETA: :eta",
            total = nrow(trans_db),
            clear = FALSE,
            width = 80
        )
        use_progress <- TRUE
    } else {
        message("Installing progress package for better user experience...")
        use_progress <- FALSE
    }
    
    # Match transformations
    trans_matched <- lapply(1:nrow(trans_db), function(i) {
        trans <- trans_db[i, ]
        name <- trans$Name
        mass_diff <- as.numeric(trans$Mass)
        mass_diff <- round(mass_diff, mass_digits)
        
        # Update progress
        if (use_progress) {
            pb$tick(tokens = list(current = i, total = nrow(trans_db)))
        } else {
            if (i %% max(1, floor(nrow(trans_db) / 10)) == 0 || i == nrow(trans_db)) {
                cat(sprintf("Progress: %d%% (%d/%d)\n", round(i/nrow(trans_db)*100), i, nrow(trans_db)))
            }
        }
        
        matches <- mass_pairs %>%
            filter(abs(PMD - mass_diff) <= error_term) %>%
            mutate(
                Trans_Name = name, 
                Trans_Mass = mass_diff,
                Mass.x = round(Mass.x, mass_digits),
                Mass.y = round(Mass.y, mass_digits)
            )
        return(matches)
    }) %>% bind_rows()
    
    # Check if any matches were found
    if (nrow(trans_matched) == 0) {
        warning("No matching transformation found.")
        return(invisible(NULL))
    }
    
    # Create edges data frame with standardized precision
    edges <- trans_matched %>%
        mutate(
            Source = round(Mass.x, mass_digits),
            Target = round(Mass.y, mass_digits),
            PMD = round(PMD, mass_digits)
        ) %>%
        select(Source, Target, PMD, Trans_Name)
    
    # Verify that all Source and Target masses exist in original data
    source_masses <- unique(edges$Source)
    target_masses <- unique(edges$Target)
    all_edge_masses <- unique(c(source_masses, target_masses))
    
    # Find any masses in edges that are not in original data
    missing_in_original <- setdiff(all_edge_masses, original_masses)
    
    # Create verification report
    verification_report <- paste0(
        "=== Mass Verification Report ===\n\n",
        "Original unique masses: ", length(original_masses), "\n",
        "Unique masses in edges (Source): ", length(source_masses), "\n",
        "Unique masses in edges (Target): ", length(target_masses), "\n",
        "Total unique masses in network: ", length(all_edge_masses), "\n",
        "Masses in edges but not in original data: ", length(missing_in_original), "\n"
    )
    
    if (length(missing_in_original) > 0) {
        verification_report <- paste0(
            verification_report,
            "\n⚠️  WARNING: Some masses in edges not found in original data!\n",
            "Missing masses: ", paste(head(missing_in_original, 20), collapse = ", "), 
            if(length(missing_in_original) > 20) " ..." else "", "\n"
        )
    } else {
        verification_report <- paste0(
            verification_report,
            "\n✅ All masses in edges match original data.\n"
        )
    }
    
    # Write verification report
    writeLines(verification_report, file.path(output_dir, "mass_verification.txt"))
    cat(verification_report)
    
    # Save edges
    write_csv(edges, file.path(output_dir, "edges.csv"))
    
    # Build network graph
    g <- graph_from_data_frame(edges, directed = TRUE)
    
    # Calculate network statistics
    num_nodes <- gorder(g)
    num_edges <- gsize(g)
    avg_degree <- mean(degree(g))
    avg_path_length <- if (is.connected(g)) average.path.length(g, directed = TRUE) else NA
    diameter <- if (is.connected(g)) diameter(g, directed = TRUE) else NA
    density <- edge_density(g)
    
    network_stats <- data.frame(
        Metric = c("Number of Nodes", "Number of Edges", "Average Degree", 
                   "Average Path Length", "Diameter", "Density"),
        Value = c(num_nodes, num_edges, avg_degree, avg_path_length, diameter, density)
    )
    
    write.table(network_stats, file.path(output_dir, "network_stats.txt"), 
                row.names = FALSE, quote = FALSE, sep = "\t")
    
    # Node attributes - ensure mass format consistency
    node_degree <- degree(g, mode = "all")
    node_attr <- data.frame(
        Mass = as.numeric(names(node_degree)), 
        Num_Transformations = node_degree
    )
    node_attr$Mass <- round(node_attr$Mass, mass_digits)
    write_csv(node_attr, file.path(output_dir, "node_attributes.csv"))
    
    # Transformation summary
    trans_summary <- edges %>%
        count(Trans_Name, name = "Count") %>%
        arrange(desc(Count))
    write_csv(trans_summary, file.path(output_dir, "transformation_summary.csv"))
    
    # Sample summary
    sample_summary <- data.frame(
        Sample = tools::file_path_sans_ext(basename(mol_file)),
        Total_Transformations = nrow(edges),
        Unique_Masses_Original = length(original_masses),
        Unique_Masses_In_Network = length(all_edge_masses),
        Mass_Match_Status = ifelse(length(missing_in_original) == 0, "All matched", "Some unmatched")
    )
    write_csv(sample_summary, file.path(output_dir, "sample_summary.csv"))
    
    message("\n✅ Transformation network analysis completed and saved to: ", output_dir)
    
    return(invisible(NULL))
}