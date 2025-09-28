#' Complete biochemical transformation analysis pipeline
#'
#' @param data Data matrix with peaks as rows and samples as columns
#' @param mol Molecular information matrix with same row names as data
#' @param trans_db Transformation database with 'Name' and 'Mass' columns
#' @param error_term Mass error tolerance (default: 0.000010)
#' @param output_dir Output directory path (default: current working directory)
#' @param sample_name Dataset name for output files (default: "Dataset")
#' @param clustering_method Method for hierarchical clustering (default: "average")
#' @param build_tree Whether to build phylogenetic tree (default: TRUE)
#' @return List containing transformation results, network, and phylogenetic tree
#' @export


complete_transformation_analysis <- function(data, mol, trans_db, 
                                           error_term = 0.000010, 
                                           output_dir = ".", 
                                           sample_name = "Dataset",
                                           clustering_method = "average",
                                           build_tree = TRUE) {
    cat("[INFO] Starting complete transformation analysis for:", sample_name, "\n")
    
    cat("[STEP] Detecting transformations...\n")
    transformation_results <- detect_transformations(
        data = data,
        mol = mol,
        trans_db = trans_db,
        error_term = error_term,
        output_dir = output_dir,
        sample_name = sample_name
    )
    cat("[DONE] Transformations detected. Total matches:", nrow(transformation_results$peak_2_peak), "\n")
    
    tree_result <- NULL
    if (build_tree && nrow(transformation_results$peak_2_peak) > 0) {
        cat("[STEP] Building phylogenetic tree...\n")
        peak2peak_file <- file.path(output_dir, paste0(sample_name, "_All-Trans_peak.2.peak.csv"))
        numtrans_file <- file.path(output_dir, paste0(sample_name, "_All-Trans_num.peak.trans.csv"))
        
        tree_result <- build_phylogenetic_tree_from_files(
            peak2peak_file = peak2peak_file,
            numtrans_file = numtrans_file,
            sample_name = sample_name,
            output_dir = output_dir,
            clustering_method = clustering_method
        )
        cat("[DONE] Phylogenetic tree built and saved.\n")
    } else {
        cat("[INFO] No tree built (build_tree =", build_tree, ", or no transformations detected).\n")
    }
    
    cat("[INFO] Analysis finished for:", sample_name, "\n")
    return(list(
        transformations = transformation_results,
        tree_analysis = tree_result,
        sample_name = sample_name,
        parameters = list(
            error_term = error_term,
            clustering_method = clustering_method
        )
    ))
}

# ----------------------------------------------------------------------

# Detect transformations

detect_transformations <- function(data, mol, trans_db, 
                                   error_term = 0.000010, 
                                   output_dir = ".", 
                                   sample_name = "Dataset") {
    cat("[INFO] Checking input data...\n")
    if (!identical(row.names(data), row.names(mol))) {
        stop("Row names of data and mol matrices do not match")
    }
    if (!all(c("Name", "Mass") %in% colnames(trans_db))) {
        stop("Transformation database must contain 'Name' and 'Mass' columns")
    }
    
    cat("[INFO] Preparing peak and transformation data...\n")
    peak_masses <- as.numeric(row.names(mol))
    n_peaks <- length(peak_masses)
    trans_names <- as.character(trans_db$Name)
    trans_masses <- trans_db$Mass
    
    cat("[INFO] Total peaks:", n_peaks, " Total transformations in DB:", length(trans_names), "\n")
    
    result_list <- vector("list", n_peaks * (n_peaks - 1) / 2)
    result_count <- 0
    
    # Main loop
    pb <- txtProgressBar(min = 0, max = n_peaks - 1, style = 3)
    for (i in 1:(n_peaks - 1)) {
        current_peak <- peak_masses[i]
        remaining_peaks <- peak_masses[(i + 1):n_peaks]
        mass_diffs <- remaining_peaks - current_peak
        
        for (j in seq_along(trans_names)) {
            trans_mass <- trans_masses[j]
            matches <- abs(mass_diffs - trans_mass) <= error_term
            if (any(matches)) {
                matched_indices <- which(matches)
                result_count <- result_count + 1
                result_list[[result_count]] <- data.frame(
                    sample = "All",
                    peak.x = remaining_peaks[matched_indices],
                    peak.y = current_peak,
                    Dist = mass_diffs[matched_indices],
                    Dist.plus = mass_diffs[matched_indices] + error_term,
                    Dist.minus = mass_diffs[matched_indices] - error_term,
                    Trans.name = trans_names[j],
                    stringsAsFactors = FALSE
                )
            }
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    if (result_count > 0) {
        cat("[INFO] Transformation matches found:", result_count, "\n")
        peak_2_peak <- do.call(rbind, result_list[1:result_count])
        all_peaks <- c(peak_2_peak$peak.x, peak_2_peak$peak.y)
        peak_counts <- table(all_peaks)
        peak_profile <- data.frame(
            peak = as.numeric(names(peak_counts)),
            num_trans_involved_in = as.numeric(peak_counts),
            sample = sample_name,
            stringsAsFactors = FALSE
        )
    } else {
        cat("[WARN] No transformations detected.\n")
        peak_2_peak <- data.frame()
        peak_profile <- data.frame()
    }
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("[INFO] Created output directory:", output_dir, "\n")
    }
    
    write.csv(peak_2_peak, 
              file.path(output_dir, paste0(sample_name, "_All-Trans_peak.2.peak.csv")),
              quote = FALSE, row.names = FALSE)
    write.csv(peak_profile,
              file.path(output_dir, paste0(sample_name, "_All-Trans_num.peak.trans.csv")),
              quote = FALSE, row.names = FALSE)
    cat("[DONE] Transformation results saved to:", output_dir, "\n")
    
    return(list(
        peak_2_peak = peak_2_peak,
        peak_profile = peak_profile
    ))
}

# ----------------------------------------------------------------------

# Build phylogenetic tree

build_phylogenetic_tree_from_files <- function(peak2peak_file, numtrans_file, 
                                              sample_name, output_dir = ".",
                                              clustering_method = "average") {
    cat("[INFO] Loading transformation files for tree building...\n")
    if (!file.exists(peak2peak_file) || !file.exists(numtrans_file)) {
        stop("Required transformation files not found")
    }
    
    peak.2.peak <- read.csv(peak2peak_file, stringsAsFactors = FALSE)
    num.trans <- read.csv(numtrans_file, stringsAsFactors = FALSE)
    cat("[INFO] Loaded:", nrow(peak.2.peak), "edges and", nrow(num.trans), "vertices\n")
    
    peak.2.peak <- standardize_peak2peak_columns_base(peak.2.peak)
    num.trans <- standardize_numtrans_columns_base(num.trans)
    
    unnecessary_cols <- c('Dist', 'Dist.plus', 'Dist.minus')
    peak.2.peak <- peak.2.peak[, !names(peak.2.peak) %in% unnecessary_cols, drop = FALSE]
    peak.2.peak$weight <- 1
    
    edges <- peak.2.peak[, c("from", "to", "type", "weight", "sample")]
    vertices <- num.trans[, c("id", "num.trans.involved.in", "sample")]
    
    cat("[STEP] Constructing network...\n")
    net <- igraph::graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)
    clus <- igraph::components(net)
    max.clus.id <- which.max(clus$csize)
    max.clus.members <- names(clus$membership)[clus$membership == max.clus.id]
    cat("[INFO] Network built with", clus$no, "clusters. Largest has size", max(clus$csize), "\n")
    
    net.dist <- igraph::distances(net)
    net.dist <- net.dist[max.clus.members, max.clus.members]
    min_dist <- min(net.dist)
    max_dist <- max(net.dist)
    if (max_dist > min_dist) {
        net.dist <- (net.dist - min_dist) / (max_dist - min_dist)
    }
    
    cat("[STEP] Performing hierarchical clustering (method:", clustering_method, ")...\n")
    tree <- ape::as.phylo(stats::hclust(stats::as.dist(net.dist), method = clustering_method))
    tree_file <- file.path(output_dir, paste0(sample_name, "_TD_UPGMA.tre"))
    ape::write.tree(tree, tree_file)
    cat("[DONE] Tree saved to:", tree_file, "\n")
    
    return(list(
        tree = tree,
        network = net,
        distance_matrix = net.dist,
        cluster_info = list(
            total_clusters = clus$no,
            largest_cluster_size = max(clus$csize),
            largest_cluster_members = max.clus.members
        )
    ))
}

# ----------------------------------------------------------------------

standardize_peak2peak_columns_base <- function(df) {
    col_mapping <- list(
        'peak' = 'id',
        'peak.x' = 'from', 
        'peak.y' = 'to',
        'Trans.name' = 'type'
    )
    for (old_name in names(col_mapping)) {
        if (old_name %in% colnames(df)) {
            colnames(df)[colnames(df) == old_name] <- col_mapping[[old_name]]
        }
    }
    return(df)
}

standardize_numtrans_columns_base <- function(df) {
    col_mapping <- list(
        'peak' = 'id',
        'num_trans_involved_in' = 'num.trans.involved.in'
    )
    for (old_name in names(col_mapping)) {
        if (old_name %in% colnames(df)) {
            colnames(df)[colnames(df) == old_name] <- col_mapping[[old_name]]
        }
    }
    return(df)
}
