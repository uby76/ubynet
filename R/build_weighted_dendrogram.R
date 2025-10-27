#' Build transformation-weighted molecular characteristics dendrogram
#' 
#' This function generates a phylogenetic tree based on both molecular characteristics
#' and transformation network distances. It multiplies molecular distance (based on 
#' elemental composition and ratios) with transformation network distance to create
#' a weighted dendrogram, following the original RED 2020 methodology.
#'
#' @param mol_file Path to the molecular CSV file (*Mol.csv) or data.frame
#' @param peak2peak_file Path to the peak.2.peak CSV file (*peak.2.peak.csv) or data.frame
#' @param numtrans_file Path to the num.peak.trans CSV file (*num.peak.trans.csv) or data.frame
#' @param sample_name Name of the dataset (used for output files, default: "Dataset")
#' @param output_dir Directory to save the resulting tree (default: current directory)
#' @param clustering_method Clustering method for hclust (default: "average")
#' @param remove_isotopes Whether to remove C13 isotopic peaks (default: TRUE)
#' @param mol_features Vector of molecular features to use for distance calculation
#' @param ratio_features Vector of ratio features to use for distance calculation
#'
#' @return A list containing the phylogenetic tree, distance matrices, and processed data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Using file paths
#' result <- build_weighted_dendrogram(
#'   mol_file = "Dataset_Mol.csv",
#'   peak2peak_file = "Dataset_All-Trans_peak.2.peak.csv",
#'   numtrans_file = "Dataset_All-Trans_num.peak.trans.csv",
#'   sample_name = "MyDataset"
#' )
#' }
build_weighted_dendrogram <- function(mol_file, 
                                      peak2peak_file, 
                                      numtrans_file,
                                      sample_name = "Dataset",
                                      output_dir = ".",
                                      clustering_method = "average",
                                      remove_isotopes = TRUE,
                                      mol_features = c("C", "H", "N", "P", "O", "S", "DBE", "AI_Mod", "kdefect"),
                                      ratio_features = c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio")) {
    
    # Set options for precision (as in original code)
    options(digits = 10)
    
    # Load required libraries
    if (!requireNamespace("phangorn", quietly = TRUE)) {
        stop("Package 'phangorn' is required but not installed.")
    }
    if (!requireNamespace("ggtree", quietly = TRUE)) {
        stop("Package 'ggtree' is required but not installed.")
    }
    if (!requireNamespace("vegan", quietly = TRUE)) {
        stop("Package 'vegan' is required but not installed.")
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required but not installed.")
    }
    
    cat("[INFO] Starting transformation-weighted dendrogram analysis for:", sample_name, "\n")
    
    # ################## #
    #### Load in data ####
    # ################## #
    
    cat("[STEP] Loading input data...\n")
    
    # Load molecular data
    if (is.character(mol_file)) {
        if (!file.exists(mol_file)) {
            stop("Molecular file not found: ", mol_file)
        }
        mol <- read.csv(mol_file, row.names = 1, stringsAsFactors = FALSE)
    } else if (is.data.frame(mol_file)) {
        mol <- mol_file
    } else {
        stop("mol_file must be either a file path or a data.frame")
    }
    
    # Load peak2peak data
    if (is.character(peak2peak_file)) {
        if (!file.exists(peak2peak_file)) {
            stop("Peak2peak file not found: ", peak2peak_file)
        }
        peak.2.peak <- read.csv(peak2peak_file, stringsAsFactors = FALSE)
    } else if (is.data.frame(peak2peak_file)) {
        peak.2.peak <- peak2peak_file
    } else {
        stop("peak2peak_file must be either a file path or a data.frame")
    }
    
    # Load numtrans data
    if (is.character(numtrans_file)) {
        if (!file.exists(numtrans_file)) {
            stop("Numtrans file not found: ", numtrans_file)
        }
        num.trans <- read.csv(numtrans_file, stringsAsFactors = FALSE)
    } else if (is.data.frame(numtrans_file)) {
        num.trans <- numtrans_file
    } else {
        stop("numtrans_file must be either a file path or a data.frame")
    }
    
    cat("[INFO] Loaded molecular data:", nrow(mol), "peaks\n")
    cat("[INFO] Loaded transformation data:", nrow(peak.2.peak), "edges,", nrow(num.trans), "vertices\n")
    
    # Before doing anything, remove peaks that have an isotope signature
    if (remove_isotopes && "C13" %in% colnames(mol)) {
        isotope_indices <- which(mol$C13 > 0)
        if (length(isotope_indices) > 0) {
            mol <- mol[-isotope_indices, ]
            cat("[INFO] Removed", length(isotope_indices), "isotopic peaks\n")
        }
    }
    
    # Remove peaks that have no formula assignments
    if ("MolForm" %in% colnames(mol)) {
        na_indices <- which(is.na(mol$MolForm))
        if (length(na_indices) > 0) {
            mol <- mol[-na_indices, ]
            cat("[INFO] Removed", length(na_indices), "peaks without molecular formulas\n")
        }
    }
    
    # Setting objects for useful parameters (following original code structure)
    available_mol_features <- intersect(mol_features, colnames(mol))
    available_ratio_features <- intersect(ratio_features, colnames(mol))
    
    if (length(available_mol_features) > 0) {
        Mol.Info <- mol[, available_mol_features, drop = FALSE]
    } else {
        Mol.Info <- data.frame(row.names = row.names(mol))
    }
    
    if (length(available_ratio_features) > 0) {
        Mol.Rat <- mol[, available_ratio_features, drop = FALSE]
        # Combine Mol.Info and Mol.Rat
        if (ncol(Mol.Info) > 0) {
            Mol.Info <- cbind(Mol.Info, Mol.Rat)
        } else {
            Mol.Info <- Mol.Rat
        }
    }
    
    if (ncol(Mol.Info) == 0) {
        stop("No requested molecular features found in the data")
    }
    
    cat("[INFO] Using molecular features:", colnames(Mol.Info), "\n")
    
    # ####################### #
    #### Cleaning the data ####
    # ####################### #
    
    cat("[STEP] Cleaning transformation data...\n")
    
    # Altering transformation to match igraph structure (following original code)
    # peak.2.peak file
    colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak')] <- 'id'
    colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak.x')] <- 'from'
    colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'peak.y')] <- 'to'
    colnames(peak.2.peak)[which(colnames(peak.2.peak) == 'Trans.name')] <- 'type'
    peak.2.peak$weight <- 1
    
    # num.trans file
    colnames(num.trans)[which(colnames(num.trans) == 'peak')] <- 'id'
    colnames(num.trans)[which(colnames(num.trans) == 'num_trans_involved_in')] <- 'num.trans.involved.in'
    
    # Ensuring the two data sets match (following original code)
    missing_from <- length(which(!peak.2.peak$from %in% num.trans$id))
    missing_to <- length(which(!peak.2.peak$to %in% num.trans$id))
    cat("[INFO] Missing 'from' peaks:", missing_from, "\n")
    cat("[INFO] Missing 'to' peaks:", missing_to, "\n")
    
    # Reordering the data (following original code)
    peak.2.peak <- peak.2.peak[, c('from', 'to', 'type', 'weight', 'sample')]
    num.trans <- num.trans[, c('id', 'num.trans.involved.in', 'sample')]
    
    # ########################### #
    #### Determining distances ####
    # ########################### #
    
    cat("[STEP] Calculating molecular distances...\n")
    
    ### Pairwise molecular distance between peaks (following original code exactly)
    # Remove rows with any NA values before scaling
    complete_rows <- complete.cases(Mol.Info)
    if (!all(complete_rows)) {
        removed_count <- sum(!complete_rows)
        cat("[INFO] Removing", removed_count, "peaks with missing molecular data\n")
        Mol.Info <- Mol.Info[complete_rows, , drop = FALSE]
    }
    
    if (nrow(Mol.Info) < 2) {
        stop("Not enough peaks with complete molecular data (need >= 2)")
    }
    
    # Scale the data
    Mol.Info <- as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info))
    
    # Check for constant columns (zero variance) after scaling
    var_check <- apply(Mol.Info, 2, var, na.rm = TRUE)
    if (any(var_check == 0 | is.na(var_check))) {
        constant_cols <- names(var_check)[var_check == 0 | is.na(var_check)]
        cat("[INFO] Removing constant columns:", paste(constant_cols, collapse = ", "), "\n")
        Mol.Info <- Mol.Info[, !names(Mol.Info) %in% constant_cols, drop = FALSE]
    }
    
    if (ncol(Mol.Info) == 0) {
        stop("No valid molecular features remaining after removing constant columns")
    }
    
    mol.dist <- as.matrix(vegan::vegdist(Mol.Info, "euclidean"))
    
    cat("[STEP] Building transformation network...\n")
    
    ### Determining transformation distance (following original code)
    # Creating the network
    net <- igraph::graph_from_data_frame(d = peak.2.peak, vertices = num.trans, directed = FALSE)
    
    # The distances command is much better than the similarity measurement
    net.dist <- igraph::distances(net)
    
    # Finding clusters and determining the distance in the largest
    clus <- igraph::clusters(net)
    max.clus <- which(clus$csize %in% max(clus$csize)) # Finding the largest cluster
    max.clus <- names(clus$membership)[which(clus$membership %in% max.clus)] # Finding the members of the largest cluster
    
    cat("[INFO] Network has", clus$no, "clusters. Largest cluster has", length(max.clus), "peaks\n")
    
    net.dist <- net.dist[max.clus, max.clus] # Setting the net dist to that size only
    
    # Need to normalize the dissimilarity to 0-1 (following original code exactly)
    net.dist <- (net.dist - min(net.dist)) / (max(net.dist) - min(net.dist))
    
    ### Parsing down the data (following original code exactly)
    q <- which(row.names(net.dist) %in% row.names(mol.dist))
    net.dist.data <- net.dist[q, q] # Matching the network distance to the molecular information
    
    q <- which(row.names(mol.dist) %in% row.names(net.dist.data))
    mol.dist.data <- mol.dist[q, q] # Matching the molecular information to the network distance
    
    cat("[INFO] Final analysis includes", nrow(mol.dist.data), "common peaks\n")
    
    if (nrow(mol.dist.data) < 2) {
        stop("Not enough common peaks between molecular and network data (need >= 2)")
    }
    
    # Weighting and tree generation (following original code - MULTIPLICATION not addition!)
    weighted.dist <- mol.dist.data * net.dist.data
    
    # ######################### #
    #### Generating the tree ####
    # ######################### #
    
    cat("[STEP] Building phylogenetic tree using", clustering_method, "linkage...\n")
    
    # Creating tree (following original code exactly)
    tree <- ape::as.phylo(stats::hclust(stats::as.dist(weighted.dist), method = clustering_method))
    
    # Save results
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        cat("[INFO] Created output directory:", output_dir, "\n")
    }
    
    # Writing the tree (following original naming convention)
    tree_file <- file.path(output_dir, paste0(sample_name, "_Weighted_All-Trans_UPGMA.tre"))
    ape::write.tree(tree, tree_file)
    
    cat("[DONE] Analysis complete. Tree saved to:", tree_file, "\n")
    
    # Return results
    result <- list(
        tree = tree,
        distance_matrices = list(
            weighted = weighted.dist,
            molecular = mol.dist.data,
            network = net.dist.data
        ),
        network_info = list(
            network = net,
            clusters = clus,
            largest_cluster_size = length(max.clus),
            largest_cluster_members = max.clus
        ),
        molecular_data = list(
            original = mol,
            processed = Mol.Info,
            features_used = colnames(Mol.Info),
            common_peaks = row.names(mol.dist.data)
        ),
        parameters = list(
            sample_name = sample_name,
            clustering_method = clustering_method,
            remove_isotopes = remove_isotopes
        )
    )
    
    cat("[INFO] Transformation-weighted dendrogram analysis completed for:", sample_name, "\n")
    return(result)
}