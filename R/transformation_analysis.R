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
    
    transformation_results <- detect_transformations(
        data = data,
        mol = mol,
        trans_db = trans_db,
        error_term = error_term,
        output_dir = output_dir,
        sample_name = sample_name
    )
    
    tree_result <- NULL
    if (build_tree && nrow(transformation_results$peak_2_peak) > 0) {
        peak2peak_file <- file.path(output_dir, paste0(sample_name, "_All-Trans_peak.2.peak.csv"))
        numtrans_file <- file.path(output_dir, paste0(sample_name, "_All-Trans_num.peak.trans.csv"))
        
        tree_result <- build_phylogenetic_tree_from_files(
            peak2peak_file = peak2peak_file,
            numtrans_file = numtrans_file,
            sample_name = sample_name,
            output_dir = output_dir,
            clustering_method = clustering_method
        )
    }
    
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

#' Detect putative biochemical transformations across dataset
#'
#' @param data Data matrix with peaks as rows and samples as columns
#' @param mol Molecular information matrix with same row names as data
#' @param trans_db Transformation database with 'Name' and 'Mass' columns
#' @param error_term Mass error tolerance (default: 0.000010)
#' @param output_dir Output directory path (default: current working directory)
#' @param sample_name Dataset name for output files (default: "Dataset")
#' @return List containing peak-to-peak transformations and peak profiles
#' @export
detect_transformations <- function(data, mol, trans_db, 
                                   error_term = 0.000010, 
                                   output_dir = ".", 
                                   sample_name = "Dataset") {
    
    if (!identical(row.names(data), row.names(mol))) {
        stop("Row names of data and mol matrices do not match")
    }
    
    if (!all(c("Name", "Mass") %in% colnames(trans_db))) {
        stop("Transformation database must contain 'Name' and 'Mass' columns")
    }
    
    peak_masses <- as.numeric(row.names(mol))
    n_peaks <- length(peak_masses)
    
    result_list <- vector("list", n_peaks * (n_peaks - 1) / 2)
    result_count <- 0
    
    trans_names <- as.character(trans_db$Name)
    trans_masses <- trans_db$Mass
    
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
    }
    
    if (result_count > 0) {
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
        peak_2_peak <- data.frame(
            sample = character(0),
            peak.x = numeric(0),
            peak.y = numeric(0),
            Dist = numeric(0),
            Dist.plus = numeric(0),
            Dist.minus = numeric(0),
            Trans.name = character(0),
            stringsAsFactors = FALSE
        )
        
        peak_profile <- data.frame(
            peak = numeric(0),
            num_trans_involved_in = numeric(0),
            sample = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    write.csv(peak_2_peak, 
              file.path(output_dir, paste0(sample_name, "_All-Trans_peak.2.peak.csv")),
              quote = FALSE, row.names = FALSE)
    
    write.csv(peak_profile,
              file.path(output_dir, paste0(sample_name, "_All-Trans_num.peak.trans.csv")),
              quote = FALSE, row.names = FALSE)
    
    return(list(
        peak_2_peak = peak_2_peak,
        peak_profile = peak_profile
    ))
}

#' Build phylogenetic tree from transformation files
#'
#' @param peak2peak_file Path to peak.2.peak CSV file
#' @param numtrans_file Path to num.peak.trans CSV file  
#' @param sample_name Sample name for output
#' @param output_dir Output directory
#' @param clustering_method Clustering method for tree construction
#' @return List containing tree and network analysis results
#' @export
build_phylogenetic_tree_from_files <- function(peak2peak_file, numtrans_file, 
                                              sample_name, output_dir = ".",
                                              clustering_method = "average") {
    
    # Load required libraries
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required but not installed. Please install it with: install.packages('ape')")
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required but not installed. Please install it with: install.packages('igraph')")
    }
    
    if (!file.exists(peak2peak_file) || !file.exists(numtrans_file)) {
        stop("Required transformation files not found")
    }
    
    # 使用基础R读取文件，避免data.table问题
    peak.2.peak <- read.csv(peak2peak_file, stringsAsFactors = FALSE)
    num.trans <- read.csv(numtrans_file, stringsAsFactors = FALSE)
    
    # 标准化列名
    peak.2.peak <- standardize_peak2peak_columns_base(peak.2.peak)
    num.trans <- standardize_numtrans_columns_base(num.trans)
    
    # 移除不必要的列并添加权重
    unnecessary_cols <- c('Dist', 'Dist.plus', 'Dist.minus')
    peak.2.peak <- peak.2.peak[, !names(peak.2.peak) %in% unnecessary_cols, drop = FALSE]
    peak.2.peak$weight <- 1
    
    # 选择需要的列
    edges <- peak.2.peak[, c("from", "to", "type", "weight", "sample")]
    vertices <- num.trans[, c("id", "num.trans.involved.in", "sample")]
    
    # 构建网络
    net <- igraph::graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)
    
    # 使用components()代替deprecated的clusters()
    clus <- igraph::components(net)
    max.clus.id <- which.max(clus$csize)
    max.clus.members <- names(clus$membership)[clus$membership == max.clus.id]
    
    # 计算距离矩阵
    net.dist <- igraph::distances(net)
    net.dist <- net.dist[max.clus.members, max.clus.members]
    
    # 标准化距离
    min_dist <- min(net.dist)
    max_dist <- max(net.dist)
    
    if (max_dist > min_dist) {
        net.dist <- (net.dist - min_dist) / (max_dist - min_dist)
    }
    
    # 构建系统发育树
    tree <- ape::as.phylo(stats::hclust(stats::as.dist(net.dist), method = clustering_method))
    
    # 保存树文件
    tree_file <- file.path(output_dir, paste0(sample_name, "_TD_UPGMA.tre"))
    ape::write.tree(tree, tree_file)
    
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

#' Analyze existing transformation files
#'
#' @param sample_name Sample name for output files
#' @param clustering_method Clustering method (default: "average") 
#' @param output_dir Output directory (default: current directory)
#' @return List containing tree analysis results
#' @export
analyze_existing_transformations <- function(sample_name = "Dataset_Name",
                                           clustering_method = "average",
                                           output_dir = ".") {
    
    peak2peak_files <- list.files(pattern = "*peak.2.peak.csv", full.names = TRUE)
    numtrans_files <- list.files(pattern = "*num.peak.trans.csv", full.names = TRUE)
    
    if (length(peak2peak_files) == 0 || length(numtrans_files) == 0) {
        stop("Required transformation files not found (*peak.2.peak.csv and *num.peak.trans.csv)")
    }
    
    if (length(peak2peak_files) > 1) {
        warning("Multiple peak.2.peak.csv files found, using first one")
    }
    if (length(numtrans_files) > 1) {
        warning("Multiple num.peak.trans.csv files found, using first one")
    }
    
    result <- build_phylogenetic_tree_from_files(
        peak2peak_file = peak2peak_files[1],
        numtrans_file = numtrans_files[1],
        sample_name = sample_name,
        output_dir = output_dir,
        clustering_method = clustering_method
    )
    
    return(result)
}

# Internal helper functions using base R (避免data.table问题)
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
    
    transformation_results <- detect_transformations(
        data = data,
        mol = mol,
        trans_db = trans_db,
        error_term = error_term,
        output_dir = output_dir,
        sample_name = sample_name
    )
    
    tree_result <- NULL
    if (build_tree && nrow(transformation_results$peak_2_peak) > 0) {
        peak2peak_file <- file.path(output_dir, paste0(sample_name, "_All-Trans_peak.2.peak.csv"))
        numtrans_file <- file.path(output_dir, paste0(sample_name, "_All-Trans_num.peak.trans.csv"))
        
        tree_result <- build_phylogenetic_tree_from_files(
            peak2peak_file = peak2peak_file,
            numtrans_file = numtrans_file,
            sample_name = sample_name,
            output_dir = output_dir,
            clustering_method = clustering_method
        )
    }
    
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

#' Detect putative biochemical transformations across dataset
#'
#' @param data Data matrix with peaks as rows and samples as columns
#' @param mol Molecular information matrix with same row names as data
#' @param trans_db Transformation database with 'Name' and 'Mass' columns
#' @param error_term Mass error tolerance (default: 0.000010)
#' @param output_dir Output directory path (default: current working directory)
#' @param sample_name Dataset name for output files (default: "Dataset")
#' @return List containing peak-to-peak transformations and peak profiles
#' @export
detect_transformations <- function(data, mol, trans_db, 
                                   error_term = 0.000010, 
                                   output_dir = ".", 
                                   sample_name = "Dataset") {
    
    if (!identical(row.names(data), row.names(mol))) {
        stop("Row names of data and mol matrices do not match")
    }
    
    if (!all(c("Name", "Mass") %in% colnames(trans_db))) {
        stop("Transformation database must contain 'Name' and 'Mass' columns")
    }
    
    peak_masses <- as.numeric(row.names(mol))
    n_peaks <- length(peak_masses)
    
    result_list <- vector("list", n_peaks * (n_peaks - 1) / 2)
    result_count <- 0
    
    trans_names <- as.character(trans_db$Name)
    trans_masses <- trans_db$Mass
    
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
    }
    
    if (result_count > 0) {
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
        peak_2_peak <- data.frame(
            sample = character(0),
            peak.x = numeric(0),
            peak.y = numeric(0),
            Dist = numeric(0),
            Dist.plus = numeric(0),
            Dist.minus = numeric(0),
            Trans.name = character(0),
            stringsAsFactors = FALSE
        )
        
        peak_profile <- data.frame(
            peak = numeric(0),
            num_trans_involved_in = numeric(0),
            sample = character(0),
            stringsAsFactors = FALSE
        )
    }
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    write.csv(peak_2_peak, 
              file.path(output_dir, paste0(sample_name, "_All-Trans_peak.2.peak.csv")),
              quote = FALSE, row.names = FALSE)
    
    write.csv(peak_profile,
              file.path(output_dir, paste0(sample_name, "_All-Trans_num.peak.trans.csv")),
              quote = FALSE, row.names = FALSE)
    
    return(list(
        peak_2_peak = peak_2_peak,
        peak_profile = peak_profile
    ))
}

#' Build phylogenetic tree from transformation files
#'
#' @param peak2peak_file Path to peak.2.peak CSV file
#' @param numtrans_file Path to num.peak.trans CSV file  
#' @param sample_name Sample name for output
#' @param output_dir Output directory
#' @param clustering_method Clustering method for tree construction
#' @return List containing tree and network analysis results
#' @export
build_phylogenetic_tree_from_files <- function(peak2peak_file, numtrans_file, 
                                              sample_name, output_dir = ".",
                                              clustering_method = "average") {
    
    # Load required libraries
    if (!requireNamespace("ape", quietly = TRUE)) {
        stop("Package 'ape' is required but not installed. Please install it with: install.packages('ape')")
    }
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required but not installed. Please install it with: install.packages('igraph')")
    }
    
    if (!file.exists(peak2peak_file) || !file.exists(numtrans_file)) {
        stop("Required transformation files not found")
    }
    
    # 使用基础R读取文件，避免data.table问题
    peak.2.peak <- read.csv(peak2peak_file, stringsAsFactors = FALSE)
    num.trans <- read.csv(numtrans_file, stringsAsFactors = FALSE)
    
    # 标准化列名
    peak.2.peak <- standardize_peak2peak_columns_base(peak.2.peak)
    num.trans <- standardize_numtrans_columns_base(num.trans)
    
    # 移除不必要的列并添加权重
    unnecessary_cols <- c('Dist', 'Dist.plus', 'Dist.minus')
    peak.2.peak <- peak.2.peak[, !names(peak.2.peak) %in% unnecessary_cols, drop = FALSE]
    peak.2.peak$weight <- 1
    
    # 选择需要的列
    edges <- peak.2.peak[, c("from", "to", "type", "weight", "sample")]
    vertices <- num.trans[, c("id", "num.trans.involved.in", "sample")]
    
    # 构建网络
    net <- igraph::graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)
    
    # 使用components()代替deprecated的clusters()
    clus <- igraph::components(net)
    max.clus.id <- which.max(clus$csize)
    max.clus.members <- names(clus$membership)[clus$membership == max.clus.id]
    
    # 计算距离矩阵
    net.dist <- igraph::distances(net)
    net.dist <- net.dist[max.clus.members, max.clus.members]
    
    # 标准化距离
    min_dist <- min(net.dist)
    max_dist <- max(net.dist)
    
    if (max_dist > min_dist) {
        net.dist <- (net.dist - min_dist) / (max_dist - min_dist)
    }
    
    # 构建系统发育树
    tree <- ape::as.phylo(stats::hclust(stats::as.dist(net.dist), method = clustering_method))
    
    # 保存树文件
    tree_file <- file.path(output_dir, paste0(sample_name, "_TD_UPGMA.tre"))
    ape::write.tree(tree, tree_file)
    
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

#' Analyze existing transformation files
#'
#' @param sample_name Sample name for output files
#' @param clustering_method Clustering method (default: "average") 
#' @param output_dir Output directory (default: current directory)
#' @return List containing tree analysis results
#' @export
analyze_existing_transformations <- function(sample_name = "Dataset_Name",
                                           clustering_method = "average",
                                           output_dir = ".") {
    
    peak2peak_files <- list.files(pattern = "*peak.2.peak.csv", full.names = TRUE)
    numtrans_files <- list.files(pattern = "*num.peak.trans.csv", full.names = TRUE)
    
    if (length(peak2peak_files) == 0 || length(numtrans_files) == 0) {
        stop("Required transformation files not found (*peak.2.peak.csv and *num.peak.trans.csv)")
    }
    
    if (length(peak2peak_files) > 1) {
        warning("Multiple peak.2.peak.csv files found, using first one")
    }
    if (length(numtrans_files) > 1) {
        warning("Multiple num.peak.trans.csv files found, using first one")
    }
    
    result <- build_phylogenetic_tree_from_files(
        peak2peak_file = peak2peak_files[1],
        numtrans_file = numtrans_files[1],
        sample_name = sample_name,
        output_dir = output_dir,
        clustering_method = clustering_method
    )
    
    return(result)
}

# Internal helper functions using base R (避免data.table问题)
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