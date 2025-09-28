#' Build molecular characteristics dendrogram (flexible feature extraction)
#'
#' @param mol_file Path to the molecular CSV file (*Mol.csv)
#' @param sample_name Name of the dataset (used for output)
#' @param output_dir Directory to save the resulting tree
#' @param clustering_method Clustering method for hclust (default: "average")
#' @return A list containing the phylogenetic tree, distance matrix, and processed molecular data
#' @export
build_molecular_dendrogram <- function(mol_file,
                                       sample_name = "Dataset",
                                       output_dir = ".",
                                       clustering_method = "average") {
  require(phangorn)
  require(ggtree)
  require(vegan)

  cat("[INFO] Loading molecular data...\n")
  mol <- read.csv(mol_file, row.names = 1, stringsAsFactors = FALSE)

  # Remove isotopic peaks if column exists
  if ("C13" %in% colnames(mol)) {
    w <- row.names(mol)[which(mol$C13 > 0)]
    if (length(w) > 0) {
      mol <- mol[-which(row.names(mol) %in% w), ]
    }
  }

  # Remove peaks with no molecular formula assignment
  mol <- mol[!is.na(mol$MolForm), ]

  # -------------------------------
  # Flexible feature extraction
  # -------------------------------
  info_cols <- c("C", "H", "O", "N", "S", "P", "DBE", "AI_Mod", "kdefect")
  ratio_cols <- c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio")

  Mol.Info <- mol[, intersect(info_cols, colnames(mol)), drop = FALSE]
  Mol.Ratio <- mol[, intersect(ratio_cols, colnames(mol)), drop = FALSE]

  cat("[INFO] Extracted", ncol(Mol.Info), "info columns and", ncol(Mol.Ratio), "ratio columns\n")

  # Ensure numeric
  Mol.Info <- as.data.frame(sapply(Mol.Info, as.numeric), row.names = row.names(Mol.Info))

  # Remove constant columns (zero variance)
  Mol.Info <- Mol.Info[, apply(Mol.Info, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

  # Scale molecular info
  Mol.Info <- as.data.frame(scale(Mol.Info), row.names = row.names(Mol.Info))

  # Remove rows with NA
  Mol.Info <- Mol.Info[complete.cases(Mol.Info), , drop = FALSE]

  if (nrow(Mol.Info) < 2) {
    stop("[ERROR] Not enough molecules to build dendrogram (need >= 2).")
  }

  # -------------------------------
  # Build dendrogram
  # -------------------------------
  mol.dist <- vegan::vegdist(Mol.Info, "euclidean")

  cat("[INFO] Building dendrogram using", clustering_method, "linkage...\n")
  tree <- ape::as.phylo(stats::hclust(mol.dist, method = clustering_method))

  # Save tree
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  tree_file <- file.path(output_dir, paste0(sample_name, "_MCD_UPGMA.tre"))
  ape::write.tree(tree, tree_file)
  cat("[DONE] Tree saved to:", tree_file, "\n")

  return(list(
    tree = tree,
    distance_matrix = mol.dist,
    mol = mol,
    Mol.Info = Mol.Info,
    Mol.Ratio = Mol.Ratio
  ))
}
