# --- Utility Function: Standardize Element Columns ---
#' Utility function to ensure required element columns exist in a dataframe.
#' Fills missing element columns with 0.
#' @param df Input data frame.
#' @param required_elements Character vector of elements required (e.g., c("C", "H", "O")).
#' @keywords internal
standardize_elements <- function(df, required_elements) {
  # All potential element columns used in matching
  element_cols_all <- c("C", "H", "O", "N", "S", "Cl", "Br", "P", "I")
  
  # Determine which of the element columns should be present based on required_elements
  elements_to_check <- base::intersect(element_cols_all, required_elements)
  
  # Find elements missing from the input dataframe
  missing_elements <- base::setdiff(elements_to_check, base::names(df))
  
  if (base::length(missing_elements) > 0) {
    for (el in missing_elements) {
      df[[el]] <- 0 # Add missing column and fill with 0
    }
    base::message(base::paste0("    [Info] Added missing element columns: ", base::paste(missing_elements, collapse = ", ")))
  }
  
  # Select only the necessary columns (Formula, intensity, and required elements)
  cols_to_select <- base::unique(c("Formula", "intensity", required_elements))
  cols_to_select <- cols_to_select[cols_to_select %in% base::names(df)]
  
  df <- df %>% dplyr::select(dplyr::all_of(cols_to_select))
  
  return(df)
}

# --- Reaction Matching Functions ---

#' @title Match reactions between two molecule datasets using database optimization
#'
#' @description This function filters precursor and product molecules based on intensity changes,
#' and matches them with a list of possible reaction deltas (elemental change) using SQLite for speed.
#'
#' @param file1 Path to the first molecular information CSV (e.g., inflow). Must contain 'MolForm', 'intensity', and elemental counts (C, H, O, N, S, etc.).
#' @param file2 Path to the second molecular information CSV (e.g., outflow). Must contain 'MolForm', 'intensity', and elemental counts (C, H, O, N, S, etc.).
#' @param reaction_delta_file Path to the reaction delta CSV file. Must contain 'reaction' and elemental delta columns (C, H, O, N, S, etc.).
#' @param out_dir Directory to save output files (default: ".").
#' @param use_memory_db Logical, whether to use in-memory database (default TRUE).
#'
#' @return Two CSV files saved in \code{out_dir}: \code{network_edge.csv} and \code{reaction_summary.csv}.
#' @export
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr select mutate rename filter bind_rows count arrange all_of
#' @importFrom RSQLite SQLite
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable dbExecute dbGetQuery
match_reactions_by_intensity <- function(file1, file2, reaction_delta_file, 
                                         out_dir = ".", use_memory_db = TRUE) {
  
  message("Step 1: Reading and standardizing molecular information files...")
  
  mol1 <- tryCatch(readr::read_csv(file1, show_col_types = FALSE), 
                   error = function(e) base::stop("Failed to read file1: ", e))
  mol2 <- tryCatch(readr::read_csv(file2, show_col_types = FALSE), 
                   error = function(e) base::stop("Failed to read file2: ", e))
  
  names(mol1) <- trimws(names(mol1))
  names(mol2) <- trimws(names(mol2))
  
  if (!"MolForm" %in% names(mol1) || !"MolForm" %in% names(mol2)) {
    base::stop("Error: 'MolForm' column is missing in input files.")
  }
  
  mol1 <- mol1 %>% dplyr::rename(Formula = MolForm)
  mol2 <- mol2 %>% dplyr::rename(Formula = MolForm)
  
  # Required element columns
  element_cols_all <- c("C", "H", "O", "N", "S", "Cl", "Br", "P", "I")
  final_element_cols <- element_cols_all # Element column names without suffix
  
  # Optimization 1: Standardize element columns (fill missing with 0)
  mol1 <- standardize_elements(mol1, element_cols_all)
  mol2 <- standardize_elements(mol2, element_cols_all)
  
  message("Step 2: Merging and filtering molecules by intensity ratio (FC < 0.5 for precursor, FC > 2.0 for product)...")
  
  # Dynamic columns needed for merge
  mol_cols_for_merge <- c("Formula", "intensity", element_cols_all)
  
  commom <- base::merge(mol1 %>% dplyr::select(dplyr::all_of(mol_cols_for_merge)), 
                        mol2 %>% dplyr::select(dplyr::all_of(mol_cols_for_merge)), 
                        by = "Formula", 
                        suffixes = c("_1", "_2")) %>%
    dplyr::rename(abundance_1 = intensity_1, abundance_2 = intensity_2)
  
  element_cols_1 <- paste0(element_cols_all, "_1")
  element_cols_2 <- paste0(element_cols_all, "_2")
  
  data <- commom %>%
    dplyr::select(Formula, abundance_1, dplyr::all_of(element_cols_1),
                  abundance_2, dplyr::all_of(element_cols_2)) %>%
    dplyr::mutate(
      pre = base::ifelse(abundance_2 / abundance_1 < 0.5, 1, 0),
      pro = base::ifelse(abundance_2 / abundance_1 > 2, 2, 0)
    )
  
  message("Step 3: Extracting precursor/product molecules...")
  
  # Extract precursors (decreasing intensity)
  data1 <- data %>% dplyr::filter(pre == 1) %>%
    dplyr::select(Formula, dplyr::all_of(element_cols_1))
  
  # Rename columns from _1 suffix to no suffix
  names(data1) <- c("Formula", final_element_cols)
  
  # Extract products (increasing intensity)
  data2 <- data %>% dplyr::filter(pro == 2) %>%
    dplyr::select(Formula, dplyr::all_of(element_cols_2))
  
  # Rename columns from _2 suffix to no suffix
  names(data2) <- c("Formula", final_element_cols)
  
  # Extract unique molecules
  unique_cols <- c("Formula", final_element_cols)
  unique1 <- mol1 %>% dplyr::filter(!(Formula %in% commom$Formula)) %>%
    dplyr::select(dplyr::all_of(unique_cols))
  unique2 <- mol2 %>% dplyr::filter(!(Formula %in% commom$Formula)) %>%
    dplyr::select(dplyr::all_of(unique_cols))
  
  mol1_filtered <- dplyr::bind_rows(unique1, data1)
  mol2_filtered <- dplyr::bind_rows(unique2, data2)
  
  message("Step 4: Reading reaction delta definitions...")
  
  reaction_delta <- tryCatch(readr::read_csv(reaction_delta_file, show_col_types = FALSE),
                             error = function(e) base::stop("Failed to read reaction delta file: ", e))
  names(reaction_delta) <- trimws(names(reaction_delta))
  
  # Optimization 2: Standardize element columns in delta file
  reaction_delta <- standardize_elements(reaction_delta, c("reaction", final_element_cols))
  
  # Ensure all element columns are numeric (for safety)
  reaction_delta[final_element_cols] <- base::lapply(reaction_delta[final_element_cols], function(x) {
    base::as.numeric(base::replace(x, base::is.na(x), 0))
  })
  
  mol1_filtered[final_element_cols] <- base::lapply(mol1_filtered[final_element_cols], function(x) {
    base::as.numeric(base::replace(x, base::is.na(x), 0))
  })
  mol2_filtered[final_element_cols] <- base::lapply(mol2_filtered[final_element_cols], function(x) {
    base::as.numeric(base::replace(x, base::is.na(x), 0))
  })
  
  message("Step 5: Creating database and loading data...")
  
  if (out_dir == ".") {
    db_path_base <- "temp_reactions"
  } else {
    db_path_base <- base::file.path(out_dir, "temp_reactions")
  }

  if (use_memory_db) {
    con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
    message("    Using in-memory database for maximum speed")
  } else {
    db_path <- paste0(db_path_base, ".db")
    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    message("    Using disk database: ", db_path)
  }
  
  base::on.exit(DBI::dbDisconnect(con), add = TRUE)
  
  DBI::dbWriteTable(con, "mol1", mol1_filtered, overwrite = TRUE)
  DBI::dbWriteTable(con, "mol2", mol2_filtered, overwrite = TRUE)
  DBI::dbWriteTable(con, "reactions", reaction_delta, overwrite = TRUE)
  
  message("Step 6: Creating indices for fast lookup...")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol1_formula ON mol1(Formula)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_formula ON mol2(Formula)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_elements ON mol2(C, H, N, O, S, Cl, Br, P, I)")
  
  message("Step 7: Matching reactions between molecules...")
  
  results_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = base::nrow(reaction_delta), style = 3)
  
  for (i in 1:base::nrow(reaction_delta)) {
    reaction_name <- reaction_delta$reaction[i]
    
    # Build SQL query (all element columns are guaranteed to exist now)
    query <- base::sprintf("
      SELECT 
        m1.Formula AS Source,
        m2.Formula AS Target,
        '%s' AS Reaction
      FROM mol1 m1
      INNER JOIN mol2 m2 ON
        ABS(m2.C - (m1.C + %f)) < 0.0001 AND
        ABS(m2.H - (m1.H + %f)) < 0.0001 AND
        ABS(m2.N - (m1.N + %f)) < 0.0001 AND
        ABS(m2.O - (m1.O + %f)) < 0.0001 AND
        ABS(m2.S - (m1.S + %f)) < 0.0001 AND
        ABS(m2.Cl - (m1.Cl + %f)) < 0.0001 AND
        ABS(m2.Br - (m1.Br + %f)) < 0.0001 AND
        ABS(m2.P - (m1.P + %f)) < 0.0001 AND
        ABS(m2.I - (m1.I + %f)) < 0.0001
    ", 
      reaction_name,
      reaction_delta$C[i],
      reaction_delta$H[i],
      reaction_delta$N[i],
      reaction_delta$O[i],
      reaction_delta$S[i],
      reaction_delta$Cl[i],
      reaction_delta$Br[i],
      reaction_delta$P[i],
      reaction_delta$I[i]
    )
    
    matches <- DBI::dbGetQuery(con, query)
    
    if (base::nrow(matches) > 0) {
      results_list[[i]] <- matches
    }
    
    utils::setTxtProgressBar(pb, i)
  }
  
  base::close(pb)
  
  if (base::length(results_list) > 0) {
    results <- dplyr::bind_rows(results_list)
  } else {
    results <- base::data.frame(Source = base::character(), Target = base::character(), 
                                Reaction = base::character(), stringsAsFactors = FALSE)
    base::warning("No reactions matched. Check your input data and delta definitions.")
  }
  
  message("Step 8: Saving outputs...")
  
  if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)
  
  out_file1 <- base::file.path(out_dir, "network_edge.csv")
  out_file2 <- base::file.path(out_dir, "reaction_summary.csv")
  
  readr::write_csv(results, out_file1)
  
  if (base::nrow(results) > 0) {
    reaction_summary <- results %>% 
      dplyr::count(Reaction, name = "Count") %>%
      dplyr::arrange(dplyr::desc(Count))
    readr::write_csv(reaction_summary, out_file2)
    
    message("Done! Found ", base::nrow(results), " reaction matches")
    message("    Results saved to:")
    message("    - ", out_file1)
    message("    - ", out_file2)
  } else {
    message("No matches found, output files created but empty")
  }
  
  return(base::invisible(results))
}
