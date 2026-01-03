#' Match reactions between two molecule datasets using database optimization
#'
#' This function filters precursor and product molecules based on intensity changes,
#' and matches them with a list of possible reaction deltas using SQLite for speed.
#'
#' @param file1 Path to the first molecular information CSV (e.g., inflow).
#' @param file2 Path to the second molecular information CSV (e.g., outflow).
#' @param reaction_delta_file Path to the reaction delta CSV file.
#' @param out_dir Directory to save output files.
#' @param use_memory_db Logical, whether to use in-memory database (default TRUE).
#'
#' @return Two CSV files saved in `out_dir`: `network_edge.csv` and `reaction_summary.csv`.
#' @export

# =========================================================
# Utility: Standardize element columns
# =========================================================
standardize_elements <- function(df, element_cols) {
  for (el in element_cols) {
    if (!el %in% names(df)) {
      df[[el]] <- 0
    }
  }
  df[element_cols] <- lapply(df[element_cols], function(x) {
    as.numeric(replace(x, is.na(x), 0))
  })
  return(df)
}

# =========================================================
# Main Function
# =========================================================
match_reactions_by_intensity <- function(file1,
                                         file2,
                                         reaction_delta_file,
                                         out_dir = ".",
                                         use_memory_db = TRUE) {

  # ---------------------------
  # Load packages
  # ---------------------------
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(DBI)
    library(RSQLite)
  })

  element_cols <- c("C", "H", "O", "N", "S", "Cl", "Br", "P", "I")

  message("Step 1: Reading and standardizing molecular information files...")

  mol1 <- tryCatch(
    readr::read_csv(file1, show_col_types = FALSE),
    error = function(e) base::stop("Failed to read file1: ", e)
  )
  mol2 <- tryCatch(
    readr::read_csv(file2, show_col_types = FALSE),
    error = function(e) base::stop("Failed to read file2: ", e)
  )

  names(mol1) <- trimws(names(mol1))
  names(mol2) <- trimws(names(mol2))

  if (!"MolForm" %in% names(mol1) || !"MolForm" %in% names(mol2)) {
    base::stop("Error: 'MolForm' column is missing in input files.")
  }

  mol1 <- mol1 %>% dplyr::rename(Formula = MolForm)
  mol2 <- mol2 %>% dplyr::rename(Formula = MolForm)

  # ---------------------------
  # Standardize molecule tables
  # ---------------------------
  mol1 <- standardize_elements(mol1, element_cols)
  mol2 <- standardize_elements(mol2, element_cols)

  if (!"intensity" %in% names(mol1) || !"intensity" %in% names(mol2)) {
    base::stop("Error: 'intensity' column is required in both input files.")
  }

  message("Step 2: Merging and filtering molecules by intensity ratio (FC < 0.5 for precursor, FC > 2.0 for product)...")

  # Dynamic columns needed for merge
  mol_cols_for_merge <- c("Formula", "intensity", element_cols)

  common <- base::merge(
    mol1 %>% dplyr::select(dplyr::all_of(mol_cols_for_merge)),
    mol2 %>% dplyr::select(dplyr::all_of(mol_cols_for_merge)),
    by = "Formula",
    suffixes = c("_1", "_2")
  ) %>%
    dplyr::rename(abundance_1 = intensity_1, abundance_2 = intensity_2)

  element_cols_1 <- paste0(element_cols, "_1")
  element_cols_2 <- paste0(element_cols, "_2")

  data <- common %>%
    dplyr::select(Formula, abundance_1, dplyr::all_of(element_cols_1),
                  abundance_2, dplyr::all_of(element_cols_2)) %>%
    dplyr::mutate(
      pre = base::ifelse(abundance_2 / abundance_1 < 0.5, 1, 0),
      pro = base::ifelse(abundance_2 / abundance_1 > 2.0, 2, 0)
    )

  message("Step 3: Extracting precursor and product molecules...")

  # Extract precursors (decreasing intensity)
  precursors <- data %>%
    dplyr::filter(pre == 1) %>%
    dplyr::select(Formula, dplyr::all_of(element_cols_1))
  names(precursors) <- c("Formula", element_cols)

  # Extract products (increasing intensity)
  products <- data %>%
    dplyr::filter(pro == 2) %>%
    dplyr::select(Formula, dplyr::all_of(element_cols_2))
  names(products) <- c("Formula", element_cols)

  # Extract unique molecules
  unique_cols <- c("Formula", element_cols)
  unique1 <- mol1 %>%
    dplyr::filter(!(Formula %in% common$Formula)) %>%
    dplyr::select(dplyr::all_of(unique_cols))
  unique2 <- mol2 %>%
    dplyr::filter(!(Formula %in% common$Formula)) %>%
    dplyr::select(dplyr::all_of(unique_cols))

  mol1_final <- dplyr::bind_rows(unique1, precursors)
  mol2_final <- dplyr::bind_rows(unique2, products)

  message("Step 4: Reading and standardizing reaction delta file...")

  reaction_delta <- tryCatch(
    readr::read_csv(reaction_delta_file, show_col_types = FALSE),
    error = function(e) base::stop("Failed to read reaction delta file: ", e)
  )
  names(reaction_delta) <- trimws(names(reaction_delta))

  if (!"reaction" %in% names(reaction_delta)) {
    base::stop("Error: 'reaction' column is required in reaction delta file.")
  }

  reaction_delta <- standardize_elements(reaction_delta, element_cols)

  # Ensure all element columns are numeric (for safety)
  reaction_delta[element_cols] <- base::lapply(reaction_delta[element_cols], function(x) {
    base::as.numeric(base::replace(x, base::is.na(x), 0))
  })

  mol1_final[element_cols] <- base::lapply(mol1_final[element_cols], function(x) {
    base::as.numeric(base::replace(x, base::is.na(x), 0))
  })

  mol2_final[element_cols] <- base::lapply(mol2_final[element_cols], function(x) {
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
    if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)
    con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    message("    Using disk database: ", db_path)
  }
  base::on.exit(DBI::dbDisconnect(con), add = TRUE)

  DBI::dbWriteTable(con, "mol1", mol1_final, overwrite = TRUE)
  DBI::dbWriteTable(con, "mol2", mol2_final, overwrite = TRUE)
  DBI::dbWriteTable(con, "reactions", reaction_delta, overwrite = TRUE)

  message("Step 6: Creating indices for fast lookup...")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol1_formula ON mol1(Formula)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_formula ON mol2(Formula)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_elements ON mol2(C, H, O, N, S, Cl, Br, P, I)")

  message("Step 7: Matching reactions between molecules...")

  results_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = base::nrow(reaction_delta), style = 3)

  for (i in seq_len(base::nrow(reaction_delta))) {

    delta <- reaction_delta[i, ]
    reaction_name <- delta$reaction

    query <- base::sprintf(
      "
      SELECT
        m1.Formula AS Source,
        m2.Formula AS Target,
        '%s' AS Reaction
      FROM mol1 m1
      INNER JOIN mol2 m2 ON
        ABS(m2.C  - (m1.C  + %f)) < 0.0001 AND
        ABS(m2.H  - (m1.H  + %f)) < 0.0001 AND
        ABS(m2.O  - (m1.O  + %f)) < 0.0001 AND
        ABS(m2.N  - (m1.N  + %f)) < 0.0001 AND
        ABS(m2.S  - (m1.S  + %f)) < 0.0001 AND
        ABS(m2.Cl - (m1.Cl + %f)) < 0.0001 AND
        ABS(m2.Br - (m1.Br + %f)) < 0.0001 AND
        ABS(m2.P  - (m1.P  + %f)) < 0.0001 AND
        ABS(m2.I  - (m1.I  + %f)) < 0.0001
      ",
      reaction_name,
      delta$C, delta$H, delta$O, delta$N,
      delta$S, delta$Cl, delta$Br, delta$P, delta$I
    )

    hits <- DBI::dbGetQuery(con, query)
    if (base::nrow(hits) > 0) {
      results_list[[i]] <- hits
    }

    utils::setTxtProgressBar(pb, i)
  }
  base::close(pb)

  if (base::length(results_list) > 0) {
    results <- dplyr::bind_rows(results_list)
  } else {
    results <- base::data.frame(
      Source = base::character(),
      Target = base::character(),
      Reaction = base::character(),
      stringsAsFactors = FALSE
    )
    base::warning("No reactions matched. Check your input data and delta definitions.")
  }

  message("Step 8: Saving outputs...")

  if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)

  out_file1 <- base::file.path(out_dir, "network_edge.csv")
  out_file2 <- base::file.path(out_dir, "reaction_summary.csv")

  readr::write_csv(results, out_file1)

  if (base::nrow(results) > 0) {
    summary <- results %>%
      dplyr::count(Reaction, name = "Count") %>%
      dplyr::arrange(dplyr::desc(Count))
    readr::write_csv(summary, out_file2)

    message("Done! Found ", base::nrow(results), " reaction matches")
    message("    Results saved to:")
    message("    - ", out_file1)
    message("    - ", out_file2)
  } else {
    message("No matches found, output files created but empty")
  }

  return(base::invisible(results))
}