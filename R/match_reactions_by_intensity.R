#' Match reactions between two molecule datasets using database optimization
#'
#' This function filters precursor and product molecules based on intensity changes,
#' and matches them with a list of possible reaction deltas using SQLite for speed.Refactor reaction matching logic with strict elemental standardization.
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

  element_cols <- c("C","H","O","N","S","Cl","Br","P","I")

  message("Step 1: Reading molecular information files...")

  mol1 <- readr::read_csv(file1, show_col_types = FALSE)
  mol2 <- readr::read_csv(file2, show_col_types = FALSE)

  names(mol1) <- trimws(names(mol1))
  names(mol2) <- trimws(names(mol2))

  if (!"MolForm" %in% names(mol1) || !"MolForm" %in% names(mol2)) {
    stop("MolForm column is required in both input files.")
  }

  mol1 <- mol1 %>% rename(Formula = MolForm)
  mol2 <- mol2 %>% rename(Formula = MolForm)

  # ---------------------------
  # Standardize molecule tables
  # ---------------------------
  mol1 <- standardize_elements(mol1, element_cols)
  mol2 <- standardize_elements(mol2, element_cols)

  if (!"intensity" %in% names(mol1) || !"intensity" %in% names(mol2)) {
    stop("intensity column is required in both input files.")
  }

  message("Step 2: Merging molecules and filtering by intensity change...")

  common <- merge(
    mol1,
    mol2,
    by = "Formula",
    suffixes = c("_1","_2")
  ) %>%
    rename(abundance_1 = intensity_1,
           abundance_2 = intensity_2)

  data <- common %>%
    mutate(
      pre = ifelse(abundance_2 / abundance_1 < 0.5, TRUE, FALSE),
      pro = ifelse(abundance_2 / abundance_1 > 2.0, TRUE, FALSE)
    )

  message("Step 3: Extracting precursor and product molecules...")

  # Precursors
  precursors <- data %>%
    filter(pre) %>%
    select(Formula, ends_with("_1"))
  names(precursors) <- c("Formula", element_cols)

  # Products
  products <- data %>%
    filter(pro) %>%
    select(Formula, ends_with("_2"))
  names(products) <- c("Formula", element_cols)

  # Unique molecules
  unique1 <- mol1 %>% filter(!Formula %in% common$Formula) %>%
    select(Formula, all_of(element_cols))
  unique2 <- mol2 %>% filter(!Formula %in% common$Formula) %>%
    select(Formula, all_of(element_cols))

  mol1_final <- bind_rows(unique1, precursors)
  mol2_final <- bind_rows(unique2, products)

  message("Step 4: Reading and standardizing reaction delta file...")

  reaction_delta <- readr::read_csv(reaction_delta_file, show_col_types = FALSE)
  names(reaction_delta) <- trimws(names(reaction_delta))

  if (!"reaction" %in% names(reaction_delta)) {
    stop("reaction column is required in reaction delta file.")
  }

  reaction_delta <- standardize_elements(reaction_delta, element_cols)

  message("Step 5: Creating SQLite database...")

  if (use_memory_db) {
    con <- dbConnect(SQLite(), ":memory:")
  } else {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    con <- dbConnect(SQLite(), file.path(out_dir, "temp_reactions.db"))
  }
  on.exit(dbDisconnect(con), add = TRUE)

  dbWriteTable(con, "mol1", mol1_final, overwrite = TRUE)
  dbWriteTable(con, "mol2", mol2_final, overwrite = TRUE)
  dbWriteTable(con, "reactions", reaction_delta, overwrite = TRUE)

  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol1_formula ON mol1(Formula)")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_formula ON mol2(Formula)")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_elements ON mol2(C,H,O,N,S,Cl,Br,P,I)")

  message("Step 6: Matching reactions...")

  results_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = nrow(reaction_delta), style = 3)

  for (i in seq_len(nrow(reaction_delta))) {

    delta <- reaction_delta[i, ]
    reaction_name <- delta$reaction

    query <- sprintf(
      "
      SELECT
        m1.Formula AS Source,
        m2.Formula AS Target,
        '%s' AS Reaction
      FROM mol1 m1
      INNER JOIN mol2 m2 ON
        ABS(m2.C  - (m1.C  + %f)) < 1e-4 AND
        ABS(m2.H  - (m1.H  + %f)) < 1e-4 AND
        ABS(m2.O  - (m1.O  + %f)) < 1e-4 AND
        ABS(m2.N  - (m1.N  + %f)) < 1e-4 AND
        ABS(m2.S  - (m1.S  + %f)) < 1e-4 AND
        ABS(m2.Cl - (m1.Cl + %f)) < 1e-4 AND
        ABS(m2.Br - (m1.Br + %f)) < 1e-4 AND
        ABS(m2.P  - (m1.P  + %f)) < 1e-4 AND
        ABS(m2.I  - (m1.I  + %f)) < 1e-4
      ",
      reaction_name,
      delta$C, delta$H, delta$O, delta$N,
      delta$S, delta$Cl, delta$Br, delta$P, delta$I
    )

    hits <- dbGetQuery(con, query)
    if (nrow(hits) > 0) {
      results_list[[length(results_list) + 1]] <- hits
    }

    utils::setTxtProgressBar(pb, i)
  }
  close(pb)

  results <- if (length(results_list) > 0) {
    bind_rows(results_list)
  } else {
    data.frame(Source = character(),
               Target = character(),
               Reaction = character(),
               stringsAsFactors = FALSE)
  }

  message("Step 7: Writing output files...")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  write_csv(results, file.path(out_dir, "network_edge.csv"))

  if (nrow(results) > 0) {
    summary <- results %>%
      count(Reaction, name = "Count") %>%
      arrange(desc(Count))
    write_csv(summary, file.path(out_dir, "reaction_summary.csv"))
  }

  message("Done. Total matched reactions: ", nrow(results))

  return(invisible(results))
}
