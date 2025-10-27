#' Calculate CAR Matrix with Bidirectional Element Changes
#'
#' @description
#' This function performs CAR (Carbon Aromaticity Ratio) matrix calculation,
#' filters molecular formulas, and generates frequency statistics.
#' Enhanced version that calculates both directions of transformations since
#' mass spectrometry reactions are bidirectional.
#' Supports C, H, O, N, S, and P elements.
#'
#' @param csv_file Character. Path to the input CSV file containing molecular formulas
#' @param db_path Character. Path to the SQLite database file (default: "FTICR_CAR_smart.db")
#' @param car_min Numeric. Minimum CAR threshold (default: 0.45)
#' @param car_max Numeric. Maximum CAR threshold (default: 1.0)
#' @param track_direction Logical. Whether to track element addition (+) or removal (-) (default: TRUE)
#' @param small_molecules Character vector. Molecules to be filtered out (default: predefined list)
#' @param export_top Logical. Whether to export top formulas to CSV (default: TRUE)
#' @param top_n Integer. Number of top formulas to export (default: 50)
#' @param top_output_file Character. Output file name for top formulas (default: "Top50_formula.csv")
#' @param export_all_summary Logical. Whether to export all formula statistics (default: TRUE)
#' @param all_summary_file Character. Output file name for all statistics (default: "All_formula_summary.csv")
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item{db_connection: }{SQLite database connection object}
#'   \item{formula_summary: }{Data frame of formula frequency statistics}
#'   \item{top_n_formulas: }{Top N most frequent formulas}
#'   \item{stats: }{Summary statistics}
#' }
#'
#' @import DBI
#' @import RSQLite
#' @import data.table
#' @import dplyr
#' @import progress
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with bidirectional tracking
#' result <- calculate_car_matrix("merged_molform_elements.csv")
#' 
#' # Without directional tracking (absolute differences only)
#' result <- calculate_car_matrix(
#'   csv_file = "data.csv",
#'   track_direction = FALSE
#' )
#' }

calculate_car_matrix <- function(
    csv_file,
    db_path = "FTICR_CAR_smart.db",
    car_min = 0.45,
    car_max = 1.0,
    track_direction = TRUE,
    small_molecules = c("C", "O", "CH2", "H2O", "CO", "O2", "CO2", "NH3", "NH"),
    export_top = TRUE,
    top_n = 50,
    top_output_file = "Top50_formula.csv",
    export_all_summary = TRUE,
    all_summary_file = "All_formula_summary.csv",
    verbose = TRUE
) {
  
  # ========== 1. Check and Load Dependencies ==========
  required_packages <- c("DBI", "RSQLite", "progress", "data.table", "dplyr")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed. Please install it first."))
    }
  }
  
  if (verbose) cat("✅ Dependency check completed\n")
  
  if (verbose && track_direction) {
    cat("📐 Bidirectional tracking enabled: calculating both MolForm_1→2 and MolForm_2→1\n")
  }
  
  # ========== 2. Read CSV Data ==========
  if (!file.exists(csv_file)) {
    stop(paste("File not found:", csv_file))
  }
  
  # Read data and handle both with/without P column
  df_raw <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Check if P column exists, if not add it with 0s
  required_cols <- c("MolForm", "C", "H", "O", "N", "S", "P")
  available_cols <- intersect(required_cols, names(df_raw))
  
  df <- df_raw[, available_cols, drop = FALSE]
  
  # Add P column if missing
  if (!"P" %in% names(df)) {
    df$P <- 0
    if (verbose) cat("ℹ️ P column not found, added with default value 0\n")
  }
  
  # Reorder columns
  df <- df[, c("MolForm", "C", "H", "O", "N", "S", "P")]
  
  num_cols <- sapply(df, is.numeric)
  df[num_cols][is.na(df[num_cols])] <- 0
  
  if (verbose) cat("✅ Data reading completed,", nrow(df), "rows in total\n")
  
  # ========== 3. Create Database and Write Origin Table ==========
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  
  if (DBI::dbExistsTable(con, "origin")) DBI::dbRemoveTable(con, "origin")
  DBI::dbWriteTable(con, "origin", df)
  
  if (verbose) cat("✅ Origin table written\n")
  
  # ========== 4. Build Matrix ==========
  data <- DBI::dbReadTable(con, "origin")
  mat <- as.matrix(data[, c("C", "H", "O", "N", "S", "P")])
  rownames(mat) <- data$MolForm
  n <- nrow(mat)
  
  if (verbose) cat("✅ Matrix construction completed,", n, "rows in total\n")
  
  # ========== 5. Define Vectorized CAR Function ==========
  calc_CAR_vec <- function(C1, H1, O1, C2, H2, O2) {
    (pmin(C1, C2) / pmax(C1, C2)) *
      (pmin(H1, H2) / pmax(H1, H2)) *
      (pmin(O1, O2) / pmax(O1, O2))
  }
  
  # ========== 6. Generate Pairs ==========
  pairs <- combn(n, 2)
  total_pairs <- ncol(pairs)
  
  if (verbose) cat("⚙️ Row pairs generated,", total_pairs, "combinations in total\n")
  
  # ========== 7. Calculate CAR with Progress Bar ==========
  if (verbose) {
    pb <- progress::progress_bar$new(
      total = total_pairs,
      format = "⏳ [:bar] :percent Remaining: :eta",
      clear = FALSE, width = 70
    )
  }
  
  C1 <- mat[pairs[1, ], "C"]; H1 <- mat[pairs[1, ], "H"]; O1 <- mat[pairs[1, ], "O"]
  C2 <- mat[pairs[2, ], "C"]; H2 <- mat[pairs[2, ], "H"]; O2 <- mat[pairs[2, ], "O"]
  
  CAR <- calc_CAR_vec(C1, H1, O1, C2, H2, O2)
  CAR[is.na(CAR)] <- 0
  keep <- which(CAR >= car_min & CAR <= car_max)
  
  if (verbose) {
    pb$tick(total_pairs)
    cat("\n✅ Filtered", length(keep), "pairs meeting CAR ∈ [", car_min, ",", car_max, "]\n")
  }
  
  # ========== 8. Calculate Differences (Bidirectional) ==========
  A <- mat[pairs[1, keep], ]
  B <- mat[pairs[2, keep], ]
  
  if (track_direction) {
    # Direction 1: B - A (MolForm_1 to MolForm_2)
    diff_mat_1to2 <- B - A
    
    # Direction 2: A - B (MolForm_2 to MolForm_1)
    diff_mat_2to1 <- A - B
    
    # Create two sets of results for bidirectional
    results_1to2 <- data.frame(
      MolForm_1 = rownames(mat)[pairs[1, keep]],
      MolForm_2 = rownames(mat)[pairs[2, keep]],
      C = diff_mat_1to2[, "C"], 
      H = diff_mat_1to2[, "H"], 
      O = diff_mat_1to2[, "O"],
      N = diff_mat_1to2[, "N"],
      S = diff_mat_1to2[, "S"],
      P = diff_mat_1to2[, "P"],
      CAR = CAR[keep]
    )
    
    results_2to1 <- data.frame(
      MolForm_1 = rownames(mat)[pairs[2, keep]],
      MolForm_2 = rownames(mat)[pairs[1, keep]],
      C = diff_mat_2to1[, "C"], 
      H = diff_mat_2to1[, "H"], 
      O = diff_mat_2to1[, "O"],
      N = diff_mat_2to1[, "N"],
      S = diff_mat_2to1[, "S"],
      P = diff_mat_2to1[, "P"],
      CAR = CAR[keep]
    )
    
    # Combine both directions
    results <- rbind(results_1to2, results_2to1)
    
  } else {
    # Absolute differences only
    diff_mat <- abs(A - B)
    
    results <- data.frame(
      MolForm_1 = rownames(mat)[pairs[1, keep]],
      MolForm_2 = rownames(mat)[pairs[2, keep]],
      C = diff_mat[, "C"], 
      H = diff_mat[, "H"], 
      O = diff_mat[, "O"],
      N = diff_mat[, "N"],
      S = diff_mat[, "S"],
      P = diff_mat[, "P"],
      CAR = CAR[keep]
    )
  }
  
  if (verbose) {
    if (track_direction) {
      cat("🧩 Bidirectional difference matrix generated,", nrow(results), 
          "records in total (", length(keep), "pairs × 2 directions )\n")
    } else {
      cat("🧩 Difference matrix generated,", nrow(results), "records in total\n")
    }
  }
  
  # ========== 9. Write car_matrix Table ==========
  if (DBI::dbExistsTable(con, "car_matrix")) DBI::dbRemoveTable(con, "car_matrix")
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "car_matrix", results)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ car_matrix table written\n")
  
  # ========== 10. Molecular Formula Composition Function (with P) ==========
  compose_formula <- function(C, H, O, N, S, P, track_dir = track_direction) {
    vals <- as.numeric(c(C, H, O, N, S, P))
    names(vals) <- c("C", "H", "O", "N", "S", "P")
    vals[is.na(vals)] <- 0
    
    if (!track_dir) {
      # Absolute mode: use absolute values
      vals <- abs(vals)
      vals <- round(vals)
      nonzero <- vals != 0
      if (!any(nonzero)) return("")
      paste0(names(vals)[nonzero],
             ifelse(vals[nonzero] == 1, "", vals[nonzero]),
             collapse = "")
    } else {
      # Directional mode: keep signs
      vals <- round(vals)
      nonzero <- vals != 0
      if (!any(nonzero)) return("")
      
      # Format with + or - prefix
      formatted <- sapply(which(nonzero), function(i) {
        val <- vals[i]
        elem <- names(vals)[i]
        sign <- ifelse(val > 0, "+", "-")
        abs_val <- abs(val)
        count_str <- ifelse(abs_val == 1, "", abs_val)
        paste0(sign, elem, count_str)
      })
      paste0(formatted, collapse = "")
    }
  }
  
  # ========== 11. Auto-select Processing Mode ==========
  row_n <- nrow(results)
  mode <- if (row_n <= 50000) "mapply" else "data.table"
  
  if (verbose) cat("🧠 Auto-selected mode:", mode, "(Data rows:", row_n, ")\n")
  
  # ========== 12. Generate Standardized Molecular Formulas ==========
  if (verbose) {
    pb2 <- progress::progress_bar$new(
      total = row_n,
      format = "🔬 [:bar] :percent Remaining: :eta",
      clear = FALSE, width = 70
    )
    tick_every <- ceiling(row_n / 100)
    counter <- 0
  }
  
  if (mode == "mapply") {
    results$Diff_Formula <- mapply(function(C, H, O, N, S, P) {
      if (verbose) {
        counter <<- counter + 1
        if (counter %% tick_every == 0) pb2$tick(tick_every)
      }
      compose_formula(C, H, O, N, S, P)
    }, results$C, results$H, results$O, results$N, results$S, results$P)
  } else {
    data.table::setDT(results)
    el_names <- c("C", "H", "O", "N", "S", "P")
    results[, Diff_Formula := apply(.SD, 1, function(x) {
      if (verbose) {
        counter <<- counter + 1
        if (counter %% tick_every == 0) pb2$tick(tick_every)
      }
      compose_formula(x[1], x[2], x[3], x[4], x[5], x[6])
    }), .SDcols = el_names]
  }
  
  if (verbose) {
    pb2$tick(row_n)
    if (track_direction) {
      cat("\n✅ Directional molecular formulas generated (e.g., +C2-H2O)\n")
    } else {
      cat("\n✅ Standardized molecular formulas generated\n")
    }
  }
  
  # ========== 13. Write car_matrix_formula Table ==========
  if (DBI::dbExistsTable(con, "car_matrix_formula")) DBI::dbRemoveTable(con, "car_matrix_formula")
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "car_matrix_formula", results)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ car_matrix_formula table written (containing all molecular formulas)\n")
  
  # ========== 14. Filter Small Molecules ==========
  # For directional mode, also filter ±small molecules
  if (track_direction) {
    small_patterns <- c(small_molecules,
                        paste0("+", small_molecules),
                        paste0("-", small_molecules))
  } else {
    small_patterns <- small_molecules
  }
  
  small_df <- subset(results, Diff_Formula %in% small_patterns)
  filtered_df <- subset(results, !(Diff_Formula %in% small_patterns))
  
  if (verbose) cat("🚫 Small molecule filtering:", nrow(small_df), "small molecules,", nrow(filtered_df), "remaining\n")
  
  # ========== 15. Write Small Molecules Table ==========
  if (DBI::dbExistsTable(con, "car_matrix_smallmolecules")) {
    DBI::dbRemoveTable(con, "car_matrix_smallmolecules")
  }
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "car_matrix_smallmolecules", small_df)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ Small molecules table car_matrix_smallmolecules written\n")
  
  # ========== 16. Write Filtered Data Table ==========
  if (DBI::dbExistsTable(con, "car_matrix_filtered")) {
    DBI::dbRemoveTable(con, "car_matrix_filtered")
  }
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "car_matrix_filtered", filtered_df)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ Filtered table car_matrix_filtered written\n")
  
  # ========== 17. Frequency Statistics with Element Details ==========
  if (verbose) cat("📊 Performing Diff_Formula frequency statistics (small molecules excluded)...\n")
  
  formula_summary <- filtered_df %>%
    dplyr::filter(Diff_Formula != "") %>%
    dplyr::group_by(Diff_Formula) %>%
    dplyr::summarise(
      count = dplyr::n(),
      avg_CAR = mean(CAR, na.rm = TRUE),
      C = mean(C, na.rm = TRUE),
      H = mean(H, na.rm = TRUE),
      O = mean(O, na.rm = TRUE),
      N = mean(N, na.rm = TRUE),
      S = mean(S, na.rm = TRUE),
      P = mean(P, na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(count))
  
  # ========== 18. Write Summary Table ==========
  if (DBI::dbExistsTable(con, "formula_summary")) DBI::dbRemoveTable(con, "formula_summary")
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "formula_summary", formula_summary)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ formula_summary table written,", nrow(formula_summary), "unique molecular formulas in total\n")
  
  # ========== 19. Export All Summary Statistics ==========
  if (export_all_summary) {
    # Reorder columns for better readability
    summary_export <- formula_summary %>%
      dplyr::select(Diff_Formula, count, avg_CAR, C, H, O, N, S, P)
    
    write.csv(summary_export, file = all_summary_file, row.names = FALSE)
    if (verbose) cat("✅ All formula statistics exported:", all_summary_file, 
                     "(", nrow(formula_summary), "formulas with CHONSP values )\n")
  }
  
  # ========== 20. Export Top N ==========
  top_n_formulas <- head(formula_summary, top_n)
  
  if (verbose) {
    cat("\n🏆 Top", top_n, "most frequent molecular formulas (small molecules excluded):\n")
    # Display in console
    print(top_n_formulas %>% dplyr::select(Diff_Formula, count, avg_CAR, C, H, O, N, S, P))
  }
  
  if (export_top) {
    # Export with all element statistics
    top_export <- top_n_formulas %>%
      dplyr::select(Diff_Formula, count, avg_CAR, C, H, O, N, S, P)
    
    write.csv(top_export, file = top_output_file, row.names = FALSE)
    if (verbose) cat("\n✅ Top", top_n, "formulas exported:", top_output_file, 
                     "with CHONSP values\n")
  }
  
  # ========== 21. Prepare Return Object ==========
  stats <- list(
    total_rows = n,
    total_pairs = total_pairs,
    filtered_pairs = length(keep),
    bidirectional_records = if(track_direction) length(keep) * 2 else length(keep),
    small_molecules_count = nrow(small_df),
    filtered_count = nrow(filtered_df),
    unique_formulas = nrow(formula_summary),
    top_n = top_n,
    track_direction = track_direction
  )
  
  result <- list(
    db_connection = con,
    formula_summary = as.data.frame(formula_summary),
    top_n_formulas = as.data.frame(top_n_formulas),
    stats = stats
  )
  
  if (verbose) cat("\n🎉 Analysis completed!\n")
  
  return(invisible(result))
}