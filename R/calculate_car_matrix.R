#' Calculate CAR Matrix and Generate Molecular Formula Statistics
#'
#' @description
#' This function performs CAR (Carbon Aromaticity Ratio) matrix calculation,
#' filters molecular formulas, and generates frequency statistics.
#'
#' @param csv_file Character. Path to the input CSV file containing molecular formulas
#' @param db_path Character. Path to the SQLite database file (default: "FTICR_CAR_smart.db")
#' @param car_min Numeric. Minimum CAR threshold (default: 0.45)
#' @param car_max Numeric. Maximum CAR threshold (default: 1.0)
#' @param small_molecules Character vector. Molecules to be filtered out (default: predefined list)
#' @param export_top50 Logical. Whether to export top 50 formulas to CSV (default: TRUE)
#' @param output_file Character. Output file name for top 50 formulas (default: "Top50_formula.csv")
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item{db_connection: }{SQLite database connection object}
#'   \item{formula_summary: }{Data frame of formula frequency statistics}
#'   \item{top50: }{Top 50 most frequent formulas}
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
#' # Basic usage
#' result <- calculate_car_matrix("merged_molform_elements.csv")
#' 
#' # Custom parameters
#' result <- calculate_car_matrix(
#'   csv_file = "data.csv",
#'   db_path = "my_database.db",
#'   car_min = 0.5,
#'   car_max = 0.95,
#'   export_top50 = TRUE
#' )
#' 
#' # Access results
#' print(result$top50)
#' print(result$stats)
#' }

calculate_car_matrix <- function(
    csv_file,
    db_path = "FTICR_CAR_smart.db",
    car_min = 0.45,
    car_max = 1.0,
    small_molecules = c("C", "O", "CH2", "H2O", "CO", "O2", "CO2", "NH3", "NH"),
    export_top50 = TRUE,
    output_file = "Top50_formula.csv",
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
  
  # ========== 2. Read CSV Data ==========
  if (!file.exists(csv_file)) {
    stop(paste("File not found:", csv_file))
  }
  
  df <- read.csv(csv_file, stringsAsFactors = FALSE)[, c("MolForm", "C", "H", "N", "O", "S")]
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
  mat <- as.matrix(data[, c("C", "H", "N", "O", "S")])
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
  
  # ========== 8. Calculate Differences ==========
  A <- mat[pairs[1, keep], ];  B <- mat[pairs[2, keep], ]
  diff_mat <- abs(A - B)
  
  results <- data.frame(
    MolForm_1 = rownames(mat)[pairs[1, keep]],
    MolForm_2 = rownames(mat)[pairs[2, keep]],
    C = diff_mat[, "C"], H = diff_mat[, "H"], N = diff_mat[, "N"],
    O = diff_mat[, "O"], S = diff_mat[, "S"], CAR = CAR[keep]
  )
  
  if (verbose) cat("🧩 Difference matrix generated,", nrow(results), "records in total\n")
  
  # ========== 9. Write car_matrix Table ==========
  if (DBI::dbExistsTable(con, "car_matrix")) DBI::dbRemoveTable(con, "car_matrix")
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "car_matrix", results)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ car_matrix table written\n")
  
  # ========== 10. Molecular Formula Composition Function ==========
  compose_formula <- function(C, H, N, O, S) {
    vals <- as.numeric(c(C, H, N, O, S))
    names(vals) <- c("C", "H", "N", "O", "S")
    vals[is.na(vals)] <- 0; vals[vals < 0] <- 0
    vals <- round(vals)
    nonzero <- vals != 0
    if (!any(nonzero)) return("")
    paste0(names(vals)[nonzero],
           ifelse(vals[nonzero] == 1, "", vals[nonzero]),
           collapse = "")
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
    results$Diff_Formula <- mapply(function(C, H, N, O, S) {
      if (verbose) {
        counter <<- counter + 1
        if (counter %% tick_every == 0) pb2$tick(tick_every)
      }
      compose_formula(C, H, N, O, S)
    }, results$C, results$H, results$N, results$O, results$S)
  } else {
    data.table::setDT(results)
    el_names <- c("C", "H", "N", "O", "S")
    results[, Diff_Formula := apply(.SD, 1, function(x) {
      if (verbose) {
        counter <<- counter + 1
        if (counter %% tick_every == 0) pb2$tick(tick_every)
      }
      compose_formula(x[1], x[2], x[3], x[4], x[5])
    }), .SDcols = el_names]
  }
  
  if (verbose) {
    pb2$tick(row_n)
    cat("\n✅ Standardized molecular formulas generated\n")
  }
  
  # ========== 13. Write car_matrix_formula Table ==========
  if (DBI::dbExistsTable(con, "car_matrix_formula")) DBI::dbRemoveTable(con, "car_matrix_formula")
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "car_matrix_formula", results)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ car_matrix_formula table written (containing all molecular formulas)\n")
  
  # ========== 14. Filter Small Molecules ==========
  small_df <- subset(results, Diff_Formula %in% small_molecules)
  filtered_df <- subset(results, !(Diff_Formula %in% small_molecules))
  
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
  
  # ========== 17. Frequency Statistics ==========
  if (verbose) cat("📊 Performing Diff_Formula frequency statistics (small molecules excluded)...\n")
  
  formula_summary <- filtered_df %>%
    dplyr::filter(Diff_Formula != "") %>%
    dplyr::group_by(Diff_Formula) %>%
    dplyr::summarise(
      count = dplyr::n(),
      avg_CAR = mean(CAR, na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(count))
  
  # ========== 18. Write Summary Table ==========
  if (DBI::dbExistsTable(con, "formula_summary")) DBI::dbRemoveTable(con, "formula_summary")
  DBI::dbExecute(con, "BEGIN TRANSACTION;")
  DBI::dbWriteTable(con, "formula_summary", formula_summary)
  DBI::dbExecute(con, "COMMIT;")
  
  if (verbose) cat("✅ formula_summary table written,", nrow(formula_summary), "unique molecular formulas in total\n")
  
  # ========== 19. Export Top 50 ==========
  top50 <- head(formula_summary, 50)
  
  if (verbose) {
    cat("\n🏆 Top 50 most frequent molecular formulas (small molecules excluded):\n")
    print(top50)
  }
  
  if (export_top50) {
    write.csv(top50, file = output_file, row.names = FALSE)
    if (verbose) cat("\n✅ Exported:", output_file, "\n")
  }
  
  # ========== 20. Prepare Return Object ==========
  stats <- list(
    total_rows = n,
    total_pairs = total_pairs,
    filtered_pairs = length(keep),
    small_molecules_count = nrow(small_df),
    filtered_count = nrow(filtered_df),
    unique_formulas = nrow(formula_summary)
  )
  
  result <- list(
    db_connection = con,
    formula_summary = as.data.frame(formula_summary),
    top50 = as.data.frame(top50),
    stats = stats
  )
  
  if (verbose) cat("\n🎉 Analysis completed!\n")
  
  return(invisible(result))
}