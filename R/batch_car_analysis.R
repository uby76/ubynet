#' Batch CAR Analysis for Multiple Samples
#'
#' @description
#' Processes multiple CSV files to generate CAR matrices in SQLite databases,
#' then queries and summarizes Top50 molecular formula occurrences across all samples.
#' Also exports Source-Target edge lists for each sample.
#'
#' @param sample_folder Character. Path to folder containing sample CSV files
#' @param db_output_folder Character. Path where SQLite databases will be saved
#' @param edge_output_folder Character. Path where edge list files will be saved
#' @param top50_file Character. Path to Top50_formula.csv file (default: "Top50_formula.csv")
#' @param car_min Numeric. Minimum CAR threshold (default: 0.45)
#' @param car_max Numeric. Maximum CAR threshold (default: 1.0)
#' @param small_molecules Character vector. Molecules to filter out
#' @param export_summary Logical. Whether to export summary tables (default: TRUE)
#' @param export_edges Logical. Whether to export edge lists (default: TRUE)
#' @param summary_wide_file Character. Output filename for wide summary table
#' @param summary_long_file Character. Output filename for long summary table
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item{db_paths: }{Named list of database file paths}
#'   \item{edge_paths: }{Named list of edge file paths}
#'   \item{wide_table: }{Data frame with samples as columns}
#'   \item{long_table: }{Data frame with one row per sample-formula combination}
#'   \item{stats: }{Summary statistics}
#' }
#'
#' @import DBI
#' @import RSQLite
#' @import data.table
#' @import dplyr
#' @import progress
#' @import tools
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- batch_car_analysis(
#'   sample_folder = "path/to/samples",
#'   db_output_folder = "path/to/databases",
#'   edge_output_folder = "path/to/edges"
#' )
#' 
#' # Access results
#' print(result$wide_table)
#' print(result$stats)
#' }

batch_car_analysis <- function(
        sample_folder,
        db_output_folder,
        edge_output_folder = "edge_lists",
        top50_file = "Top50_formula.csv",
        car_min = 0.45,
        car_max = 1.0,
        small_molecules = c("C", "O", "CH2", "H2O", "CO", "O2", "CO2", "NH3", "NH"),
        export_summary = TRUE,
        export_edges = TRUE,
        summary_wide_file = "Top50_All_Samples_Summary.csv",
        summary_long_file = "Top50_All_Samples_Long.csv",
        verbose = TRUE
) {
    
    # ========== 1. Check and Load Dependencies ==========
    required_packages <- c("DBI", "RSQLite", "progress", "data.table", "dplyr", "tools")
    
    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("Package", pkg, "is required but not installed. Please install it first."))
        }
    }
    
    if (verbose) cat("✅ Dependency check completed\n")
    
    # ========== 2. Validate Paths ==========
    if (!dir.exists(sample_folder)) {
        stop(paste("Sample folder not found:", sample_folder))
    }
    
    if (!dir.exists(db_output_folder)) {
        dir.create(db_output_folder, recursive = TRUE)
        if (verbose) cat("✅ Created database folder:", db_output_folder, "\n")
    }
    
    if (export_edges && !dir.exists(edge_output_folder)) {
        dir.create(edge_output_folder, recursive = TRUE)
        if (verbose) cat("✅ Created edge list folder:", edge_output_folder, "\n")
    }
    
    # ========== 3. Load Top50 ==========
    if (!file.exists(top50_file)) {
        stop(paste("Top50 file not found:", top50_file))
    }
    
    top50 <- read.csv(top50_file, stringsAsFactors = FALSE)
    target_formulas <- top50$Diff_Formula
    
    if (verbose) cat("✅ Loaded Top50,", length(target_formulas), "formulas in total\n")
    
    # ========== 4. Get All CSV Files ==========
    csv_files <- list.files(sample_folder, pattern = "\\.csv$", full.names = TRUE)
    
    if (length(csv_files) == 0) {
        stop("No CSV files found in sample folder")
    }
    
    if (verbose) cat("📂 Found", length(csv_files), "sample files\n\n")
    
    # ========== 5. Helper Function: Process Single Sample ==========
    process_sample_to_db <- function(csv_path) {
        sample_name <- tools::file_path_sans_ext(basename(csv_path))
        db_path <- file.path(db_output_folder, paste0(sample_name, ".db"))
        
        if (verbose) {
            cat(paste0(rep("=", 62), collapse = ""), "\n")
            cat("🔵 Processing sample:", sample_name, "\n")
            cat(paste0(rep("=", 62), collapse = ""), "\n")
        }
        
        tryCatch({
            df <- read.csv(csv_path, stringsAsFactors = FALSE)
            
            required_cols <- c("C", "H", "N", "O", "S")
            molform_col <- intersect(c("MolForm", "sumFormula", "Formula", "Molecular_Formula"), names(df))
            
            if (length(molform_col) == 0) {
                if (verbose) cat("⚠️ Skipped: Molecular formula column not found\n\n")
                return(NULL)
            }
            
            if (!all(required_cols %in% names(df))) {
                if (verbose) cat("⚠️ Skipped: Missing required columns (C,H,N,O,S)\n\n")
                return(NULL)
            }
            
            df <- df[, c(molform_col[1], required_cols)]
            names(df)[1] <- "MolForm"
            num_cols <- sapply(df, is.numeric)
            df[num_cols][is.na(df[num_cols])] <- 0
            
            if (verbose) cat("✅ Data loaded,", nrow(df), "rows in total\n")
            
            con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
            
            if (DBI::dbExistsTable(con, "origin")) DBI::dbRemoveTable(con, "origin")
            DBI::dbWriteTable(con, "origin", df)
            if (verbose) cat("✅ Origin table written\n")
            
            mat <- as.matrix(df[, c("C", "H", "N", "O", "S")])
            rownames(mat) <- df$MolForm
            n <- nrow(mat)
            
            pairs <- combn(n, 2)
            total_pairs <- ncol(pairs)
            if (verbose) cat("⚙️ Generated", total_pairs, "combinations\n")
            
            if (verbose) {
                pb <- progress::progress_bar$new(
                    total = total_pairs,
                    format = "⏳ [:bar] :percent Remaining: :eta",
                    clear = FALSE, width = 70
                )
            }
            
            C1 <- mat[pairs[1, ], "C"]; H1 <- mat[pairs[1, ], "H"]; O1 <- mat[pairs[1, ], "O"]
            C2 <- mat[pairs[2, ], "C"]; H2 <- mat[pairs[2, ], "H"]; O2 <- mat[pairs[2, ], "O"]
            
            CAR <- (pmin(C1, C2) / pmax(C1, C2)) *
                (pmin(H1, H2) / pmax(H1, H2)) *
                (pmin(O1, O2) / pmax(O1, O2))
            CAR[is.na(CAR)] <- 0
            keep <- which(CAR >= car_min & CAR <= car_max)
            
            if (verbose) {
                pb$tick(total_pairs)
                cat("\n✅ Filtered", length(keep), "pairs meeting CAR ∈ [", car_min, ",", car_max, "]\n")
            }
            
            A <- mat[pairs[1, keep], ];  B <- mat[pairs[2, keep], ]
            diff_mat <- abs(A - B)
            
            results <- data.frame(
                MolForm_1 = rownames(mat)[pairs[1, keep]],
                MolForm_2 = rownames(mat)[pairs[2, keep]],
                C = diff_mat[, "C"], H = diff_mat[, "H"], N = diff_mat[, "N"],
                O = diff_mat[, "O"], S = diff_mat[, "S"], CAR = CAR[keep]
            )
            
            if (verbose) cat("🧩 Difference matrix generated,", nrow(results), "records in total\n")
            
            if (DBI::dbExistsTable(con, "car_matrix")) DBI::dbRemoveTable(con, "car_matrix")
            DBI::dbExecute(con, "BEGIN TRANSACTION;")
            DBI::dbWriteTable(con, "car_matrix", results)
            DBI::dbExecute(con, "COMMIT;")
            if (verbose) cat("✅ car_matrix table written\n")
            
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
            
            row_n <- nrow(results)
            mode <- if (row_n <= 50000) "mapply" else "data.table"
            if (verbose) cat("🧠 Auto-selected mode:", mode, "(Data rows:", row_n, ")\n")
            
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
            
            if (DBI::dbExistsTable(con, "car_matrix_formula")) DBI::dbRemoveTable(con, "car_matrix_formula")
            DBI::dbExecute(con, "BEGIN TRANSACTION;")
            DBI::dbWriteTable(con, "car_matrix_formula", results)
            DBI::dbExecute(con, "COMMIT;")
            if (verbose) cat("✅ car_matrix_formula table written\n")
            
            small_df <- subset(results, Diff_Formula %in% small_molecules)
            filtered_df <- subset(results, !(Diff_Formula %in% small_molecules))
            if (verbose) cat("🚫 Small molecule filtering:", nrow(small_df), "small molecules,", nrow(filtered_df), "remaining\n")
            
            if (DBI::dbExistsTable(con, "car_matrix_smallmolecules")) DBI::dbRemoveTable(con, "car_matrix_smallmolecules")
            DBI::dbExecute(con, "BEGIN TRANSACTION;")
            DBI::dbWriteTable(con, "car_matrix_smallmolecules", small_df)
            DBI::dbExecute(con, "COMMIT;")
            
            if (DBI::dbExistsTable(con, "car_matrix_filtered")) DBI::dbRemoveTable(con, "car_matrix_filtered")
            DBI::dbExecute(con, "BEGIN TRANSACTION;")
            DBI::dbWriteTable(con, "car_matrix_filtered", filtered_df)
            DBI::dbExecute(con, "COMMIT;")
            
            DBI::dbDisconnect(con)
            if (verbose) cat("✅ Database generation completed:", db_path, "\n\n")
            
            return(db_path)
            
        }, error = function(e) {
            if (verbose) cat("❌ Error:", e$message, "\n\n")
            return(NULL)
        })
    }
    
    # ========== 6. Helper Function: Export Edge List for Top50 Formulas ==========
    export_edge_list <- function(db_path, target_formulas, edge_folder) {
        sample_name <- tools::file_path_sans_ext(basename(db_path))
        edge_file <- file.path(edge_folder, paste0(sample_name, "_edges.csv"))
        
        tryCatch({
            con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
            
            if (!DBI::dbExistsTable(con, "car_matrix_formula")) {
                DBI::dbDisconnect(con)
                if (verbose) cat("⚠️ No car_matrix_formula table found in", sample_name, "\n")
                return(NULL)
            }
            
            # Query all rows where Diff_Formula matches Top50
            query <- sprintf(
                "SELECT MolForm_1 AS Source, MolForm_2 AS Target, Diff_Formula, CAR 
         FROM car_matrix_formula 
         WHERE Diff_Formula IN (%s)",
                paste0("'", target_formulas, "'", collapse = ", ")
            )
            
            edge_data <- DBI::dbGetQuery(con, query)
            DBI::dbDisconnect(con)
            
            if (nrow(edge_data) == 0) {
                if (verbose) cat("⚠️ No matching edges found for", sample_name, "\n")
                return(NULL)
            }
            
            # Export edge list
            write.csv(edge_data, edge_file, row.names = FALSE)
            
            if (verbose) cat("📝 Exported", nrow(edge_data), "edges to:", edge_file, "\n")
            
            return(edge_file)
            
        }, error = function(e) {
            if (verbose) cat("❌ Edge export error (", sample_name, "):", e$message, "\n")
            return(NULL)
        })
    }
    
    # ========== 7. Helper Function: Query Top50 from Database ==========
    query_top50_from_db <- function(db_path, target_formulas) {
        sample_name <- tools::file_path_sans_ext(basename(db_path))
        
        tryCatch({
            con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
            
            if (!DBI::dbExistsTable(con, "car_matrix_formula")) {
                DBI::dbDisconnect(con)
                return(NULL)
            }
            
            query <- sprintf("SELECT Diff_Formula, CAR FROM car_matrix_formula WHERE Diff_Formula IN (%s)",
                             paste0("'", target_formulas, "'", collapse = ", "))
            matched_data <- DBI::dbGetQuery(con, query)
            
            DBI::dbDisconnect(con)
            
            if (nrow(matched_data) == 0) {
                return(data.frame(
                    Diff_Formula = target_formulas,
                    count = 0,
                    avg_CAR = 0
                ))
            }
            
            freq_result <- matched_data %>%
                dplyr::group_by(Diff_Formula) %>%
                dplyr::summarise(
                    count = dplyr::n(),
                    avg_CAR = mean(CAR, na.rm = TRUE),
                    .groups = "drop"
                )
            
            result <- data.frame(Diff_Formula = target_formulas) %>%
                dplyr::left_join(freq_result, by = "Diff_Formula") %>%
                dplyr::mutate(
                    count = ifelse(is.na(count), 0, count),
                    avg_CAR = ifelse(is.na(avg_CAR), 0, avg_CAR)
                )
            
            return(result)
            
        }, error = function(e) {
            if (verbose) cat("❌ Query error (", sample_name, "):", e$message, "\n")
            return(NULL)
        })
    }
    
    # ========== 8. Phase 1: Generate All Databases ==========
    if (verbose) {
        cat("\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("🚀 Phase 1: Generate Databases\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n\n")
    }
    
    db_paths <- list()
    
    for (csv_file in csv_files) {
        db_path <- process_sample_to_db(csv_file)
        if (!is.null(db_path)) {
            sample_name <- tools::file_path_sans_ext(basename(db_path))
            db_paths[[sample_name]] <- db_path
        }
    }
    
    if (length(db_paths) == 0) {
        stop("No databases were successfully generated")
    }
    
    if (verbose) cat("✅ Successfully generated", length(db_paths), "databases\n\n")
    
    # ========== 9. Phase 2: Export Edge Lists (NEW) ==========
    if (export_edges) {
        if (verbose) {
            cat(paste0(rep("=", 62), collapse = ""), "\n")
            cat("📊 Phase 2: Export Edge Lists\n")
            cat(paste0(rep("=", 62), collapse = ""), "\n\n")
        }
        
        edge_paths <- list()
        
        for (sample_name in names(db_paths)) {
            if (verbose) cat("📝 Exporting edges for:", sample_name, "\n")
            edge_file <- export_edge_list(db_paths[[sample_name]], target_formulas, edge_output_folder)
            
            if (!is.null(edge_file)) {
                edge_paths[[sample_name]] <- edge_file
            }
        }
        
        if (verbose) cat("\n✅ Successfully exported", length(edge_paths), "edge files\n\n")
    } else {
        edge_paths <- list()
    }
    
    # ========== 10. Phase 3: Query All Databases for Top50 ==========
    if (verbose) {
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("🔍 Phase 3: Query Top50 Matches\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n\n")
    }
    
    all_results <- list()
    
    for (sample_name in names(db_paths)) {
        if (verbose) cat("🔍 Querying sample:", sample_name, "\n")
        result <- query_top50_from_db(db_paths[[sample_name]], target_formulas)
        
        if (!is.null(result)) {
            result$Sample <- sample_name
            all_results[[sample_name]] <- result
            if (verbose) cat("   ✅ Detected", sum(result$count > 0), "formulas\n\n")
        }
    }
    
    # ========== 11. Phase 4: Merge Results and Generate Output ==========
    if (verbose) {
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("📊 Phase 4: Generate Summary Table\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n\n")
    }
    
    wide_table <- data.frame(Diff_Formula = target_formulas)
    
    for (sample_name in names(all_results)) {
        sample_data <- all_results[[sample_name]]
        wide_table <- wide_table %>%
            dplyr::left_join(
                sample_data %>% dplyr::select(Diff_Formula, count),
                by = "Diff_Formula"
            )
        names(wide_table)[ncol(wide_table)] <- sample_name
    }
    
    wide_table[is.na(wide_table)] <- 0
    
    # Handle single or multiple samples
    if (ncol(wide_table) == 2) {
        # Only one sample
        wide_table$Total_Count <- wide_table[, 2]
        wide_table$Detection_Rate <- ifelse(wide_table[, 2] > 0, 100, 0)
    } else {
        # Multiple samples
        wide_table$Total_Count <- rowSums(wide_table[, -1, drop = FALSE])
        wide_table$Detection_Rate <- round(
            rowSums(wide_table[, -c(1, ncol(wide_table)), drop = FALSE] > 0) / (ncol(wide_table) - 2) * 100, 2
        )
    }
    
    wide_table <- wide_table %>% dplyr::arrange(dplyr::desc(Total_Count))
    
    if (verbose) {
        cat("🏆 Top50 formula occurrence across all samples:\n")
        print(as.data.frame(head(wide_table, 20)), row.names = FALSE)
    }
    
    long_table <- do.call(rbind, all_results)
    
    if (export_summary) {
        write.csv(wide_table, summary_wide_file, row.names = FALSE)
        if (verbose) cat("\n✅ Summary table saved to:", summary_wide_file, "\n")
        
        write.csv(long_table, summary_long_file, row.names = FALSE)
        if (verbose) cat("✅ Detailed table saved to:", summary_long_file, "\n")
    }
    
    # ========== 12. Statistics Summary ==========
    stats <- list(
        databases_generated = length(db_paths),
        edge_files_exported = length(edge_paths),
        samples_processed = length(all_results),
        top50_formulas = nrow(wide_table),
        detected_at_least_once = sum(wide_table$Total_Count > 0),
        not_detected = sum(wide_table$Total_Count == 0),
        avg_detection_rate = round(mean(wide_table$Detection_Rate), 2)
    )
    
    if (verbose) {
        cat("\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("📈 Statistics Summary\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("   🗄️ Databases generated:", stats$databases_generated, "\n")
        cat("   📝 Edge files exported:", stats$edge_files_exported, "\n")
        cat("   📦 Samples processed:", stats$samples_processed, "\n")
        cat("   🎯 Top50 formulas:", stats$top50_formulas, "\n")
        cat("   ✅ Detected at least once:", stats$detected_at_least_once, "\n")
        cat("   ❌ Not detected in any sample:", stats$not_detected, "\n")
        cat("   📊 Average detection rate:", stats$avg_detection_rate, "%\n")
        
        cat("\n🥇 Top 10 most frequent formulas:\n")
        top10 <- head(wide_table, 10)
        print(as.data.frame(top10), row.names = FALSE)
        
        cat("\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("🎯 All tasks completed!\n")
        cat(paste0(rep("=", 62), collapse = ""), "\n")
        cat("\n📂 Generated files:\n")
        cat("   📁 Database folder:", db_output_folder, "\n")
        if (export_edges) {
            cat("   📁 Edge list folder:", edge_output_folder, "\n")
        }
        cat("   📄", summary_wide_file, "(Summary table)\n")
        cat("   📄", summary_long_file, "(Detailed table)\n")
    }
    
    # ========== 13. Return Results ==========
    result <- list(
        db_paths = db_paths,
        edge_paths = edge_paths,
        wide_table = as.data.frame(wide_table),
        long_table = as.data.frame(long_table),
        stats = stats
    )
    
    return(invisible(result))
}