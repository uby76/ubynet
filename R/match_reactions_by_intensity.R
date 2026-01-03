# --- Utility Function: Standardize Columns (Elements + Required) ---
#' Utility function to ensure required columns exist in a dataframe and
#' to fill missing element columns with 0.
#' @param df Input data frame.
#' @param required_cols Character vector of columns that should be kept (e.g., c("Formula","intensity","C","H","O")).
#' @param element_cols_all Character vector of supported element columns.
#' @keywords internal
# ============================================================
# 0. Global settings
# ============================================================
element_cols_all <- c("C","H","O","N","S","Cl","Br","P","I")

# ============================================================
# 1. Utility: Standardize columns
#    - No element stop
#    - others 0
# ============================================================
standardize_elements <- function(df, required_cols) {

  names(df) <- base::trimws(names(df))
  required_cols <- base::unique(required_cols)

  missing <- base::setdiff(required_cols, names(df))
  if (length(missing) > 0) {

    missing_elements <- base::intersect(missing, element_cols_all)
    missing_non_elements <- base::setdiff(missing, element_cols_all)

    if (length(missing_non_elements) > 0) {
      base::stop(
        "Error: Missing required non-element columns: ",
        base::paste(missing_non_elements, collapse = ", ")
      )
    }

    if (length(missing_elements) > 0) {
      for (el in missing_elements) df[[el]] <- 0L
      base::message(
        "    [Info] Added missing element columns: ",
        base::paste(missing_elements, collapse = ", ")
      )
    }
  }

  df <- dplyr::select(df, dplyr::all_of(required_cols))
  df
}

# ============================================================
# 2. Utility: Parse MolForm
# ============================================================
parse_formula_elements <- function(molforms, element_cols) {

  molforms <- base::gsub("\\s+", "", base::as.character(molforms))
  counts <- base::matrix(0L, nrow = length(molforms), ncol = length(element_cols))
  base::colnames(counts) <- element_cols
  unknown <- character()

  for (i in seq_along(molforms)) {
    f <- molforms[i]
    if (base::is.na(f) || f == "") next

    tokens <- base::regmatches(
      f,
      base::gregexpr("[A-Z][a-z]?[0-9]*", f, perl = TRUE)
    )[[1]]

    if (length(tokens) == 0) next

    for (tk in tokens) {
      el <- base::sub("([A-Z][a-z]?)[0-9]*", "\\1", tk, perl = TRUE)
      n  <- base::sub("^[A-Za-z]+", "", tk, perl = TRUE)
      n  <- if (n == "") 1L else base::suppressWarnings(base::as.integer(n))

      if (base::is.na(n)) next

      if (el %in% element_cols) {
        counts[i, el] <- counts[i, el] + n
      } else {
        unknown <- c(unknown, el)
      }
    }
  }

  list(
    counts = base::as.data.frame(counts),
    unknown_elements = base::unique(unknown)
  )
}

# ============================================================
# 3. Utility: Ensure elements from MolForm + mismatch export
# ============================================================
ensure_elements_from_molform <- function(df) {

  if (!"MolForm" %in% names(df)) {
    base::stop("Error: 'MolForm' column is required.")
  }

  molforms <- base::trimws(base::as.character(df$MolForm))
  if (base::any(base::is.na(molforms) | molforms == "")) {
    base::stop("Error: MolForm column contains missing or empty values.")
  }

  parsed  <- parse_formula_elements(molforms, element_cols_all)
  derived <- parsed$counts

  if (length(parsed$unknown_elements) > 0) {
    base::warning(
      "Unsupported elements in MolForm ignored: ",
      base::paste(parsed$unknown_elements, collapse = ", ")
    )
  }

  # 确保元素列存在
  for (el in element_cols_all) {
    if (!el %in% names(df)) df[[el]] <- 0L
    df[[el]] <- base::suppressWarnings(base::as.integer(df[[el]]))
  }

  mismatch_flag <- rep(FALSE, nrow(df))
  mismatch_detail <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    diff_el <- character()
    for (el in element_cols_all) {
      if (!is.na(df[[el]][i]) && df[[el]][i] != derived[[el]][i]) {
        diff_el <- c(diff_el, el)
      }
    }
    if (length(diff_el) > 0) {
      mismatch_flag[i] <- TRUE
      mismatch_detail[[i]] <- diff_el
    }
  }

  mismatch_rows <- NULL
  if (any(mismatch_flag)) {

    raw_part <- df[mismatch_flag, c("MolForm", element_cols_all), drop = FALSE]
    colnames(raw_part) <- c("MolForm", paste0(element_cols_all, "_raw"))

    form_part <- derived[mismatch_flag, , drop = FALSE]
    colnames(form_part) <- paste0(element_cols_all, "_molform")

    mismatch_rows <- base::cbind(
      raw_part,
      form_part,
      mismatch_elements = sapply(
        mismatch_detail[mismatch_flag],
        function(x) base::paste(x, collapse = ";")
      )
    )

    base::warning(
      "Element columns differ from MolForm and were replaced: ",
      base::paste(base::unique(unlist(mismatch_detail)), collapse = ", ")
    )
  }

  # 强制以 MolForm 推导值为准
  for (el in element_cols_all) df[[el]] <- derived[[el]]

  list(
    df_clean = df,
    mismatch_rows = mismatch_rows
  )
}

# ============================================================
# 4. Main function
# ============================================================
match_reactions_by_intensity <- function(file1, file2, reaction_delta_file,
                                         out_dir = ".", use_memory_db = TRUE) {

  if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)

  # -------------------------
  # Step 1. Read data
  # -------------------------
  base::message("Step 1: Reading data...")
  mol1 <- readr::read_csv(file1, show_col_types = FALSE)
  mol2 <- readr::read_csv(file2, show_col_types = FALSE)

  if (!"MolForm" %in% names(mol1) || !"MolForm" %in% names(mol2)) {
    base::stop("Error: 'MolForm' column is missing in input files.")
  }
  if (!"intensity" %in% names(mol1) || !"intensity" %in% names(mol2)) {
    base::stop("Error: 'intensity' column is missing in input files.")
  }

  mol1 <- standardize_elements(mol1, c("MolForm","intensity",element_cols_all))
  mol2 <- standardize_elements(mol2, c("MolForm","intensity",element_cols_all))

  mol1$intensity <- base::suppressWarnings(base::as.numeric(mol1$intensity))
  mol2$intensity <- base::suppressWarnings(base::as.numeric(mol2$intensity))

  if (any(!is.finite(mol1$intensity)) || any(!is.finite(mol2$intensity))) {
    base::stop("Error: intensity contains NA/Inf values.")
  }

  res1 <- ensure_elements_from_molform(mol1)
  res2 <- ensure_elements_from_molform(mol2)

  mol1 <- res1$df_clean
  mol2 <- res2$df_clean

  if (!is.null(res1$mismatch_rows))
    readr::write_csv(
      res1$mismatch_rows,
      base::file.path(out_dir, "mol1_molform_element_mismatch.csv")
    )

  if (!is.null(res2$mismatch_rows))
    readr::write_csv(
      res2$mismatch_rows,
      base::file.path(out_dir, "mol2_molform_element_mismatch.csv")
    )

  # -------------------------
  # Step 2. Merge & FC filter
  # -------------------------
  base::message("Step 2: Filtering by intensity change...")

  common <- base::merge(
    mol1, mol2,
    by = "MolForm",
    suffixes = c("_1","_2")
  )

  ratio <- base::ifelse(
    common$intensity_1 > 0,
    common$intensity_2 / common$intensity_1,
    NA_real_
  )

  pre_flag <- is.finite(ratio) & ratio < 0.5
  pro_flag <- is.finite(ratio) & ratio > 2

  precursors <- common[pre_flag, c("MolForm", paste0(element_cols_all,"_1"))]
  products   <- common[pro_flag, c("MolForm", paste0(element_cols_all,"_2"))]

  colnames(precursors) <- c("MolForm", element_cols_all)
  colnames(products)   <- c("MolForm", element_cols_all)

  mol1_unique <- mol1[!(mol1$MolForm %in% common$MolForm),
                      c("MolForm", element_cols_all)]

  mol2_unique <- mol2[!(mol2$MolForm %in% common$MolForm),
                      c("MolForm", element_cols_all)]

  mol1f <- base::rbind(mol1_unique, precursors)
  mol2f <- base::rbind(mol2_unique, products)

  mol1f <- ensure_elements_from_molform(mol1f)$df_clean
  mol2f <- ensure_elements_from_molform(mol2f)$df_clean

  # -------------------------
  # Step 3. Reaction delta
  # -------------------------
  base::message("Step 3: Reading reaction delta...")

  reaction_delta <- readr::read_csv(reaction_delta_file, show_col_types = FALSE)
  if (!"reaction" %in% names(reaction_delta)) {
    base::stop("Error: 'reaction' column missing in reaction delta file.")
  }

  reaction_delta <- standardize_elements(
    reaction_delta,
    c("reaction", element_cols_all)
  )

  for (el in element_cols_all) {
    reaction_delta[[el]] <- base::replace(
      base::suppressWarnings(base::as.integer(reaction_delta[[el]])),
      is.na(reaction_delta[[el]]),
      0L
    )
  }

  # -------------------------
  # Step 4. SQLite matching
  # -------------------------
  base::message("Step 4: Matching reactions...")

  con <- if (use_memory_db)
    DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  else
    DBI::dbConnect(RSQLite::SQLite(), base::file.path(out_dir, "temp.db"))

  base::on.exit(DBI::dbDisconnect(con), add = TRUE)

  DBI::dbWriteTable(con,"mol1",mol1f,overwrite=TRUE)
  DBI::dbWriteTable(con,"mol2",mol2f,overwrite=TRUE)
  DBI::dbWriteTable(con,"rxn",reaction_delta,overwrite=TRUE)

  DBI::dbExecute(con,"CREATE INDEX idx_mol1 ON mol1(C,H,O,N,S,Cl,Br,P,I)")
  DBI::dbExecute(con,"CREATE INDEX idx_mol2 ON mol2(C,H,O,N,S,Cl,Br,P,I)")

  results_list <- list()
  pb <- utils::txtProgressBar(0,nrow(reaction_delta),style=3)

  for (i in seq_len(nrow(reaction_delta))) {

    d <- reaction_delta[i,]

    q <- base::sprintf(
      "SELECT m1.MolForm AS Source, m2.MolForm AS Target, '%s' AS Reaction
       FROM mol1 m1 JOIN mol2 m2 ON
       m2.C  = m1.C  + %d AND m2.H  = m1.H  + %d AND
       m2.O  = m1.O  + %d AND m2.N  = m1.N  + %d AND
       m2.S  = m1.S  + %d AND m2.Cl = m1.Cl + %d AND
       m2.Br = m1.Br + %d AND m2.P  = m1.P  + %d AND
       m2.I  = m1.I  + %d",
      d$reaction,
      d$C,d$H,d$O,d$N,d$S,d$Cl,d$Br,d$P,d$I
    )

    res <- DBI::dbGetQuery(con,q)
    if (nrow(res) > 0) results_list[[length(results_list) + 1]] <- res
    utils::setTxtProgressBar(pb,i)
  }
  close(pb)

  results <- if (length(results_list) > 0)
    dplyr::bind_rows(results_list)
  else
    data.frame(Source=character(),Target=character(),Reaction=character())

  # -------------------------
  # Step 5. Output
  # -------------------------
  readr::write_csv(results, base::file.path(out_dir,"network_edge.csv"))
  readr::write_csv(
    dplyr::count(results,Reaction,name="Count"),
    base::file.path(out_dir,"reaction_summary.csv")
  )

  base::message("Done.")
  invisible(results)
}
