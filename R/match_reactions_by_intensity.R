#' Match reactions between two molecule datasets based on intensity filtering
#'
#' This function filters precursor and product molecules based on intensity changes,
#' and matches them with a list of possible reaction deltas.
#'
#' @param file1 Path to the first molecular information CSV (e.g., inflow).
#' @param file2 Path to the second molecular information CSV (e.g., outflow).
#' @param reaction_delta_file Path to the reaction delta CSV file.
#' @param out_dir Directory to save output files.
#'
#' @return Two CSV files saved in `out_dir`: `network_edge.csv` and `reaction_summary.csv`.
#' @export
match_reactions_by_intensity <- function(file1, file2, reaction_delta_file, out_dir = ".") {
  library(dplyr)
  library(readr)

  message("📥 Step 1: Reading molecular information files...")

  mol1 <- tryCatch(read_csv(file1), error = function(e) stop("❌ Failed to read file1: ", e))
  mol2 <- tryCatch(read_csv(file2), error = function(e) stop("❌ Failed to read file2: ", e))

  names(mol1) <- trimws(names(mol1))
  names(mol2) <- trimws(names(mol2))

  if (!"MolForm" %in% names(mol1) || !"MolForm" %in% names(mol2)) {
    stop("❌ Error: 'MolForm' column is missing in input files.")
  }

  mol1 <- mol1 %>% rename(Formula = MolForm)
  mol2 <- mol2 %>% rename(Formula = MolForm)

  message("🧮 Step 2: Merging and filtering molecules by intensity ratio...")

  commom <- merge(mol1, mol2, by = "Formula", suffixes = c("_1", "_2")) %>%
    rename(abundance_1 = intensity_1, abundance_2 = intensity_2)

  data <- commom %>%
    select(Formula, abundance_1, C_1, H_1, O_1, N_1, S_1, Cl_1, Br_1, P_1, I_1,
           abundance_2, C_2, H_2, O_2, N_2, S_2, Cl_2, Br_2, P_2, I_2) %>%
    mutate(
      pre = ifelse(abundance_2 / abundance_1 < 0.5, 1, 0),
      pro = ifelse(abundance_2 / abundance_1 > 2, 2, 0)
    )

  message("📊 Step 3: Extracting precursor/product molecules...")

  data1 <- data %>% filter(pre == 1) %>%
    select(Formula, C = C_1, H = H_1, O = O_1, N = N_1, S = S_1, Cl = Cl_1, Br = Br_1, P = P_1, I = I_1)
  data2 <- data %>% filter(pro == 2) %>%
    select(Formula, C = C_2, H = H_2, O = O_2, N = N_2, S = S_2, Cl = Cl_2, Br = Br_2, P = P_2, I = I_2)

  unique1 <- mol1 %>% filter(!(Formula %in% commom$Formula)) %>%
    select(Formula, C, H, O, N, S, Cl, Br, P, I)
  unique2 <- mol2 %>% filter(!(Formula %in% commom$Formula)) %>%
    select(Formula, C, H, O, N, S, Cl, Br, P, I)

  mol1_filtered <- bind_rows(unique1, data1)
  mol2_filtered <- bind_rows(unique2, data2)

  message("🧪 Step 4: Reading reaction delta definitions...")

  reaction_delta <- tryCatch(read_csv(reaction_delta_file),
                             error = function(e) stop("❌ Failed to read reaction delta file: ", e))
  names(reaction_delta) <- trimws(names(reaction_delta))

  element_cols <- c("C", "H", "N", "O", "S", "Cl", "Br", "P", "I")

  mol1_filtered[element_cols] <- lapply(mol1_filtered[element_cols], as.numeric)
  mol2_filtered[element_cols] <- lapply(mol2_filtered[element_cols], as.numeric)
  reaction_delta[element_cols] <- lapply(reaction_delta[element_cols], as.numeric)

  message("🔍 Step 5: Matching reactions between molecules...")
  results <- data.frame()

  pb <- txtProgressBar(min = 0, max = nrow(reaction_delta), style = 3)

  for (i in 1:nrow(reaction_delta)) {
    reaction_name <- reaction_delta$reaction[i]
    delta <- as.numeric(reaction_delta[i, element_cols])

    for (j in 1:nrow(mol1_filtered)) {
      mol1_vals <- as.numeric(mol1_filtered[j, element_cols])
      transformed <- mol1_vals + delta

      match_idx <- which(apply(mol2_filtered[, element_cols], 1, function(x)
        all(abs(as.numeric(x) - transformed) < 1e-8)))

      if (length(match_idx) > 0) {
        for (k in match_idx) {
          results <- rbind(results, data.frame(
            Source = mol1_filtered$Formula[j],
            Target = mol2_filtered$Formula[k],
            Reaction = reaction_name,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  if (nrow(results) == 0) {
    warning("⚠️ No reactions matched. Check your input data and delta definitions.")
  }

  message("💾 Step 6: Saving outputs...")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_file1 <- file.path(out_dir, "network_edge.csv")
  out_file2 <- file.path(out_dir, "reaction_summary.csv")

  write_csv(results, out_file1)

  reaction_summary <- results %>% count(Reaction, name = "Count")
  write_csv(reaction_summary, out_file2)

  message("✅ Done! Results saved to:\n- ", out_file1, "\n- ", out_file2)
  message("📚 Reference: 10.1016/j.watres.2020.116484")
}