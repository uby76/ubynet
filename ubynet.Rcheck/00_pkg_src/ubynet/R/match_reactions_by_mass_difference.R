#' Match reactions between two molecule datasets based on mass difference
#'
#' This function filters precursor and product molecules based on intensity changes,
#' and matches them with a list of possible reaction mass differences.
#'
#' @param file1 Path to the first molecular information CSV (e.g., inflow).
#' @param file2 Path to the second molecular information CSV (e.g., outflow).
#' @param reaction_delta_file Path to the reaction delta CSV file (with mass differences).
#' @param out_dir Directory to save output files.
#' @param mass_tolerance Mass tolerance for matching (default: 0.005 Da).
#'
#' @return Two CSV files saved in `out_dir`: `network_edge.csv` and `reaction_summary.csv`.
#' @export
match_reactions_by_mass_difference <- function(file1, file2, reaction_delta_file, out_dir = ".", mass_tolerance = 0.005) {
  library(dplyr)
  library(readr)

  message("📥 Step 1: Reading molecular information files...")

  mol1 <- tryCatch(read_csv(file1), error = function(e) stop("❌ Failed to read file1: ", e))
  mol2 <- tryCatch(read_csv(file2), error = function(e) stop("❌ Failed to read file2: ", e))

  names(mol1) <- trimws(names(mol1))
  names(mol2) <- trimws(names(mol2))

  # 检查必要的列
  required_cols_mol <- c("MolForm", "Mass")
  if (!all(required_cols_mol %in% names(mol1)) || !all(required_cols_mol %in% names(mol2))) {
    stop("❌ Error: 'MolForm' and 'Mass' columns are required in input files.")
  }

  mol1 <- mol1 %>% rename(Formula = MolForm)
  mol2 <- mol2 %>% rename(Formula = MolForm)

  message("🧮 Step 2: Merging and filtering molecules by intensity ratio...")

  common <- merge(mol1, mol2, by = "Formula", suffixes = c("_1", "_2")) %>%
    rename(abundance_1 = intensity_1, abundance_2 = intensity_2)

  data <- common %>%
    select(Formula, abundance_1, Mass_1, abundance_2, Mass_2) %>%
    mutate(
      pre = ifelse(abundance_2 / abundance_1 < 0.5, 1, 0),
      pro = ifelse(abundance_2 / abundance_1 > 2, 2, 0)
    )

  message("📊 Step 3: Extracting precursor/product molecules...")

  # 提取前体分子 (强度降低的分子)
  data1 <- data %>% filter(pre == 1) %>%
    select(Formula, Mass = Mass_1)
  
  # 提取产物分子 (强度增加的分子)
  data2 <- data %>% filter(pro == 2) %>%
    select(Formula, Mass = Mass_2)

  # 添加只在一个文件中出现的分子
  unique1 <- mol1 %>% filter(!(Formula %in% common$Formula)) %>%
    select(Formula, Mass)
  unique2 <- mol2 %>% filter(!(Formula %in% common$Formula)) %>%
    select(Formula, Mass)

  mol1_filtered <- bind_rows(unique1, data1)
  mol2_filtered <- bind_rows(unique2, data2)

  message("🧪 Step 4: Reading reaction mass difference definitions...")

  reaction_delta <- tryCatch(read_csv(reaction_delta_file),
                             error = function(e) stop("❌ Failed to read reaction delta file: ", e))
  names(reaction_delta) <- trimws(names(reaction_delta))

  # 检查反应差异文件必要的列
  required_cols_delta <- c("reaction", "mass_difference")
  if (!all(required_cols_delta %in% names(reaction_delta))) {
    stop("❌ Error: 'reaction' and 'mass_difference' columns are required in reaction delta file.")
  }

  # 确保质量相关列是数值型
  mol1_filtered$Mass <- as.numeric(mol1_filtered$Mass)
  mol2_filtered$Mass <- as.numeric(mol2_filtered$Mass)
  reaction_delta$mass_difference <- as.numeric(reaction_delta$mass_difference)

  message("🔍 Step 5: Matching reactions based on mass differences...")
  results <- data.frame()

  pb <- txtProgressBar(min = 0, max = nrow(reaction_delta), style = 3)

  for (i in 1:nrow(reaction_delta)) {
    reaction_name <- reaction_delta$reaction[i]
    mass_diff <- reaction_delta$mass_difference[i]

    for (j in 1:nrow(mol1_filtered)) {
      precursor_mass <- mol1_filtered$Mass[j]
      expected_product_mass <- precursor_mass + mass_diff

      # 在产物分子中寻找匹配的质量
      mass_matches <- which(abs(mol2_filtered$Mass - expected_product_mass) <= mass_tolerance)

      if (length(mass_matches) > 0) {
        for (k in mass_matches) {
          actual_mass_diff <- mol2_filtered$Mass[k] - precursor_mass
          results <- rbind(results, data.frame(
            Source = mol1_filtered$Formula[j],
            Target = mol2_filtered$Formula[k],
            Reaction = reaction_name,
            Precursor_Mass = precursor_mass,
            Product_Mass = mol2_filtered$Mass[k],
            Expected_Mass_Diff = mass_diff,
            Actual_Mass_Diff = actual_mass_diff,
            Mass_Error = abs(actual_mass_diff - mass_diff),
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  if (nrow(results) == 0) {
    warning("⚠️ No reactions matched. Check your input data, mass differences, and tolerance settings.")
  }

  message("💾 Step 6: Saving outputs...")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_file1 <- file.path(out_dir, "network_edge.csv")
  out_file2 <- file.path(out_dir, "reaction_summary.csv")

  write_csv(results, out_file1)

  reaction_summary <- results %>% 
    group_by(Reaction) %>%
    summarise(
      Count = n(),
      Avg_Mass_Error = mean(Mass_Error),
      Max_Mass_Error = max(Mass_Error),
      .groups = 'drop'
    )
  write_csv(reaction_summary, out_file2)

  message("✅ Done! Results saved to:\n- ", out_file1, "\n- ", out_file2)
  message("📊 Summary:")
  message("- Total reactions found: ", nrow(results))
  message("- Mass tolerance used: ±", mass_tolerance, " Da")
  message("📚 Reference: 10.1016/j.watres.2020.116484")
  
  return(list(
    network_edges = results,
    reaction_summary = reaction_summary
  ))
}