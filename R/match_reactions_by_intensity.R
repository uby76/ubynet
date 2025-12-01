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
match_reactions_by_intensity <- function(file1, file2, reaction_delta_file, 
                                        out_dir = ".", use_memory_db = TRUE) {
  
  # 加载必要的包
  if (!require("dplyr")) install.packages("dplyr")
  if (!require("readr")) install.packages("readr")
  if (!require("RSQLite")) install.packages("RSQLite")
  if (!require("DBI")) install.packages("DBI")
  
  library(dplyr)
  library(readr)
  library(RSQLite)
  library(DBI)
  
  message("📥 Step 1: Reading molecular information files...")
  
  mol1 <- tryCatch(read_csv(file1, show_col_types = FALSE), 
                   error = function(e) stop("❌ Failed to read file1: ", e))
  mol2 <- tryCatch(read_csv(file2, show_col_types = FALSE), 
                   error = function(e) stop("❌ Failed to read file2: ", e))
  
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
    select(Formula, C = C_1, H = H_1, O = O_1, N = N_1, S = S_1, 
           Cl = Cl_1, Br = Br_1, P = P_1, I = I_1)
  data2 <- data %>% filter(pro == 2) %>%
    select(Formula, C = C_2, H = H_2, O = O_2, N = N_2, S = S_2, 
           Cl = Cl_2, Br = Br_2, P = P_2, I = I_2)
  
  unique1 <- mol1 %>% filter(!(Formula %in% commom$Formula)) %>%
    select(Formula, C, H, O, N, S, Cl, Br, P, I)
  unique2 <- mol2 %>% filter(!(Formula %in% commom$Formula)) %>%
    select(Formula, C, H, O, N, S, Cl, Br, P, I)
  
  mol1_filtered <- bind_rows(unique1, data1)
  mol2_filtered <- bind_rows(unique2, data2)
  
  message("🧪 Step 4: Reading reaction delta definitions...")
  
  reaction_delta <- tryCatch(read_csv(reaction_delta_file, show_col_types = FALSE),
                             error = function(e) stop("❌ Failed to read reaction delta file: ", e))
  names(reaction_delta) <- trimws(names(reaction_delta))
  
  element_cols <- c("C", "H", "N", "O", "S", "Cl", "Br", "P", "I")
  
  # 确保所有元素列都是数值型
  mol1_filtered[element_cols] <- lapply(mol1_filtered[element_cols], function(x) {
    as.numeric(replace(x, is.na(x), 0))
  })
  mol2_filtered[element_cols] <- lapply(mol2_filtered[element_cols], function(x) {
    as.numeric(replace(x, is.na(x), 0))
  })
  reaction_delta[element_cols] <- lapply(reaction_delta[element_cols], function(x) {
    as.numeric(replace(x, is.na(x), 0))
  })
  
  message("💾 Step 5: Creating database and loading data...")
  
  # 创建数据库连接（内存或磁盘）
  if (use_memory_db) {
    con <- dbConnect(RSQLite::SQLite(), ":memory:")
    message("   Using in-memory database for maximum speed")
  } else {
    db_path <- file.path(out_dir, "temp_reactions.db")
    con <- dbConnect(RSQLite::SQLite(), db_path)
    message("   Using disk database: ", db_path)
  }
  
  # 确保断开连接
  on.exit(dbDisconnect(con), add = TRUE)
  
  # 将数据写入数据库
  dbWriteTable(con, "mol1", mol1_filtered, overwrite = TRUE)
  dbWriteTable(con, "mol2", mol2_filtered, overwrite = TRUE)
  dbWriteTable(con, "reactions", reaction_delta, overwrite = TRUE)
  
  # 创建索引以加速查询
  message("🔧 Step 6: Creating indices for fast lookup...")
  # 为元素列创建索引以加速JOIN查询
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol1_formula ON mol1(Formula)")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_formula ON mol2(Formula)")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mol2_elements ON mol2(C, H, N, O, S, Cl, Br, P, I)")
  
  message("🔍 Step 7: Matching reactions between molecules...")
  
  results_list <- list()
  pb <- txtProgressBar(min = 0, max = nrow(reaction_delta), style = 3)
  
  for (i in 1:nrow(reaction_delta)) {
    reaction_name <- reaction_delta$reaction[i]
    
    # 构建SQL查询，直接在数据库中计算转换后的元素组成
    query <- sprintf("
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
    
    # 执行查询
    matches <- dbGetQuery(con, query)
    
    if (nrow(matches) > 0) {
      results_list[[i]] <- matches
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # 合并所有结果
  if (length(results_list) > 0) {
    results <- bind_rows(results_list)
  } else {
    results <- data.frame(Source = character(), Target = character(), 
                         Reaction = character(), stringsAsFactors = FALSE)
    warning("⚠️ No reactions matched. Check your input data and delta definitions.")
  }
  
  message("💾 Step 8: Saving outputs...")
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_file1 <- file.path(out_dir, "network_edge.csv")
  out_file2 <- file.path(out_dir, "reaction_summary.csv")
  
  write_csv(results, out_file1)
  
  if (nrow(results) > 0) {
    reaction_summary <- results %>% 
      count(Reaction, name = "Count") %>%
      arrange(desc(Count))
    write_csv(reaction_summary, out_file2)
    
    message("✅ Done! Found ", nrow(results), " reaction matches")
    message("   Results saved to:")
    message("   - ", out_file1)
    message("   - ", out_file2)
  } else {
    message("⚠️ No matches found, output files created but empty")
  }
  
  message("📚 Reference: 10.1016/j.watres.2020.116484")
  
  return(invisible(results))
}