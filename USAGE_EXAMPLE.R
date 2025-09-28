
# 使用新添加的质量差异匹配函数

library(ubynet)

# 查看函数帮助
?match_reactions_by_mass_difference

# 创建示例数据
# 分子文件1 (inflow)
mol1 <- data.frame(
  MolForm = c("C6H12O6", "C12H22O11", "C2H4O2"),
  Mass = c(180.0634, 342.1162, 60.0211),
  intensity = c(1000, 500, 800)
)
write.csv(mol1, "example_mol1.csv", row.names = FALSE)

# 分子文件2 (outflow)  
mol2 <- data.frame(
  MolForm = c("C6H12O6", "C12H22O11", "C2H4O2", "C6H10O6"),
  Mass = c(180.0634, 342.1162, 60.0211, 178.0477),
  intensity = c(200, 1500, 400, 300)
)
write.csv(mol2, "example_mol2.csv", row.names = FALSE)

# 反应数据
reactions <- data.frame(
  reaction = c("Dehydration", "Oxidation"),
  mass_difference = c(-18.0106, -2.0156)
)
write.csv(reactions, "example_reactions.csv", row.names = FALSE)

# 运行分析
results <- match_reactions_by_mass_difference(
  file1 = "example_mol1.csv",
  file2 = "example_mol2.csv", 
  reaction_delta_file = "example_reactions.csv",
  out_dir = "results"
)

# 查看结果
print(results$reaction_summary)
head(results$network_edges)

