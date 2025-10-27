
# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c(
  # Column names used in dplyr/tidyr operations
  "Mass", "MolForm", "intensity", "C", "H", "O", "N", "S", "P",
  "Mass.x", "Mass.y", "PMD", "Trans_Name", "Count", "Error",
  "Diff_Formula", "CAR", "avg_CAR", "Total_Count",
  "status", "value", "fc", "intensity_before", "intensity_after",
  "Formula", "abundance_1", "abundance_2",
  "C_1", "H_1", "O_1", "N_1", "S_1", "P_1",
  "C_2", "H_2", "O_2", "N_2", "S_2", "P_2",
  "Cl", "Br", "I_1", "Cl_1", "Br_1", "Cl_2", "Br_2", "I_2",
  "pre", "pro", "Reaction", "Mass_1", "Mass_2", "Mass_Error",
  "intensity_1", "intensity_2"
))
