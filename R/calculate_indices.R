#' Calculate Chemical Indices from Molecular Formula Elements
#'
#' This function calculates various chemical indices including aromaticity (AI, AI_Mod),
#' double bond equivalents (DBE, DBE_O), Gibbs free energy (GFE), Kendrick mass and defect,
#' nominal oxidation state of carbon (NOSC), and elemental ratios (O:C, H:C, P:C, N:C, N:P).
#'
#' @param input_file Character string. Path to the input CSV file containing molecular formula elements.
#'   Default is "merged_molform_elements.csv".
#' @param output_file Character string. Path to save the output CSV file with calculated indices.
#'   Default is "chemical_indices_results.csv".
#' @param kendrick_base Character string. Base for Kendrick mass calculation.
#'   Options: "CH2" (default), "O", "H2", "COO", "H2O".
#'
#' @return A data frame containing the original data plus calculated chemical indices.
#'
#' @details
#' The function calculates the following indices:
#' \itemize{
#'   \item \strong{AI}: Aromaticity Index = (1 + C - O - S - 0.5*(N + P + H)) / (C - O - S - N - P)
#'   \item \strong{AI_Mod}: Modified Aromaticity Index = (1 + C - 0.5*O - S - 0.5*(N + P + H)) / (C - 0.5*O - S - N - P)
#'   \item \strong{DBE}: Double Bond Equivalent = 1 + C - 0.5*(N + P + H)
#'   \item \strong{DBE_O}: DBE minus Oxygen = 1 + C - 0.5*(N + P + H) - O
#'   \item \strong{GFE}: Cox Gibbs Free Energy = 60.3 - 28.5*(N + O + S)/C
#'   \item \strong{Kendrick_Mass}: (Observed Mass) * (Nominal Mass of base / Exact Mass of base)
#'   \item \strong{Kendrick_Defect}: ceiling(Kendrick_Mass) - Kendrick_Mass
#'   \item \strong{NOSC}: Nominal Oxidation State of Carbon = -((4*C + H - 3*N - 2*O + 5*P - 2*S) / C) + 4
#'   \item \strong{O_C, H_C, P_C, N_C, N_P}: Elemental ratios
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- calculate_chemical_indices("merged_molform_elements.csv")
#'
#' # With custom output file
#' results <- calculate_chemical_indices(
#'   input_file = "my_data.csv",
#'   output_file = "my_results.csv"
#' )
#'
#' # With different Kendrick base
#' results <- calculate_chemical_indices(
#'   input_file = "my_data.csv",
#'   kendrick_base = "O"
#' )
#' }
#'
#' @importFrom utils read.csv write.csv
#' @export
calculate_chemical_indices <- function(input_file = "merged_molform_elements.csv",
                                       output_file = "chemical_indices_results.csv",
                                       kendrick_base = "CH2") {
  
  # Read input data
  cat("Reading input file:", input_file, "\n")
  data <- utils::read.csv(input_file, stringsAsFactors = FALSE)
  
  # Check required columns exist
  required_cols <- c("C", "H", "O", "N", "S", "P")
  missing_cols <- setdiff(required_cols, names(data))
  if(length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Replace NA with 0 for element counts
  data[required_cols][is.na(data[required_cols])] <- 0
  
  # ===== calc_aroma: Aromaticity indices =====
  cat("Calculating aromaticity indices...\n")
  data$AI_denominator <- data$C - data$O - data$S - data$N - data$P
  data$AI <- ifelse(data$AI_denominator != 0,
                    (1 + data$C - data$O - data$S - 0.5 * (data$N + data$P + data$H)) / data$AI_denominator,
                    NA)
  
  data$AI_Mod_denominator <- data$C - 0.5 * data$O - data$S - data$N - data$P
  data$AI_Mod <- ifelse(data$AI_Mod_denominator != 0,
                        (1 + data$C - 0.5 * data$O - data$S - 0.5 * (data$N + data$P + data$H)) / data$AI_Mod_denominator,
                        NA)
  
  # ===== calc_dbe: Double Bond Equivalent =====
  cat("Calculating DBE indices...\n")
  data$DBE <- 1 + data$C - 0.5 * (data$N + data$P + data$H)
  data$DBE_O <- 1 + data$C - 0.5 * (data$N + data$P + data$H) - data$O
  
  # ===== calc_gibbs: Cox Gibbs Free Energy =====
  cat("Calculating Gibbs Free Energy...\n")
  data$NOSC_temp <- data$N + data$O + data$S
  data$GFE <- ifelse(data$C != 0,
                     60.3 - 28.5 * (data$NOSC_temp / data$C),
                     NA)
  
  # ===== calc_kendrick: Kendrick Mass and Defect =====
  cat("Calculating Kendrick mass and defect...\n")
  
  # Define exact masses for common Kendrick bases
  kendrick_bases <- list(
    "CH2" = list(nominal = 14, exact = 14.01565),
    "O" = list(nominal = 16, exact = 15.99491),
    "H2" = list(nominal = 2, exact = 2.01565),
    "COO" = list(nominal = 44, exact = 43.98983),
    "H2O" = list(nominal = 18, exact = 18.01056)
  )
  
  if(!kendrick_base %in% names(kendrick_bases)) {
    warning("Kendrick base not recognized, using CH2. Available: ",
            paste(names(kendrick_bases), collapse = ", "))
    kendrick_base <- "CH2"
  }
  
  base_info <- kendrick_bases[[kendrick_base]]
  kendrick_factor <- base_info$nominal / base_info$exact
  
  # Check if observed mass column exists
  mass_col <- NULL
  possible_mass_cols <- c("mass", "Mass", "ObservedMass", "observed_mass", "m/z", "mz")
  for(col in possible_mass_cols) {
    if(col %in% names(data)) {
      mass_col <- col
      break
    }
  }
  
  if(is.null(mass_col)) {
    warning("No mass column found. Kendrick calculations skipped.")
    data$Kendrick_Mass <- NA
    data$Kendrick_Defect <- NA
  } else {
    data$Kendrick_Mass <- data[[mass_col]] * kendrick_factor
    data$Kendrick_Defect <- ceiling(data$Kendrick_Mass) - data$Kendrick_Mass
  }
  
  # ===== calc_nosc: Nominal Oxidation State of Carbon =====
  cat("Calculating NOSC...\n")
  data$NOSC <- ifelse(data$C != 0,
                      -(4 * data$C + data$H - 3 * data$N - 2 * data$O + 5 * data$P - 2 * data$S) / data$C + 4,
                      NA)
  
  # ===== calc_element_ratios: Elemental ratios =====
  cat("Calculating elemental ratios...\n")
  data$O_C <- ifelse(data$C != 0, data$O / data$C, NA)
  data$H_C <- ifelse(data$C != 0, data$H / data$C, NA)
  data$P_C <- ifelse(data$C != 0, data$P / data$C, NA)
  data$N_C <- ifelse(data$C != 0, data$N / data$C, NA)
  data$N_P <- ifelse(data$P != 0, data$N / data$P, NA)
  
  # Remove temporary columns
  data$AI_denominator <- NULL
  data$AI_Mod_denominator <- NULL
  data$NOSC_temp <- NULL
  
  # Write output
  cat("Writing results to:", output_file, "\n")
  utils::write.csv(data, output_file, row.names = FALSE)
  
  # Print summary
  cat("\n===== Calculation Summary =====\n")
  cat("Total rows processed:", nrow(data), "\n")
  cat("Columns added:", "AI, AI_Mod, DBE, DBE_O, GFE, Kendrick_Mass, Kendrick_Defect, NOSC, O_C, H_C, P_C, N_C, N_P\n")
  cat("Output saved to:", output_file, "\n")
  
  return(invisible(data))
}