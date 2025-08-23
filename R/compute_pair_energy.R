#' @useDynLib RpairEnergy, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import data.table
NULL

#' Compute pairwise energy using C++ cell list
#'
#' This function wraps the C++ implementation for nearest-neighbor energy calculations.
#'
#' @param dt_s data.table or data.frame with the structure prepared by preprocess_structure()
#' @param cutoff numeric, cutoff distance in Angstroms for neighbor search
#'
#' @return data.table with pairwise energies (LJ and Coulomb) for each residue pair
#' @export
compute_pair_energy <- function(dt_s, cutoff = 12.0) {
  # Ensure dt_s is data.table
  if (!is.data.table(dt_s)) dt_s <- as.data.table(dt_s)

  # Check required columns
  required_cols <- c("x", "y", "z", "eps", "r_min2", "q", "residue")
  missing_cols <- setdiff(required_cols, colnames(dt_s))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in input dt_s: ", paste(missing_cols, collapse = ", "))
  }

  # Load Rcpp if not already
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package 'Rcpp' is required for compute_pair_energy.")
  }

  # Call the C++ function
  # Assumes compute_pair_energy_cpp is exported via Rcpp::export in the C++ file
  res <- compute_pair_energy_cpp(dt_s, cutoff)

  # Return as data.table
  return(as.data.table(res))
}
