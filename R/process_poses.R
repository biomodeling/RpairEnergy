#' Preprocess biological structure (pdb)
#'
#' This function standardizes atom names, identifies N- and C-terminal atoms,
#' and merges the pose data with force field parameters.
#'
#' @param df_structure data.frame or data.table with columns:
#'   id_col, elety, resid, resno, insert, chain, x, y, z
#' @param id_col character string specifying the column used as frame ID (default "pdb_id")
#'
#' @return data.table merged with force field parameters, ready for energy calculations
#' @export
#' @import data.table
preprocess_structure <- function(df_structure, id_col) {

  # Check that id_col is provided
  if (missing(id_col) || !nzchar(id_col)) {
    stop("Error: `id_col` must be specified and cannot be empty.")
  }

  # Required columns
  required_cols <- c(id_col, "elety", "resid", "resno", "insert", "chain", "x", "y", "z")
  missing_cols <- setdiff(required_cols, colnames(df_structure))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing in df_structure: ", paste(missing_cols, collapse = ", "))
  }

  atoms_nter <- c("N", "H1", "H2", "H3", "HT1", "HT2", "HT3", "CA", "HA")
  atoms_cter <- c("C", "OT1", "OT2")

  # Convert to data.table if necessary
  if (!is.data.table(df_structure)) df_structure <- as.data.table(df_structure)

  # Load force field parameters
  if (!is.data.table(df_params)) df_params <- as.data.table(df_params)

  # to avoid note in devtools::check()
  is_nter <- NULL
  elety <- NULL
  ..id_col <- NULL
  is_cter <- NULL
  charmm_resid <- NULL
  resid <- NULL
  residue <- NULL
  resno <- NULL
  insert <- NULL
  chain <- NULL
  . <- NULL

  df_structure_prepared <- df_structure |>
    _[, is_nter := any(elety %in% c("H1", "HT1", "H2", "HT2", "H3", "HT3")), by = ..id_col] |>
    _[, is_cter := any(elety %in% c("OT1", "OT2")), by = ..id_col] |>
    _[, charmm_resid := fcase(
      is_nter == TRUE & elety %in% atoms_nter, "NTER",
      is_cter == TRUE & elety %in% atoms_cter, "CTER",
      rep(TRUE, .N), resid
    )] |>
    _[charmm_resid == "NTER",
      elety := fcase(
        elety == "H1", "HT1",
        elety == "H2", "HT2",
        elety == "H3", "HT3",
        rep(TRUE, .N), elety
      )] |>
    _[, residue := paste(resid, resno, insert, chain, sep = '_')] |>
    _[, `:=` (is_nter = NULL, is_cter = NULL)] |>
    _[elety != "OXT"] |>
    _[charmm_resid == 'HIS',
      charmm_resid := ifelse(any(elety == "HD1") & any(elety == "HE2"), "HSP",
                             ifelse(any(elety == "HD1"), "HSD",
                                    ifelse(any(elety == "HE2"), "HSE", charmm_resid))),
      by = .(residue, get(id_col))
    ]

  df_merged <- merge(df_structure_prepared, df_params, by.x = c('resid', 'elety'), by.y = c('resid', 'atom'))

  return(df_merged)
}
