#' Compute Statistical Potentials for Interface Contacts
#'
#' Computes the sum of statistical potentials (both symmetric and asymmetric)
#' for antibody-antigen interface contacts, based on input data frames.
#'
#' @param df_contacts A data.table or data.frame containing the contact information with at least
#'   three columns: `resid_ab`, `resid_ag`, and `id_col` (default: 'pdb_id'). This data frame
#'   should already be filtered based on distance threshold (e.g., distance < 8.5).
#' @param parts A character vector specifying the parts of the protein for which to calculate the potentials.
#'   Default is 'all', but can also include 'h3', 'l3', or other cdr's loops.
#' @param id_col The name of the column to use as the identifier (e.g., 'pdb_id'). Default is 'pdb_id'.
#'
#' @return A data.table containing the summed statistical potentials for each identifier,
#'   with columns for symmetric (sp_sym_<part>) and asymmetric (sp_asym_<part>) potentials for each part.
#'   The dataframe is grouped by the specified identifier (e.g., 'pdb_id').
#'
#' @import data.table
#' @export
get_statistical_potentials <- function(df_contacts, parts = 'all', id_col = 'pdb_id') {

  complete_sym_matrix <- function(df, part_filter) {
    df1 <- df[part == part_filter, .(resid_ab, resid_ag, potential)]
    df2 <- df1[resid_ab != resid_ag, .(resid_ab2 = resid_ag, resid_ag2 = resid_ab, potential = potential)]
    setnames(df2, c("resid_ab2", "resid_ag2"), c("resid_ab", "resid_ag"))
    rbindlist(list(df1, df2))
  }

  df_sp_interface <- df_contacts[, c(id_col, "resid_ab", "resid_ag"), with = FALSE]

  # run dynamic merges
  for (p in parts) {
    if (p == 'all') {
      sp_sym_col <- 'sp_sym'
      sp_asym_col <- 'sp_asym'
    } else {
      sp_sym_col <- paste0("sp_sym_", p)
      sp_asym_col <- paste0("sp_asym_", p)
    }

    # Merging for sp_sym
    df_sp_interface <- merge(
      df_sp_interface, complete_sym_matrix(df_sp_s, p),
      by = c("resid_ab", "resid_ag"), all.x = TRUE)
    setnames(df_sp_interface, "potential", sp_sym_col)

    # Merging for sp_asym
    df_sp_interface <- merge(
      df_sp_interface, df_sp_a[part == p],
      by = c("resid_ab", "resid_ag"), all.x = TRUE)
    setnames(df_sp_interface, "potential", sp_asym_col)
  }

  # Sum potential for each id_col (default: pdb_id)
  df_sp_interface <- df_sp_interface |>
    _[, lapply(.SD, sum), by = id_col, .SDcols = patterns("^sp_")]

  return(df_sp_interface)
}


#' Construct geometric triangles from residue contacts
#'
#' @description
#' Internal helper function used by `get_triangles_potentials()`.
#' Given a contact table, it identifies all valid base–vertex triangles
#' according to geometric constraints (base pairs, distances, and angles).
#'
#' @param df_contacts A `data.table` containing residue–residue contact data.
#' @param base_side Character string, either `"ab"` or `"ag"`, selecting which
#'   coordinate set is used as the triangle base.
#' @param col_id Character string giving the name of the identifier column
#'   (e.g., model, frame, structure id).
#' @param id Optional. If provided, the function is restricted to rows where
#'   `df_contacts[[col_id]] == id`. If `NULL`, all rows are used.
#'
#' @return A `data.table` containing all valid triangles for that model/id,
#'   including:
#'   - base endpoints (`res_b1`, `res_b2`)
#'   - vertex residue (`res_v`)
#'   - base–vertex distances (`d1`, `d2`)
#'   - base length
#'   - angles (`alpha`, `beta`, `gamma`)
#'
#' @details
#' A triangle is kept only if:
#' - the two base residues are different,
#' - both vertex distances are < 8.5 Å,
#' - angles `beta` and `gamma` are acute.
#'
#' @keywords internal
#' @noRd
get_triangles <- function(df_contacts, base_side, col_id, id = NULL) {
  df_i <- if (is.null(id)) {
    df_contacts
  } else {
    df_contacts[get(col_id) == id]
  }

  if (base_side == "ab") {
    base_coords <- df_i |>
      _[, res_b := paste(resno_ab, resid_ab, insert_ab, chain_ab, sep = '_')] |>
      _[, c(col_id, "res_b", "x_ab", "y_ab", "z_ab"), with = FALSE] |>
      setnames(c("x_ab", "y_ab", "z_ab"), c("x_b", "y_b", "z_b")) |>
      unique()
    df_vertex <- df_i |>
      _[, res_v := paste(resno_ag, resid_ag, insert_ag, chain_ag, sep = '_')] |>
      _[, c(col_id, "res_v", "x_ag", "y_ag", "z_ag"), with = FALSE] |>
      setnames(c("x_ag", "y_ag", "z_ag"), c("x_v", "y_v", "z_v")) |>
      unique()
  } else if (base_side == "ag") {
    base_coords <- df_i |>
      _[, res_b := paste(resno_ag, resid_ag, insert_ag, chain_ag, sep = '_')] |>
      _[, c(col_id, "res_b", "x_ag", "y_ag", "z_ag"), with = FALSE] |>
      setnames(c("x_ag", "y_ag", "z_ag"), c("x_b", "y_b", "z_b")) |>
      unique()
    df_vertex <- df_i |>
      _[, res_v := paste(resno_ab, resid_ab, insert_ab, chain_ab, sep = '_')] |>
      _[, c(col_id, "res_v", "x_ab", "y_ab", "z_ab"), with = FALSE] |>
      setnames(c("x_ab", "y_ab", "z_ab"), c("x_v", "y_v", "z_v")) |>
      unique()
  } else {
    stop("Invalid base_side: choose 'ab' or 'ag'")
  }

  df_base <- merge(
    base_coords, base_coords,
    by = col_id, allow.cartesian = TRUE, suffixes = c("1", "2")
  ) |>
    _[res_b1 != res_b2] |>
    _[, base_distance := sqrt((x_b1 - x_b2)^2 + (y_b1 - y_b2)^2 + (z_b1 - z_b2)^2)]

  df_triangle <- merge(
    df_base,
    df_vertex,
    by = col_id,
    allow.cartesian = TRUE
  ) |>
    _[, `:=`(d1 = sqrt((x_b1 - x_v)^2 + (y_b1 - y_v)^2 + (z_b1 - z_v)^2),
             d2 = sqrt((x_b2 - x_v)^2 + (y_b2 - y_v)^2 + (z_b2 - z_v)^2))] |>
    _[d1 < 8.5 & d2 < 8.5] |>
    _[, `:=`(alpha = acos((d1^2 + d2^2 - base_distance^2) / (2*d1*d2)),
             beta  = acos((base_distance^2 + d1^2 - d2^2) / (2*base_distance*d1)),
             gamma = acos((base_distance^2 + d2^2 - d1^2) / (2*base_distance*d2)))] |>
    _[beta < pi/2 & gamma < pi/2]

  return(df_triangle)
}


#' Compute triangle-based interaction potentials
#'
#' @description
#' For each structure or model, this function identifies all residue triangles
#' formed by an antibody-antigen base pair and a third residue. The function
#' matches these triangles against a reference table of triangle potentials
#' and returns the total log-potential score for the interactions.
#'
#' @importFrom pbapply pblapply
#'
#' @param df_contacts A `data.table` containing contact information and
#'   Cartesian coordinates for the residues involved.
#' @param base_side Character string, `"ab"` or `"ag"`, selecting which residue
#'   pair (from antibody or antigen) forms the base of the triangle.
#' @param col_id Character string; the column used to identify distinct models
#'   or structural units (e.g., `"pdb_id"`).
#' @param split_col_id Logical.
#'   - If `FALSE`, triangles are computed using all rows together.
#'   - If `TRUE`, triangles are computed separately for each unique value in `col_id`.
#'
#' @return A `data.table` with two columns:
#'   - `col_id`: model identifier
#'   - `score`: sum of log-potentials for all triangles in that model
#'
#' @details
#' The function:
#' 1. Generates antibody-antigen residue triangles using `get_triangles()`.
#' 2. Splits residue identifiers and filters to standard 20 amino acids.
#' 3. Normalizes base residue order (alphabetical).
#' 4. Builds triangle keys and merges them with the reference potential tables
#'    (`df_triang_ab` or `df_triang_ag`).
#' 5. Sums log-potentials per model.
#'
#' @seealso
#' `get_triangles()` — geometric triangle construction (internal)
#'
#' @export
get_triangles_potentials <- function(df_contacts, base_side, col_id, split_col_id) {
  if (base_side == 'ab') {
    df_triang <- copy(df_triang_ab)
    df_triang <- df_triang[, log_r := log(r)]
  } else if (base_side == 'ag') {
    df_triang <- copy(df_triang_ag)
    df_triang <- df_triang[, log_r := log(r)]
  } else {
    stop("Invalid base_side: choose 'ab' or 'ag'")
  }

  aa3 <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")

  if (isFALSE(split_col_id)) {
    df_triangles_base <- get_triangles(df_contacts, base_side, col_id)
  } else {
    models <- df_contacts[[col_id]] |> unique()
    df_triangles_base <- rbindlist(
      pblapply(models, function(id) {
        get_triangles(df_contacts, base_side, col_id, id)
      })
    )
  }

  df_pair <- df_triangles_base |>
    _[, c("resno_1", "resid_1", "insert_1", "chain_1") := tstrsplit(res_b1, "_")] |>
    _[, c("resno_2", "resid_2", "insert_2", "chain_2") := tstrsplit(res_b2, "_")] |>
    _[, c("resno_v", "resid_v", "insert_v", "chain_v") := tstrsplit(res_v, "_")] |>
    _[resid_1 %in% aa3 & resid_2 %in% aa3 & resid_v %in% aa3] |>
    _[, `:=`(
      resid_a = pmin(resid_1, resid_2),
      resid_b = pmax(resid_1, resid_2)
    )] |>
    _[, base := paste(resid_a, resid_b, sep = "-")] |>
    _[, c(col_id, "base", "resid_v", "base_distance", "d1", "d2", "alpha", "beta", "gamma"), with = FALSE]

  df_score <- merge(
    df_pair[, triangle := paste(base, resid_v, sep = '--')],
    df_triang,
    by = 'triangle'
  ) |>
    _[, .(score = sum(log_r)), by = col_id]

  if (base_side == 'ab') {
    setnames(df_score, 'score', 'r_ab')
  } else if (base_side == 'ag') {
    setnames(df_score, 'score', 'r_ag')
  } else {
    stop("Invalid base_side: choose 'ab' or 'ag'")
  }

  return(df_score)
}
