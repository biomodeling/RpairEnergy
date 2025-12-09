#' Compute Lennard-Jones interface energy between two chain groups
#'
#' This function computes the average and total Lennard-Jones energy between
#' two groups of chains (g1 and g2) across residue-residue interactions,
#' optionally filtered by a distance cutoff.
#'
#' @param dt_pair_energy A data.table containing pairwise residue energy data,
#'                       typically from compute_pair_energy().
#' @param cutoff Numeric distance cutoff (Å) to filter interactions (e.g., 8.5).
#' @param id_col String. Column name used to group results (e.g., `"pdb_id"`).
#' @param g1 Character vector of chain IDs for group 1. If `NA`, inferred automatically.
#' @param g2 Character vector of chain IDs for group 2. If `NA`, inferred automatically.
#'
#' @details
#' If both `g1` and `g2` are not provided (i.e., `NA`), the function will attempt
#' to automatically infer the two groups of chains from the data. This requires
#' exactly two unique chains to be present in the input, otherwise an error is thrown.
#'
#' Chains in `g1` and `g2` must be disjoint and present in the input data.
#' Interactions within the same group are excluded.
#'
#' @return A data.table with one row per unique value in `id_col`, containing:
#' \itemize{
#'   \item `lj_avg`: Average Lennard-Jones energy at the interface
#'   \item `lj_interface`: Total Lennard-Jones energy at the interface
#'   \item `n_interface`: Number of residue-residue interactions at the interface
#' }
#'
#' @examples
#' \dontrun{
#' # With explicit groups
#' get_lennard_jones(dt, cutoff = 8.5, id_col = "pdb_id", g1 = c("H", "L"), g2 = c("A"))
#'
#' # With automatic group detection (only if exactly 2 unique chains)
#' get_lennard_jones(dt, cutoff = 8.5, id_col = "pdb_id")
#' }
#'
#' @export
get_lennard_jones <- function(dt_pair_energy, cutoff, id_col, g1 = NA, g2 = NA) {
  df <- copy(dt_pair_energy)

  df_energy_inter <- df[r_avg < cutoff & resid_i != resid_j] |>
    _[, c("resid_clean_i", "resno_i", "insert_i", "ch_i") := tstrsplit(resid_i, "_", fixed = TRUE)] |>
    _[, c("resid_clean_j", "resno_j", "insert_j", "ch_j") := tstrsplit(resid_j, "_", fixed = TRUE)]

  if (all(is.na(g1)) || all(is.na(g2))) {
    all_chains <- unique(c(df_energy_inter$ch_i, df_energy_inter$ch_j))
    if (length(all_chains) != 2) {
      stop(paste0("Cannot deduce g1 and g2: expected exactly 2 unique chains, found: ", paste(all_chains, collapse = ", ")))
    }
    g1 <- all_chains[1]
    g2 <- all_chains[2]
  }

  all_chains <- unique(c(df_energy_inter$ch_i, df_energy_inter$ch_j))

  missing_g1 <- setdiff(g1, all_chains)
  missing_g2 <- setdiff(g2, all_chains)

  if (length(missing_g1) > 0 || length(missing_g2) > 0) {
    stop(paste0(
      "Some chains in g1 or g2 are not present in the data:\n",
      if (length(missing_g1) > 0) paste0(" - Missing in g1: ", paste(missing_g1, collapse = ", "), "\n") else "",
      if (length(missing_g2) > 0) paste0(" - Missing in g2: ", paste(missing_g2, collapse = ", "), "\n") else ""
    ))
  }

  df_lj <- df_energy_inter |>
    _[, group_i := fcase(ch_i %in% g1, "g1", ch_i %in% g2, "g2")] |>
    _[, group_j := fcase(ch_j %in% g1, "g1", ch_j %in% g2, "g2")] |>
    _[!is.na(group_i) & !is.na(group_j) & group_i != group_j] |>
    _[, .(lj_avg = mean(lj_res), lj_interface = sum(lj_res), n_interface = .N), by = c(id_col)]

  return(df_lj)
}


#' Compute Graph Metrics from Pairwise Energy Data
#'
#' Given a subset of pairwise energy data for a single structure or frame,
#' this function constructs an undirected residue interaction graph and computes
#' network metrics for each residue (node).
#'
#' The graph can be weighted using the absolute Lennard-Jones energy values (`|lj_res|`)
#' or treated as unweighted (binary edges).
#'
#' @param df_energy_sub A `data.table` containing pairwise interaction data,
#'   typically filtered from a larger dataset by structure/frame.
#'   Must contain the columns: `resid_i`, `resid_j`, and optionally `lj_res` if `weighted = TRUE`.
#' @param weighted Logical. If `TRUE` (default), edges are weighted by `abs(lj_res)`.
#'   If `FALSE`, edges are unweighted.
#'
#' @return A `data.table` with one row per residue (node), including the following metrics:
#' \describe{
#'   \item{vertex}{Residue identifier.}
#'   \item{strength}{Sum of edge weights (or degree, if unweighted).}
#'   \item{closeness}{Closeness centrality.}
#'   \item{betweenness}{Betweenness centrality.}
#'   \item{clustering_coefficient}{Local clustering coefficient.}
#'   \item{degree}{Node degree.}
#' }
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame strength closeness betweenness transitivity degree V E
#'
#' @keywords internal
#' @noRd
compute_graph_metrics <- function(df_energy_sub, weighted = TRUE) {
  if (weighted) {
    df_edges <- df_energy_sub[lj_res < 0, .(resid_i, resid_j, weight = abs(lj_res))]
  } else {
    df_edges <- unique(df_energy_sub[, .(resid_i, resid_j)])
  }

  df_edges <- as.data.frame(df_edges)
  g_vertices <- unique(c(df_edges$resid_i, df_edges$resid_j))

  graph <- if (weighted) {
    igraph::graph_from_data_frame(df_edges, directed = FALSE, vertices = g_vertices)
  } else {
    igraph::graph_from_data_frame(df_edges, directed = FALSE, vertices = g_vertices)
  }

  df_graph <- data.table(
    vertex = names(V(graph)),
    strength = if (weighted) igraph::strength(graph, weights = E(graph)$weight) else igraph::degree(graph),
    closeness = igraph::closeness(graph, weights = if (weighted) 1/E(graph)$weight else NULL),
    betweenness = igraph::betweenness(graph, weights = if (weighted) 1/E(graph)$weight else NULL, directed = FALSE),
    clustering_coefficient = igraph::transitivity(graph, type = "local"),
    degree = igraph::degree(graph)
  )
  return(df_graph)
}


#' Compute Network Parameters for Residue Interaction Graphs
#'
#' Constructs residue interaction graphs from a pairwise energy data.table and computes
#' various graph-based centrality measures for each residue. Optionally, returns
#' metrics summarized only for residues at the interface between two chain groups.
#'
#' @param dt_pair_energy A `data.table` containing pairwise interaction data, typically
#'   produced by `compute_pair_energy()`. Must include columns: `resid_i`, `resid_j`, `r_avg`,
#'   and `lj_res`.
#' @param cutoff Numeric. Distance cutoff (in Ångström) for considering interactions.
#'   Only pairs with `r_avg < cutoff` are included.
#' @param id_col Character string. Name of the column used to group structures (e.g., `"pdb_id"`).
#' @param g1 Character vector of chain IDs representing group 1. If not provided (default `NA`),
#'   the function will attempt to infer groups automatically, assuming exactly two chains.
#' @param g2 Character vector of chain IDs representing group 2. If not provided (default `NA`),
#'   the function will attempt to infer groups automatically, assuming exactly two chains.
#' @param interface_only Logical. If `TRUE` (default), returns mean centrality metrics
#'   for residues at the interface (i.e., those interacting between g1 and g2).
#'   If `FALSE`, returns centrality metrics for all residues in each structure.
#' @param weighted Logical. If `TRUE` (default), edges are weighted by `abs(lj_res)`.
#'   If `FALSE`, the interaction graph is unweighted.
#' @param progress Logical. If `TRUE`, displays a simple text-based progress bar.
#'
#' @details
#' The function processes each structure (identified by `id_col`) independently.
#' For each, it builds an interaction graph where nodes are residues and edges represent
#' non-covalent interactions below the given distance cutoff. Centrality measures
#' (e.g., strength, closeness, betweenness) are computed using the `igraph` package.
#'
#' If `interface_only = TRUE`, the function identifies interface residues as those
#' involved in interactions between chains in `g1` and `g2`, and averages the
#' network metrics over these residues.
#'
#' If `g1` and `g2` are not specified, the function will try to infer them from the data,
#' assuming exactly two unique chain identifiers are present.
#'
#' @return A `data.table` or a `list`:
#' \describe{
#'   \item{If `interface_only = TRUE`}{
#'     Returns a `data.table` with one row per structure and mean values of
#'     network metrics computed over:
#'     \itemize{
#'       \item `*_all`: all interface residues (g1 + g2)
#'       \item `*_g1`: interface residues from group 1
#'       \item `*_g2`: interface residues from group 2
#'     }
#'   }
#'   \item{If `interface_only = FALSE`}{
#'     Returns a `list` with two elements:
#'     \itemize{
#'       \item `full`: a `data.table` of per-residue network metrics for all structures
#'       \item `interface`: a `data.table` of mean interface metrics, as described above
#'     }
#'   }
#' }
#'
#' @seealso [compute_pair_energy()], [compute_graph_metrics()], [igraph]
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame strength closeness betweenness transitivity degree V E
#' @export
get_network_params <- function(
    dt_pair_energy, cutoff, id_col, g1 = NA, g2 = NA, interface_only = TRUE, weighted = TRUE, progress = FALSE) {
  if (!id_col %in% names(dt_pair_energy)) {
    stop("Error: `id_col` must be a column in `dt_pair_energy`.")
  }

  df_energy_for_network <- dt_pair_energy[r_avg < cutoff & resid_i != resid_j]
  unique_ids <- unique(df_energy_for_network[[id_col]])

  out_list <- vector("list", length(unique_ids))

  n_ids <- length(unique_ids)
  bar_len <- 40

  for (i in seq_along(unique_ids)) {
    current_id <- unique_ids[[i]]
    df_sub <- df_energy_for_network[get(id_col) == current_id]

    graph_params <- compute_graph_metrics(df_sub, weighted = weighted)
    graph_params[[id_col]] <- current_id

    out_list[[i]] <- graph_params

    # PROGRESS BAR
    if (progress && (i %% max(1, floor(n_ids / 100)) == 0 || i == n_ids)) {
      filled_len <- round(i / n_ids * bar_len)
      cat("\r[", paste(rep("=", filled_len), collapse = ""),
          paste(rep(" ", bar_len - filled_len), collapse = ""),
          "] ", sprintf("%3d", round(i / n_ids * 100)), "% (", i, "/", n_ids, ")", sep = "")
      flush.console()
      if (i == n_ids) cat("\n")
    }
  }

  df_network_params <- rbindlist(out_list)

  # Prepare chain columns
  df_energy_for_network <- df_energy_for_network |>
    _[, c("resid_clean_i", "resno_i", "insert_i", "ch_i") := tstrsplit(resid_i, "_", fixed = TRUE)] |>
    _[, c("resid_clean_j", "resno_j", "insert_j", "ch_j") := tstrsplit(resid_j, "_", fixed = TRUE)]

  # Infer g1/g2 if not provided
  if (all(is.na(g1)) || all(is.na(g2))) {
    all_chains <- unique(c(df_energy_for_network$ch_i, df_energy_for_network$ch_j))
    if (length(all_chains) != 2) {
      stop(paste0("Cannot deduce g1 and g2: expected exactly 2 unique chains, found: ", paste(all_chains, collapse = ", ")))
    }
    g1 <- all_chains[1]
    g2 <- all_chains[2]
  }

  # Check for chain presence
  all_chains <- unique(c(df_energy_for_network$ch_i, df_energy_for_network$ch_j))
  missing_g1 <- setdiff(g1, all_chains)
  missing_g2 <- setdiff(g2, all_chains)
  if (length(missing_g1) > 0 || length(missing_g2) > 0) {
    stop(paste0(
      "Some chains in g1 or g2 are not present in the data:\n",
      if (length(missing_g1) > 0) paste0(" - Missing in g1: ", paste(missing_g1, collapse = ", "), "\n") else "",
      if (length(missing_g2) > 0) paste0(" - Missing in g2: ", paste(missing_g2, collapse = ", "), "\n") else ""
    ))
  }

  df_interface_residues_expanded <- df_energy_for_network |>
    _[(ch_i %in% g1 & ch_j %in% g2) | (ch_i %in% g2 & ch_j %in% g1)] |>
    _[, .(all_res = list(unique(c(resid_i, resid_j))),
          g1_res = list(unique(c(resid_i[ch_i %in% g1], resid_j[ch_j %in% g1]))),
          g2_res = list(unique(c(resid_i[ch_i %in% g2], resid_j[ch_j %in% g2])))), by = id_col]

  graph_params <- c("strength", "closeness", "betweenness", "clustering_coefficient", "degree")

  df_all <- merge(df_network_params, df_interface_residues_expanded, by = id_col, allow.cartesian = TRUE) |>
    _[vertex %in% unlist(all_res), lapply(.SD, mean, na.rm = TRUE), by = id_col, .SDcols = graph_params] |>
    setnames(old = graph_params, new = paste0(graph_params, "_all"))

  df_g1 <- merge(df_network_params, df_interface_residues_expanded, by = id_col, allow.cartesian = TRUE) |>
    _[vertex %in% unlist(g1_res), lapply(.SD, mean, na.rm = TRUE), by = id_col, .SDcols = graph_params] |>
    setnames(old = graph_params, new = paste0(graph_params, "_g1"))

  df_g2 <- merge(df_network_params, df_interface_residues_expanded, by = id_col, allow.cartesian = TRUE) |>
    _[vertex %in% unlist(g2_res), lapply(.SD, mean, na.rm = TRUE), by = id_col, .SDcols = graph_params] |>
    setnames(old = graph_params, new = paste0(graph_params, "_g2"))

  df_interface_summary <- Reduce(function(x, y) merge(x, y, by = c(id_col)), list(df_all, df_g1, df_g2))

  if (identical(interface_only, TRUE)) {
    return(df_interface_summary)
  } else {
    return(list(
      full = df_network_params,
      interface = df_interface_summary
    ))
  }
}
