utils::globalVariables(c(
  ".", "flush.console", "lj_res", "r_avg", "resid_i", "resid_j", "vertex",
  "unique_residues", "group_i", "group_j", "ch_i", "ch_j", "part",
  "all_res", "g1_res", "g2_res", "resid_ab", "resid_ag", "potential",
  # new data.table variables from get_triangles()
  "res_b", "resno_ab", "resid_ab", "insert_ab", "chain_ab",
  "res_v", "resno_ag", "insert_ag", "chain_ag",
  "res_b1", "res_b2", "base_distance",
  "x_b1", "x_b2", "y_b1", "y_b2", "z_b1", "z_b2",
  "x_v", "y_v", "z_v", "d1", "d2",
  # new variables from get_triangles_potentials()
  "log_r", "r",
  "resid_1", "resid_2", "resid_v",
  "resid_a", "resid_b",
  "base", "triangle"
))

