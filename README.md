# RpairEnergy

[![Version](man/figures/version.svg)](https://github.com/biomodeling/RpairEnergy)
[![License: GPL v3](man/figures/license.svg)](https://github.com/biomodeling/RpairEnergy/blob/master/LICENSE.md)

`RpairEnergy` is an R package designed to compute pairwise energies (Coulomb and Lennard-Jones) for biological systems using the CHARMM36 force field parameters. The package leverages `Rcpp` for efficient computation and `data.table` for fast data manipulation.

## Functions

### `preprocess_structure`

Preprocesses biological structure data (PDB format) by standardizing atom names, identifying N- and C-terminal atoms, and merging the pose data with force field parameters.

#### Usage

```R
preprocess_structure(df_structure, id_col = "pdb_id")
```

#### Arguments

* `df_structure`: A data.frame or data.table with the following columns:

  * `id_col`: Frame ID (default is `"pdb_id"`)
  * `elety`: Atom type
  * `resid`: Residue ID
  * `resno`: Residue number
  * `insert`: Insert code (optional)
  * `chain`: Chain identifier
  * `x, y, z`: Atomic coordinates (in Angstroms)

* `id_col`: Character string specifying the column used as frame ID (default is `"pdb_id"`).

#### Returns

A `data.table` merged with force field parameters, ready for energy calculations.

---

### `compute_pair_energy`

Computes pairwise energies (Lennard-Jones and Coulomb) for a structure using a C++ implementation with a cell list for efficient nearest-neighbor search.

#### Usage

```R
compute_pair_energy(dt_s, cutoff = 12.0)
```

#### Arguments

* `dt_s`: A data.table or data.frame with the structure prepared by `preprocess_structure()`.
* `cutoff`: Numeric value specifying the cutoff distance (in Angstroms) for the neighbor search (default is `12.0`).

#### Returns

A `data.table` containing the pairwise energies (LJ and Coulomb) for each residue pair.

---

### `get_lennard_jones`

Calculates the Lennard-Jones energy at the interface between two groups of chains, where the interface is defined by the interactions between residues of the two groups, filtered by a distance cutoff.

#### Usage

```R
get_lennard_jones(dt, cutoff = 8.5, id_col = "pdb_id", g1 = c("H", "L"), g2 = "A")
```

#### Arguments

* `dt`: `data.table` with pairwise residue energy data obtained from `compute_pair_energy`.
* `cutoff`: Numeric distance cutoff (Å) to filter interactions. Default is 8.5.
* `id_col`: Column name to group results (e.g., `"pdb_id"`).
* `g1`: Character vector of chain IDs for group 1 (optional).
* `g2`: Character vector of chain IDs for group 2 (optional).

#### Returns

A `data.table` with:
* `id_col`: Frame ID.
* `lj_avg`: Average Lennard-Jones energy at the interface.
* `lj_interface`: Total Lennard-Jones energy at the interface.
* `n_interface`: Number of residue-residue interactions.

---

### `get_network_params`

Computes network parameters (e.g., centrality measures) for residue interaction graphs, optionally summarizing results only for residues at the interface between two chain groups, filtered by a distance cutoff.

#### Usage

```R
get_network_params(dt, cutoff = 8.5, id_col = "pdb_id", g1 = c("H", "L"), g2 = "A")
```

#### Arguments

* `dt`: `data.table` with pairwise residue interaction data (e.g., from `compute_pair_energy`).
* `cutoff`: Numeric distance cutoff (Å). Default is 8.5.
* `id_col`: Column name to group results (e.g., `"pdb_id"`).
* `g1`: Character vector of chain IDs for group 1 (optional).
* `g2`: Character vector of chain IDs for group 2 (optional).
* `interface_only`: Logical. If `TRUE` (default), returns metrics for residues at the interface (g1 + g2). If `FALSE`, returns metrics for all residues.
* `weighted`: Logical. If `TRUE` (default), edges are weighted by `abs(lj_res)`. If `FALSE`, the graph is unweighted.
* `progress`: Logical. If `TRUE`, shows a progress bar.

#### Returns

* If `interface_only = TRUE`: A `data.table` with mean network metrics for interface residues (g1 + g2).
* If `interface_only = FALSE`: A `list` containing:

  * `full`: `data.table` of per-residue network metrics.
  * `interface`: `data.table` of mean interface metrics.

---

### `get_statistical_potentials`

Computes the sum of symmetric and asymmetric statistical potentials for antibody-antigen interface contacts, based on contact data, optionally filtered by specific parts (e.g., CDR loops).

#### Usage

```R
get_statistical_potentials(df_contacts, parts = 'all', id_col = 'pdb_id')
```

#### Arguments

* `df_contacts`: `data.table` or `data.frame` containing contact information with columns: `resid_ab`, `resid_ag`, and `id_col` (e.g., `"pdb_id"`). Should be pre-filtered by a distance threshold (e.g., distance < 8.5 Å).
* `parts`: Character vector specifying which parts of the protein to calculate potentials for. Default is `"all"`, but can also include specific CDR loops like `'h3'`, `'l3'`, etc.
* `id_col`: The column name used to group results (e.g., `"pdb_id"`). Default is `"pdb_id"`.

#### Returns

A `data.table` with summed symmetric (`sp_sym_<part>`) and asymmetric (`sp_asym_<part>`) statistical potentials for each identifier, grouped by the specified `id_col`.

---

### `get_triangles`

Computes the total log-potential score for antibody-antigen residue triangles, formed by a base pair and a third residue, matching them against a reference table of triangle potentials.

#### Usage

```R
get_triangles_potentials(df_contacts, base_side, col_id, split_col_id)
```

#### Arguments

* `df_contacts`: `data.table` containing contact information and Cartesian coordinates for residues.
* `base_side`: Character string, `"ab"` or `"ag"`, selecting which residue pair (antibody or antigen) forms the base of the triangle.
* `col_id`: Column name used to identify distinct models or structural units (e.g., `"model"`).
* `split_col_id`: Logical. If `TRUE`, computes triangles separately for each unique `col_id`. If `FALSE`, computes using all rows together.

#### Returns

A `data.table` with:

* `col_id`: Frame ID (e.g., `"pdb_id"`).
* `score`: Total sum of log-potentials for all antibody-antigen triangles in that model.

---

## Dependencies

* `data.table`: For efficient data manipulation.
* `Rcpp`: For wrapping the C++ implementation for energy calculations.
* `igraph`: For graph-based centrality measures and network analysis.
* `pbapply`: For progress bars in batch processing.

---

## Installation

You can install the latest version of the package directly from the GitHub repository:

```R
devtools::install_github("biomodeling/RpairEnergy")
```

Make sure you have the `devtools` package installed. If you don’t have it yet, install it using:

```R
install.packages("devtools")
```

---

## Example

```R
library(RpairEnergy)

# Load your structure data
structure_data <- fread("your_structure_data.csv")

# Preprocess the structure
preprocessed_data <- preprocess_structure(structure_data)

# Compute pairwise energies
pairwise_energies <- compute_pair_energy(preprocessed_data)
```
