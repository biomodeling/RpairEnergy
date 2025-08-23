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

## Dependencies

* `data.table`: For efficient data manipulation.
* `Rcpp`: For wrapping the C++ implementation for energy calculations.

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

