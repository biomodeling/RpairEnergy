#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <string>
#include <iostream>

// Atom structure
struct Atom {
    double x, y, z;
    double eps;      // Lennard-Jones epsilon
    double r_min2;   // Lennard-Jones r_min^2
    double q;        // charge
    std::string resid; // residue identifier
};

// Pairwise energy result
struct PairEnergy {
    std::string resid_i;
    std::string resid_j;
    double lj_res;
    double coulomb_res;
    double r_avg;
};

// Cell list container
struct Cell {
    std::vector<int> atom_indices;
};

// Distance function
inline double distance(const Atom& a, const Atom& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// Compute pair energies using cell list
std::vector<PairEnergy> compute_pair_energy_core(
    const std::vector<Atom>& atoms,
    double cutoff,           // cutoff distance for neighbors
    double eps0 = 8.85418782e-12,  // C^2 N^-1 m^-2
    double qe = -1.602176634e-19,  // C
    double e_conv = 1.439327774e20, // conversion factor
    double d_conv = 1e-10          // Angstrom to m
) {
    if (atoms.empty()) throw std::runtime_error("Atoms vector is empty");

    double k = 1.0 / (4.0 * M_PI * eps0);

    // Determine bounding box
    double xmin = atoms[0].x, xmax = atoms[0].x;
    double ymin = atoms[0].y, ymax = atoms[0].y;
    double zmin = atoms[0].z, zmax = atoms[0].z;

    for (const auto& a : atoms) {
        xmin = std::min(xmin, a.x); xmax = std::max(xmax, a.x);
        ymin = std::min(ymin, a.y); ymax = std::max(ymax, a.y);
        zmin = std::min(zmin, a.z); zmax = std::max(zmax, a.z);
    }

    // Compute number of cells along each dimension
    int nx = std::max(1, int((xmax - xmin) / cutoff) + 1);
    int ny = std::max(1, int((ymax - ymin) / cutoff) + 1);
    int nz = std::max(1, int((zmax - zmin) / cutoff) + 1);

    // Initialize 3D grid of cells
    std::vector<std::vector<std::vector<Cell>>> grid(nx, std::vector<std::vector<Cell>>(ny, std::vector<Cell>(nz)));

    // Assign atoms to cells
    for (size_t i = 0; i < atoms.size(); ++i) {
        int ix = std::min(nx-1, int((atoms[i].x - xmin)/cutoff));
        int iy = std::min(ny-1, int((atoms[i].y - ymin)/cutoff));
        int iz = std::min(nz-1, int((atoms[i].z - zmin)/cutoff));
        grid[ix][iy][iz].atom_indices.push_back(i);
    }

    // Compute pairwise energies
    std::vector<PairEnergy> results;
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                Cell& cell = grid[ix][iy][iz];

                // Loop over neighbor cells (including self)
                for (int dx = -1; dx <= 1; ++dx) {
                    int jx = ix + dx; if (jx < 0 || jx >= nx) continue;
                    for (int dy = -1; dy <= 1; ++dy) {
                        int jy = iy + dy; if (jy < 0 || jy >= ny) continue;
                        for (int dz = -1; dz <= 1; ++dz) {
                            int jz = iz + dz; if (jz < 0 || jz >= nz) continue;

                            Cell& neighbor = grid[jx][jy][jz];

                            for (int i_idx : cell.atom_indices) {
                                for (int j_idx : neighbor.atom_indices) {
                                    if (i_idx >= j_idx) continue; // avoid double counting

                                    const Atom& ai = atoms[i_idx];
                                    const Atom& aj = atoms[j_idx];
                                    double r = distance(ai, aj);
                                    if (r > cutoff) continue;

                                    double r_min_ij = ai.r_min2 + aj.r_min2;
                                    double eps_ij = std::sqrt(ai.eps * aj.eps);
                                    double lj = eps_ij * (std::pow(r_min_ij / r, 12) - 2.0 * std::pow(r_min_ij / r, 6));
                                    double coulomb = k * (ai.q * aj.q) / (r * d_conv) * qe * qe * e_conv;

                                    // Push both directions for symmetric contribution
                                    results.push_back({ai.resid, aj.resid, lj, coulomb, r});
                                    results.push_back({aj.resid, ai.resid, lj, coulomb, r});
                                }
                            }
                        }
                    }
                }

            }
        }
    }

    return results;
}


#include <Rcpp.h>
using namespace Rcpp;

struct SumCount {
  double lj_sum = 0.0;
  double coul_sum = 0.0;
  double r_sum = 0.0;
  int count = 0;
};

// [[Rcpp::export]]
DataFrame compute_pair_energy_cpp(DataFrame df, double cutoff) {
  NumericVector x = df["x"];
  NumericVector y = df["y"];
  NumericVector z = df["z"];
  NumericVector eps = df["eps"];
  NumericVector r_min2 = df["r_min2"];
  NumericVector q = df["q"];
  CharacterVector resid = df["residue"];

  int n = x.size();
  std::vector<Atom> atoms(n);
  for (int i = 0; i < n; ++i) {
    atoms[i].x = x[i];
    atoms[i].y = y[i];
    atoms[i].z = z[i];
    atoms[i].eps = eps[i];
    atoms[i].r_min2 = r_min2[i];
    atoms[i].q = q[i];
    atoms[i].resid = std::string(resid[i]);
  }

  std::vector<PairEnergy> results = compute_pair_energy_core(atoms, cutoff);

  std::unordered_map<std::string, SumCount> map;
  for (const auto& p : results) {
    std::string key = p.resid_i + "|" + p.resid_j;
    map[key].lj_sum += p.lj_res;
    map[key].coul_sum += p.coulomb_res;
    map[key].r_sum += p.r_avg;
    map[key].count += 1;
  }

  CharacterVector resid_i_map(map.size());
  CharacterVector resid_j_map(map.size());
  NumericVector lj_avg(map.size());
  NumericVector coul_avg(map.size());
  NumericVector r_avg(map.size());

  int idx = 0;
  for (const auto& kv : map) {
    std::string key = kv.first;
    size_t sep = key.find('|');
    resid_i_map[idx] = key.substr(0, sep);
    resid_j_map[idx] = key.substr(sep + 1);
    lj_avg[idx] = kv.second.lj_sum / 2.0; /// kv.second.count;
    coul_avg[idx] = kv.second.coul_sum / 2.0; /// kv.second.count;
    r_avg[idx] = kv.second.r_sum / kv.second.count;
    ++idx;
  }

  return DataFrame::create(
    Named("resid_i") = resid_i_map,
    Named("resid_j") = resid_j_map,
    Named("lj_res") = lj_avg,
    Named("coulomb_res") = coul_avg,
    Named("r_avg") = r_avg
  );
}
