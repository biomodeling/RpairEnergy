#include "compute_pair_energy_core.cpp" // include core C++ logic
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame compute_pair_energy_cpp(DataFrame df, double cutoff) {
    // Estrai colonne
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

    // Chiama la funzione C++ esistente
    std::vector<PairEnergy> results = compute_pair_energy_core(atoms, cutoff);

    // Prepara vettori per R
    CharacterVector resid_i(results.size());
    CharacterVector resid_j(results.size());
    NumericVector lj_res(results.size());
    NumericVector coulomb_res(results.size());
    NumericVector r_avg(results.size());

    for (size_t i = 0; i < results.size(); ++i) {
        resid_i[i] = results[i].resid_i;
        resid_j[i] = results[i].resid_j;
        lj_res[i] = results[i].lj_res;
        coulomb_res[i] = results[i].coulomb_res;
        r_avg[i] = results[i].r_avg;
    }

    return DataFrame::create(
        Named("resid_i") = resid_i,
        Named("resid_j") = resid_j,
        Named("lj_res") = lj_res,
        Named("coulomb_res") = coulomb_res,
        Named("r_avg") = r_avg
    );
}

