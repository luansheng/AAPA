// topk.cpp — Rcpp acceleration for top-k selection
// Placeholder for Phase 2 (Rcpp optimization)
//
// This file will contain an optimized C++ implementation for top-k
// family selection using partial sorting, avoiding full sort overhead.
//
// Target function signature:
//   List topk_select_cpp(
//     NumericMatrix score_mat,    // N_test x F
//     NumericMatrix conflict_mat, // N_test x F
//     int k
//   )
//
// Returns: List of per-individual top-k family indices and scores
//
// TODO (Phase 2):
// - Use std::partial_sort or std::nth_element
// - Consider early termination with conflict pre-filtering
