#include <Rcpp.h>
#include <tidysq.h>
#include "tidysqadv/computing_kernel.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::NumericVector CPP_computing_kernel(const Rcpp::List& x,
                                         const R_xlen_t& max_kmer_length,
                                         const double exponential,
                                         const std::basic_string<char>& NA_letter) {
    Sq<RCPP_IT> sq = import_sq_from_R(x, NA_letter);
    return Rcpp::wrap(correlation_kernel_3(sq[0].get(), sq[1].get(), sq.alphabet().alphabet_size(), max_kmer_length, exponential));
}
