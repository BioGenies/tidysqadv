#include <Rcpp.h>
#include <tidysq.h>
#include "tidysqadv/computing_kernel.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::NumericMatrix CPP_computing_kernel(const Rcpp::List& x,
                                         const R_xlen_t& max_kmer_length,
                                         const double exponential,
                                         const std::basic_string<char>& NA_letter) {
    Sq<RCPP_IT> sq = import_sq_from_R(x, NA_letter);
    std::vector<std::vector<double>> result = correlation_kernel_3(sq, max_kmer_length, exponential);
    std::vector<double> result_flattened = std::accumulate(result.begin(), result.end(), std::vector<double>(),
            [](std::vector<double> column_1, std::vector<double> column_2) {
        column_1.insert(column_1.end(), column_2.begin(), column_2.end());
        return column_1;
    });
    return Rcpp::NumericMatrix(sq.size(), sq.size(), result_flattened.cbegin());
}
