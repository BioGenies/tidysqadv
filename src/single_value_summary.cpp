#include <Rcpp.h>
#include <tidysq.h>
#include "tidysqadv/single_value_summary.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::NumericVector CPP_molecular_weight(const Rcpp::List& x,
                                         const std::string& type,
                                         const std::string& NA_letter) {
    return Rcpp::wrap(molecular_weight<RCPP_IT>(import_sq_from_R(x, NA_letter), type));
}
