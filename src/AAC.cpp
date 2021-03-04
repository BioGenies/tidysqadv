#include <Rcpp.h>
#include <tidysq.h>
#include "tidysqadv/AAC.h"
#include "tidysqadv/Rcpp-export.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::DataFrame CPP_AAC(const Rcpp::List& x,
                        const std::basic_string<char>& NA_letter) {
  return export_to_R(AAC<RCPP_IT>(import_sq_from_R(x, NA_letter)));
}
