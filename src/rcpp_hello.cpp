#include <Rcpp.h>
#include <tidysq.h>

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List rcpp_hello() {
  Alphabet x(DNA_BSC);
  Rcpp::NumericVector y = Rcpp::NumericVector::create(0.0, 1.0);
  return Rcpp::List::create(export_to_R(x), y);
}
