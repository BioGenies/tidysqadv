#pragma once

#include <Rcpp.h>
#include <tidysq.h>
#include "tidysqadv/internal/AACFrame.h"

namespace tidysq {
    inline Rcpp::DataFrame export_to_R(const internal::AACFrame<RCPP_IT> &aac) {
        Rcpp::DataFrame ret = Rcpp::DataFrame::create(
                Rcpp::Named("A", aac.A_frequency()),
                Rcpp::Named("C", aac.C_frequency()),
                Rcpp::Named("D", aac.D_frequency()),
                Rcpp::Named("E", aac.E_frequency()),
                Rcpp::Named("F", aac.F_frequency()),
                Rcpp::Named("G", aac.G_frequency()),
                Rcpp::Named("H", aac.H_frequency()),
                Rcpp::Named("I", aac.I_frequency()),
                Rcpp::Named("K", aac.K_frequency()),
                Rcpp::Named("L", aac.L_frequency()),
                Rcpp::Named("M", aac.M_frequency()),
                Rcpp::Named("N", aac.N_frequency()),
                Rcpp::Named("P", aac.P_frequency()),
                Rcpp::Named("Q", aac.Q_frequency()),
                Rcpp::Named("R", aac.R_frequency()),
                Rcpp::Named("S", aac.S_frequency()),
                Rcpp::Named("T", aac.T_frequency()),
                Rcpp::Named("V", aac.V_frequency()),
                Rcpp::Named("W", aac.W_frequency()),
                Rcpp::Named("Y", aac.Y_frequency()));
        ret.attr("class") = Rcpp::StringVector{"tbl_df", "tbl", "data.frame"};
        return ret;
    }
}
