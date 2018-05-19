#include "geeq.h"

using namespace arma;

RO::RO(vec Beta, vec Alpha, double Phi, bool Converged):
        beta(Beta),
        alpha(Alpha),
        phi(Phi),
        coveraged(Converged)
{}

Family::Family(Rcpp::List family_obj):
        link_funr((SEXP) family_obj["linkfun"]),
        link_invr((SEXP) family_obj["linkinv"]),
        variancer((SEXP) family_obj["variance"]),
        mu_etar((SEXP) family_obj["mu.eta"])
{}

vec Family::link_fun(const vec mu) {
    return Rcpp::as<vec>(link_funr(mu));
}

vec Family::link_inv(const vec eta) {
    return Rcpp::as<vec>(link_invr(eta));
}

vec Family::mu_eta(const vec eta) {
    return Rcpp::as<vec>(mu_etar(eta));
}

vec Family::variance(const vec mu) {
    return Rcpp::as<vec>(variancer(mu));
}