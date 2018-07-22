#include "Family.h"

Family::Family(Rcpp::List family_obj):
        link_fun_r((SEXP) family_obj["linkfun"]),
        link_inv_r((SEXP) family_obj["linkinv"]),
        variance_r((SEXP) family_obj["variance"]),
        mu_eta_r((SEXP) family_obj["mu.eta"]),
        dev_resids_r((SEXP) family_obj["dev.resids"])
{}

vec Family::link_fun(const vec &mu) {
    return Rcpp::as<vec>(link_fun_r(mu));
}

vec Family::link_inv(const vec &eta) {
    return Rcpp::as<vec>(link_inv_r(eta));
}

vec Family::variance(const vec &mu) {
    return Rcpp::as<vec>(variance_r(mu));
}

vec Family::derivative(const vec &eta) {
    return Rcpp::as<vec>(mu_eta_r(eta));
}

vec Family::likelyhood(const vec &y, const vec &mu, const vec &wt) {
    vec dev_resid = Rcpp::as<vec>(dev_resids_r(y, mu, wt));
    return  - dev_resid / 2;
}
