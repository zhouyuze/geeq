#include "gee.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List pgee_c(const arma::vec Y, const arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
                  const Rcpp::List family_objs, const std::string corstr,
                  const arma::vec init_beta, const arma::vec init_alpha, double init_phi, bool scale_fix,
                  double lambda, const arma::uvec pindex, double eps, int maxit, double tol) {
    Family family(family_objs);
    WorkCor type;
    if (corstr == "independence") {
        type = Independence;
    } else if (corstr == "ar1") {
        type = AR1;
    } else if (corstr == "exchangeable") {
        type = Exchangable;
    } else {
        return Rcpp::List::create(Rcpp::Named("error") = "Unsupported type");
    }

    RO result = gee_iteration(Y, X, offset, cluster_sizes, family, type, arma::mat(),
                              init_beta, init_alpha, init_phi, scale_fix,
                              true, lambda, pindex, eps, maxit, tol);
    return Rcpp::List::create(Rcpp::Named("beta")=result.beta,
                              Rcpp::Named("alpha")=result.alpha,
                              Rcpp::Named("phi")=result.phi);

}

// [[Rcpp::export]]
Rcpp::List gee_c(const arma::vec Y, const arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
                 const Rcpp::List family_objs, const std::string corstr,
                 const arma::vec init_beta, const arma::vec init_alpha, double init_phi, bool scale_fix,
                 int maxit, double tol) {
    Family family(family_objs);
    WorkCor type;
    if (corstr == "independence") {
        type = Independence;
    } else if (corstr == "ar1") {
        type = AR1;
    } else if (corstr == "exchangeable") {
        type = Exchangable;
    } else {
        return Rcpp::List::create(Rcpp::Named("error") = "Unsupported type");
    }
    RO result = gee_iteration(Y, X, offset, cluster_sizes, family, type, arma::mat(),
                              init_beta, init_alpha, init_phi, scale_fix,
                              false, 0, arma::uvec(), 0, maxit, tol);
    return Rcpp::List::create(Rcpp::Named("beta")=result.beta,
                              Rcpp::Named("alpha")=result.alpha,
                              Rcpp::Named("phi")=result.phi);
}
