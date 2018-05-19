#include "gee.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List pgee_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
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

    GEE_Control ctl(scale_fix, maxit, tol, false);
    Penalty_Options op(lambda, pindex, eps);
    GEE_Para gee_para(Y, X, offset, cluster_sizes, family, type, ctl, init_beta, init_alpha, init_phi);
    gee_para.iterator_penalty(op);

    return Rcpp::List::create(Rcpp::Named("beta")=gee_para.get_beta(),
                              Rcpp::Named("alpha")=gee_para.get_alpha(),
                              Rcpp::Named("phi")=gee_para.get_phi(),
                              Rcpp::Named("sandwich")=gee_para.get_sandwich(),
                              Rcpp::Named("gaussian pseudolikelihood")=gee_para.gaussian_pseudolikelihood(),
                              Rcpp::Named("geodesic distance")=gee_para.geodesic_distance());

}

// [[Rcpp::export]]
Rcpp::List gee_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
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

    GEE_Control ctl(scale_fix, maxit, tol, false);
    GEE_Para gee_para(Y, X, offset, cluster_sizes, family, type, ctl, init_beta, init_alpha, init_phi);
    gee_para.iterator();

    return Rcpp::List::create(Rcpp::Named("beta")=gee_para.get_beta(),
                              Rcpp::Named("alpha")=gee_para.get_alpha(),
                              Rcpp::Named("phi")=gee_para.get_phi(),
                              Rcpp::Named("sandwich")=gee_para.get_sandwich(),
                              Rcpp::Named("gaussian.pseudolikelihood")=gee_para.gaussian_pseudolikelihood(),
                              Rcpp::Named("geodesic.distance")=gee_para.geodesic_distance());
}
