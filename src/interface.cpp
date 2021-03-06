#include "GEE.h"
#include "QIF.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List pgee_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::vec weight,
                  const arma::uvec cluster_sizes, const Rcpp::List family_objs, const std::string corstr,
                  const arma::vec init_beta, const arma::vec init_alpha, double init_phi, bool scale_fix,
                  double lambda, const arma::uvec pindex, double eps, int maxit, double tol,
                  arma::mat cor_mat, int Mv = 0) {
    WorkCor type = workCorMap.at(corstr);
    Family family(family_objs);


    Control ctl(maxit, tol, false);
    Penalty_Options op(lambda, pindex, eps);
    GEE gee(Y, X, offset, weight, cluster_sizes, family, type, ctl, init_beta, init_alpha, init_phi, scale_fix, cor_mat, Mv);
    gee.iterator_penalty(op);

    return gee.get_result();
}

// [[Rcpp::export]]
Rcpp::List gee_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::vec weight,
                 const arma::uvec cluster_sizes, const Rcpp::List family_objs, const std::string corstr,
                 const arma::vec init_beta, const arma::vec init_alpha, double init_phi, bool scale_fix,
                 int maxit, double tol, arma::mat cor_mat, int Mv = 0) {
    WorkCor type = workCorMap.at(corstr);
    Family family(family_objs);

    Control ctl(maxit, tol, false);
    GEE gee(Y, X, offset, weight, cluster_sizes, family, type, ctl, init_beta, init_alpha, init_phi, scale_fix, cor_mat, Mv);
    gee.iterator();

    return gee.get_result();
}

// [[Rcpp::export]]
Rcpp::List qif_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::vec weight,
                 const arma::uvec cluster_sizes, const Rcpp::List family_objs, const std::string corstr,
                 const arma::vec init_beta, int maxit, double tol) {
    WorkCor type = workCorMap.at(corstr);
    Family family(family_objs);


    Control ctl(maxit, tol, false);
    QIF qif(Y, X, offset, weight, cluster_sizes, family, type, ctl, init_beta);
    qif.iterator();

    return qif.get_result();
}

