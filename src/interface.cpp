#include "GEE.h"
#include "QIF.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List pgee_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
                  const Rcpp::List family_objs, const std::string corstr,
                  const arma::vec init_beta, const arma::vec init_alpha, double init_phi, bool scale_fix,
                  double lambda, const arma::uvec pindex, double eps, int maxit, double tol) {
    string familystr = Rcpp::as<string>(family_objs["family"]);
    string linkstr = Rcpp::as<string>(family_objs["link"]);
    if (workCorMap.count(corstr) == 0 ||
        familyTypeMap.count(familystr) == 0 ||
        linkTypeMap.count(linkstr) == 0) {
        return Rcpp::List::create(Rcpp::Named("error") = "Unsupported type");
    }

    WorkCor type = workCorMap.at(corstr);
    Family family(familyTypeMap.at(familystr), linkTypeMap.at(linkstr));


    Control ctl(maxit, tol, false);
    Penalty_Options op(lambda, pindex, eps);
    GEE gee(Y, X, offset, cluster_sizes, family, type, ctl, init_beta, init_alpha, init_phi, scale_fix);
    gee.iterator_penalty(op);

    return gee.get_result();
}

// [[Rcpp::export]]
Rcpp::List gee_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
                 const Rcpp::List family_objs, const std::string corstr,
                 const arma::vec init_beta, const arma::vec init_alpha, double init_phi, bool scale_fix,
                 int maxit, double tol) {
    string familystr = Rcpp::as<string>(family_objs["family"]);
    string linkstr = Rcpp::as<string>(family_objs["link"]);
    if (workCorMap.count(corstr) == 0 ||
            familyTypeMap.count(familystr) == 0 ||
            linkTypeMap.count(linkstr) == 0) {
        return Rcpp::List::create(Rcpp::Named("error") = "Unsupported type");
    }

    WorkCor type = workCorMap.at(corstr);
    Family family(familyTypeMap.at(familystr), linkTypeMap.at(linkstr));

    Control ctl(maxit, tol, false);
    GEE gee(Y, X, offset, cluster_sizes, family, type, ctl, init_beta, init_alpha, init_phi, scale_fix);
    gee.iterator();

    return gee.get_result();
}

// [[Rcpp::export]]
Rcpp::List qif_c(const arma::vec Y, arma::mat X, const arma::vec offset, const arma::uvec cluster_sizes,
                 const Rcpp::List family_objs, const std::string corstr,
                 const arma::vec init_beta, int maxit, double tol) {
    string familystr = Rcpp::as<string>(family_objs["family"]);
    string linkstr = Rcpp::as<string>(family_objs["link"]);
    if (workCorMap.count(corstr) == 0 ||
        familyTypeMap.count(familystr) == 0 ||
        linkTypeMap.count(linkstr) == 0) {
        return Rcpp::List::create(Rcpp::Named("error") = "Unsupported type");
    }

    WorkCor type = workCorMap.at(corstr);
    Family family(familyTypeMap.at(familystr), linkTypeMap.at(linkstr));


    Control ctl(maxit, tol, false);
    QIF qif(Y, X, offset, cluster_sizes, family, type, ctl, init_beta);
    qif.iterator();

    return qif.get_result();
}
