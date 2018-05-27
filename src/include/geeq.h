#ifndef SRC_GEEQ_H
#define SRC_GEEQ_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

enum WorkCor {
    Independence,
    AR1,
    Exchangable
};

class Family {
private:
    Rcpp::Function link_funr;
    Rcpp::Function link_invr;
    Rcpp::Function variancer;
    Rcpp::Function mu_etar;
public:
    explicit Family(Rcpp::List family_obj);

    vec link_fun(const vec mu);

    vec link_inv(const vec eta);

    vec variance(const vec mu);

    vec mu_eta(const vec eta);
};

#endif //SRC_GEE_H
