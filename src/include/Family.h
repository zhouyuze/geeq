#ifndef SRC_FAMILY_H
#define SRC_FAMILY_H

#include "utils.h"

class Family {
private:
    Rcpp::Function link_fun_r;
    Rcpp::Function link_inv_r;
    Rcpp::Function variance_r;
    Rcpp::Function mu_eta_r;
    Rcpp::Function dev_resids_r;
public:
    explicit Family(Rcpp::List family_obj);

    vec link_fun(const vec &mu);

    vec link_inv(const vec &eta);

    vec variance(const vec &mu);

    vec derivative(const vec &eta);

    vec likelyhood(const vec &y, const vec &mu, const vec &wt);
};




#endif //SRC_FAMILY_H
