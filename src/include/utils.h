#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

class Control {
public:
    int maxit;
    double tol;
    bool trace;
    explicit Control(int maxit = 30, double tol = 10^-6, bool trace = false):
            maxit(maxit), tol(tol), trace(trace) {}
};

enum WorkCor {
    Independence,
    AR1,
    Exchangable
};

enum FamilyType {
    Gaussian,
    Poisson,
    Binomial,
    Gamma,
    InvGaussian
};

enum LinkType {
    Identity,
    Log,
    Logit,
    Reciprocal,
    SquareReciprocal
};

const static map<string, WorkCor> workCorMap = {
        {"independence", Independence},
        {"ar1", AR1},
        {"exchangable", Exchangable}
};

const static map<string, FamilyType> familyTypeMap = {
        {"gaussian", Gaussian},
        {"poisson", Poisson},
        {"binomial", Binomial},
        {"gamma", Gamma},
        {"inverse.gaussian", InvGaussian}
};

const static map<string, LinkType> linkTypeMap = {
        {"identity", Identity},
        {"log", Log},
        {"logit", Logit},
        {"reciprocal", Reciprocal},
        {"1/mu^2", SquareReciprocal}
};

vec identity_link(const vec &mu);
vec identity_link_inv(const vec &eta);
vec identity_deriv(const vec &eta);

vec log_link(const vec &mu);
vec log_link_inv(const vec &eta);
vec log_deriv(const vec &eta);

vec logit_link(const vec &mu);
vec logit_link_inv(const vec &eta);
vec logit_deriv(const vec &eta);

vec reciprocal_link(const vec &mu);
vec reciprocal_link_inv(const vec &eta);
vec reciprocal_deriv(const vec &eta);

vec square_reciprocal_link(const vec &mu);
vec square_reciprocal_link_inv(const vec &eta);
vec square_reciprocal_deriv(const vec &eta);

vec gaussian_variance(const vec &mu);
vec gaussian_likelyhood(const vec &y, const vec &mu);

vec poisson_variance(const vec &mu);
vec poisson_likelyhood(const vec &y, const vec &mu);

vec binomial_variance(const vec &mu);
vec binomial_likelyhood(const vec &y, const vec &mu);

vec gamma_variance(const vec &mu);
vec gamma_likelyhood(const vec &y, const vec &mu);

vec inverse_gaussian_variance(const vec &mu);
vec inverse_gaussian_likelyhood(const vec &y, const vec &mu);

#endif //SRC_UTILS_H
