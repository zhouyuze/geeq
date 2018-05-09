#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

enum WorkCor {
    Independent,
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

class RO {
public:
    vec beta;
    double alpha;
    double phi;
    RO(vec Beta, double Alpha, double Phi);
};

RO gee_iteration(const mat X, const vec Y, const uvec cluster_sizes,
                 Family funcs, WorkCor type, int maxit);