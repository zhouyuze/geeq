#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

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

class RO {
public:
    vec beta;
    vec alpha;
    double phi;
    bool coveraged;
    RO(vec Beta, vec Alpha, double Phi, bool Converged);
};

RO gee_iteration(const vec &Y, const mat &X, const vec &offset, const uvec &cluster_sizes,
                 Family &funcs, WorkCor type, const mat &cor_mat,
                 const vec &init_beta, const vec &init_alpha, double init_phi, bool scale_fix,
                 bool penalty, double lambda, const uvec &pindex, double eps,
                 int maxit, double tol);