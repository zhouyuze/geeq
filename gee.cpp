#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]
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
    explicit Family(Rcpp::List family_obj):
            link_funr((SEXP)family_obj["linkfun"]),
            link_invr((SEXP)family_obj["linkinv"]),
            variancer((SEXP)family_obj["variance"]),
            mu_etar((SEXP)family_obj["mu.eta"])
    {}
    vec link_fun(const vec mu) {
        return Rcpp::as<vec>(link_funr(mu));
    }
    vec link_inv(const vec eta) {
        return Rcpp::as<vec>(link_invr(eta));
    }
    vec variance(const vec mu) {
        return Rcpp::as<vec>(variancer(mu));
    }
    vec mu_eta(const vec eta) {
        return Rcpp::as<vec>(mu_etar(eta));
    }
};

class RO {
public:
    vec beta;
    double alpha;
    double phi;
    RO(vec Beta, double Alpha, double Phi);
};

RO::RO(vec Beta, double Alpha, double Phi) {
    beta = Beta;
    alpha = Alpha;
    phi = Phi;
}

double update_Phi(const vec &resid, int denominator) {
    return denominator / sum(resid % resid);
}

double update_WorkCor(mat &Cor, const vec &resid, const uvec cluster_bound,
                    int p, double phi, WorkCor type) {
    int N = resid.n_elem;
    int n = cluster_bound.n_elem;
    double alpha;
    switch (type) {
        case Independent: {
            break;
        }
        case Exchangable: {
            double numerator = 0;
            double denominator = -p;
            vec r = resid;
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;
                int scope = end - start + 1;

                vec tmp = r.subvec(start, end);
                vec R = tmp * tmp.t();
                numerator += sum(R) - sum(diagvec(R));
                denominator += scope * (scope - 1);

                start = end + 1;
            }

            alpha = numerator / denominator;
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;

                Cor.submat(start, start, end, end).fill(alpha);

                start = end + 1;
            }
            Cor.diag().fill(1);
            break;
        }
        case AR1: {
            vec next_resid = resid;
            next_resid.head(N - 1) = next_resid.tail(N - 1);
            next_resid.elem(cluster_bound - 1) = zeros<vec>(n);

            alpha = sum(resid % next_resid) / (N - n) * phi;

            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;
                Cor.submat(start, start, end, end).diag(1).fill(alpha);
                Cor.submat(start, start, end, end).diag(-1).fill(alpha);
            }
            break;
        }
    }
    return alpha;
}

void update_Beta(vec &Beta, const mat &X, const vec &Y,
                 const mat &Cor, const vec &mu, const vec &deriv,
                 const vec &var, double phi) {
    mat D = diagmat(deriv) * X;
    mat std_err = diagmat(sqrt(var));
    mat invV = (std_err * Cor * std_err / phi).i();

    mat hessian = D.t() * invV * D;
    mat score = D.t() * invV * (Y - mu);
    vec delta_beta = hessian.i() * score;

    Beta = Beta + delta_beta;
}

RO gee_iteration(const mat X, const vec Y, const uvec clusterSizes,
                        Family funcs, WorkCor type, int maxit) {

    uvec cluster_bound = cumsum(clusterSizes);
    mat cor = eye<mat>(sum(clusterSizes), sum(clusterSizes));
    vec beta(X.n_cols, fill::zeros);
    beta[0] = mean(Y);

    bool stop = false;
    bool converged = false;
    int count = 0;

    double phi = 0, alpha = 0;
    while(!stop) {
        count++;
        if (maxit <= count) {
            stop = true;
        }

        vec Beta_old = beta;
        vec eta = X * beta;
        vec mu = funcs.link_inv(eta);
        vec var = funcs.variance(mu);
        vec deriv = funcs.mu_eta(eta);
        vec resid = (Y - mu) / sqrt(var);

        phi = update_Phi(resid, X.n_rows - X.n_cols);
        alpha = update_WorkCor(cor, resid, cluster_bound, X.n_cols, phi, type);
        update_Beta(beta, X, Y, cor, mu, deriv, var, phi);
    }

    RO result(beta, alpha, phi);
    return result;
}

// [[Rcpp::export]]
Rcpp::List gee(SEXP ys, SEXP Xs, SEXP clusterSizes, SEXP family_objs) {
    mat X = Rcpp::as<mat>(Xs);
    vec y = Rcpp::as<vec>(ys);
    uvec cs = Rcpp::as<uvec>(clusterSizes);
    Family family(Rcpp::as<Rcpp::List>(family_objs));

    RO result = gee_iteration(X, y, cs, family, AR1, 20);
    return Rcpp::List::create(Rcpp::Named("beta")=result.beta,
                              Rcpp::Named("alpha")=result.alpha,
                              Rcpp::Named("phi")=result.phi);
}

