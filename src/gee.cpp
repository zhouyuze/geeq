#include "gee.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

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
            // moment estimation for alpha
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;
                int size = end - start + 1;

                vec tmp = r.subvec(start, end);
                vec R = tmp * tmp.t();
                numerator += sum(R) - sum(diagvec(R));
                denominator += size * (size - 1);

                start = cluster_bound[i];
            }
            alpha = numerator / denominator;

            // update the correlation matrix
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;
                Cor.submat(start, start, end, end).fill(alpha);
                start = cluster_bound[i];
            }
            Cor.diag().fill(1);
            break;
        }
        case AR1: {
            vec next_resid = resid;
            next_resid.head(N - 1) = next_resid.tail(N - 1);
            // set the resid value at cluster bound to zero
            next_resid.elem(cluster_bound - 1) = zeros<vec>(n);

            alpha = sum(resid % next_resid) / (N - n) * phi;

            // update the correlation matrix
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
                 mat &Cor, const vec &mu, const vec &deriv,
                 const vec &var, const uvec &cluster_bound, double phi) {
    mat D = diagmat(deriv) * X;
    vec std_err = sqrt(var);

    mat hessian = zeros<mat>(X.n_cols, X.n_cols);
    vec score = zeros<vec>(X.n_cols);
    vec err = Y - mu;

    for (int i = 0, start = 0; i < cluster_bound.n_elem; i++) {
        int end = cluster_bound[i] - 1;

        mat sub_sqrt_A = diagmat(std_err.subvec(start, end));
        mat sub_sqrt_cor = Cor.submat(start, start, end, end);
        mat sub_inverse_var = (sub_sqrt_A * sub_sqrt_cor * sub_sqrt_A / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        score += sub_D.t() * sub_inverse_var * (err.subvec(start, end));

        start = cluster_bound[i];
    }

    vec delta_beta = hessian.i() * score;

    Beta = Beta + delta_beta;
}

RO gee_iteration(const mat X, const vec Y, const uvec cluster_sizes,
                        Family funcs, WorkCor type, int maxit) {

    uvec cluster_bound = cumsum(cluster_sizes);
    mat cor = eye<mat>(sum(cluster_sizes), sum(cluster_sizes));
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
        update_Beta(beta, X, Y, cor, mu, deriv, var, cluster_bound, phi);
    }

    RO result(beta, alpha, phi);
    return result;
}

// [[Rcpp::export]]
Rcpp::List gee(arma::vec y, arma::mat X, arma::uvec clusterSizes, Rcpp::List family_objs) {
    Family family(family_objs);
    RO result = gee_iteration(X, y, clusterSizes, family, AR1, 20);
    return Rcpp::List::create(Rcpp::Named("beta")=result.beta,
                              Rcpp::Named("alpha")=result.alpha,
                              Rcpp::Named("phi")=result.phi);
 }



