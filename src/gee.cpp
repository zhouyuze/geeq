#include "gee.h"

using namespace arma;
using namespace std;

double update_Phi(const vec &resid, int denominator) {
    return denominator / sum(resid % resid);
}

vec update_WorkCor(vector<mat> &cluster_cor, const vec &resid, const uvec cluster_bound,
                    int p, double phi, WorkCor type) {
    int N = resid.n_elem;
    int n = cluster_bound.n_elem;
    double alpha;
    switch (type) {
        case Independence: {
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
                mat R = tmp * tmp.t();
                numerator += accu(R) - sum(diagvec(R));
                denominator += size * (size - 1);

                start = cluster_bound[i];
            }
            alpha = numerator / denominator;

            // update the correlation matrix
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;
                cluster_cor[i].fill(alpha);
                cluster_cor[i].diag().fill(1);
                start = cluster_bound[i];
            }
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
                cluster_cor[i].diag(1).fill(alpha);
                cluster_cor[i].diag(-1).fill(alpha);
            }
            break;
        }
    }
    vec result = {alpha};
    return result;
}

vec q_scad(vec beta, double lambda, double a=3.7) {
    beta.transform([&](double val) {
        if (val < lambda) {
            val = lambda;
        } else if (val < a * lambda) {
            val = (a * lambda - val) / (a - 1) * lambda;
        }
        return val;
    });
    return beta;
}

void update_Beta(vec &beta, mat &X, const vec &err,
                 vector<mat> &cluster_cor, const vec &deriv, const vec &var,
                 const uvec &cluster_bound, double phi, bool penalty = false,
                 const uvec &pindex = uvec(), double lambda = 0, double eps = 0) {
    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);

    mat hessian = zeros<mat>(X.n_cols, X.n_cols);
    vec score = zeros<vec>(X.n_cols);
    vec delta_beta = zeros<vec>(X.n_cols);
    vec E = vec();
    if (penalty) {
        E = q_scad(beta, lambda) / (abs(beta) + eps);
        if (pindex.n_elem != 0) {
            E.elem(pindex) = zeros<vec>(pindex.n_cols);
        }
    }

    for (int i = 0, start = 0; i < cluster_bound.n_elem; i++) {
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        mat sub_inverse_var = ((cluster_cor[i] % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        score += sub_D.t() * sub_inverse_var * (err.subvec(start, end));

        start = cluster_bound[i];
    }

    if (penalty) {
        int N = cluster_bound.n_elem;
        delta_beta = solve(hessian + N * diagmat(E), score - N * (E % beta));
    } else {
        delta_beta = solve(hessian, score);
    }

    beta = beta + delta_beta;
}

void init_correlation(vector<mat> &cluster_cors, const uvec &cluster_sizes, WorkCor type,
                      const vec &init_alpha, const mat &cor_mat) {
    for (int i = 0; i < cluster_sizes.n_elem; i++) {
        cluster_cors.emplace_back(cluster_sizes[i], cluster_sizes[i], fill::eye);
    }
}

RO gee_iteration(const vec &Y, mat &X, const vec &offset, const uvec &cluster_sizes,
                 Family &funcs, WorkCor type, const mat &cor_mat,
                 const vec &init_beta, const vec &init_alpha, double init_phi, bool scale_fix,
                 bool penalty, double lambda, const uvec &pindex, double eps,
                 int maxit, double tol) {
    uvec cluster_bound = cumsum(cluster_sizes);
    vector<mat> cluster_cors;
    init_correlation(cluster_cors, cluster_sizes, type, init_alpha, cor_mat);

    vec beta = init_beta;
    vec alpha = init_alpha;
    double phi = init_phi;

    bool stop = false;
    bool converged = false;
    int count = 0;

    while(!stop) {
        count++;
        vec beta_old = beta;
        vec eta = X * beta;
        vec mu = funcs.link_inv(eta);
        vec var = funcs.variance(mu);
        vec deriv = funcs.mu_eta(eta);
        vec resid = (Y - mu) / sqrt(var);

        if (!scale_fix) {
            phi = update_Phi(resid, X.n_rows - X.n_cols);
        }

        alpha = update_WorkCor(cluster_cors, resid, cluster_bound, X.n_cols, phi, type);

        if (penalty) {
            update_Beta(beta, X, Y - mu, cluster_cors, deriv, var, cluster_bound, phi,
                        penalty, pindex, lambda, eps);
        } else {
            update_Beta(beta, X, Y - mu, cluster_cors, deriv, var, cluster_bound, phi);
        }

        if(sum(abs(beta_old - beta)) < tol) {
            converged = true;
            stop = true;
        } else if (count >= maxit){
            stop = true;
        }
    }

    RO result(beta, alpha, phi, converged);
    return result;
}