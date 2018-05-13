#include "gee.h"

using namespace arma;

double update_Phi(const vec &resid, int denominator) {
    return denominator / sum(resid % resid);
}

vec update_WorkCor(mat &Cor, const vec &resid, const uvec cluster_bound,
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

void update_Beta(vec &beta, const mat &X, const vec &err,
                 mat &Cor, const vec &deriv, const vec &var,
                 const uvec &cluster_bound, double phi, bool penalty = false,
                 const uvec &pindex = uvec(), double lambda = 0, double eps = 0) {
    mat D = diagmat(deriv) * X;
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

        mat sub_sqrt_A = diagmat(std_err.subvec(start, end));
        mat sub_sqrt_cor = Cor.submat(start, start, end, end);
        mat sub_inverse_var = (sub_sqrt_A * sub_sqrt_cor * sub_sqrt_A / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        score += sub_D.t() * sub_inverse_var * (err.subvec(start, end));


        start = cluster_bound[i];
    }

    if (penalty) {
        int N = cluster_bound.n_elem;
        delta_beta = (hessian + N * diagmat(E)).i() * (score - N * diagmat(E) * beta);
    } else {
        delta_beta = hessian.i() * score;
    }

    beta = beta + delta_beta;
}

void init_correlation(mat &cor, WorkCor, const vec &init_alpha, const mat &cor_mat) {
    // TODO initialize correlation matrix;
}

RO gee_iteration(const vec &Y, const mat &X, const vec &offset, const uvec &cluster_sizes,
                 Family &funcs, WorkCor type, const mat &cor_mat,
                 const vec &init_beta, const vec &init_alpha, double init_phi, bool scale_fix,
                 bool penalty, double lambda, const uvec &pindex, double eps,
                 int maxit, double tol) {
    uvec cluster_bound = cumsum(cluster_sizes);
    mat cor = eye<mat>(sum(cluster_sizes), sum(cluster_sizes));
    init_correlation(cor, type, init_alpha, cor_mat);

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

        alpha = update_WorkCor(cor, resid, cluster_bound, X.n_cols, phi, type);

        if (penalty) {
            update_Beta(beta, X, Y - mu, cor, deriv, var, cluster_bound, phi,
                        penalty, pindex, lambda, eps);
        } else {
            update_Beta(beta, X, Y - mu, cor, deriv, var, cluster_bound, phi);
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