#include "GEE.h"

GEE::GEE(vec y, mat X, vec offset, uvec cluster_sizes,
         Family family, WorkCor cor_type, Control ctl,
         vec beta, vec alpha, double phi, bool fix):
        Model(std::move(y), std::move(X), std::move(offset),
              cluster_sizes, std::move(family),
              cor_type, ctl, std::move(beta)),
        alpha(std::move(alpha)), phi(phi), scale_fix(fix) {
    for (int i = 0; i < n; i++) {
        cluster_cor.emplace_back(cluster_sizes[i], cluster_sizes[i], fill::eye);
    }
}

int GEE::iterator() {
    bool stop = false;
    int count = 0;

    while(!stop) {
        count++;
        if (!scale_fix) {
            update_phi();
        }
        update_alpha();
        double diff = update_beta();

        if(diff < ctl.tol) {
            converged = true;
            stop = true;
        } else if (count >= ctl.maxit){
            stop = true;
        }
    }
    return count;
}

int GEE::iterator_penalty(Penalty_Options op) {
    bool stop = false;
    int count = 0;

    while(!stop) {
        count++;
        if (!scale_fix) {
            update_phi();
        }
        update_alpha();
        double diff = update_beta_penalty(op);

        if(diff < ctl.tol) {
            converged = true;
            stop = true;
        } else if (count >= ctl.maxit){
            stop = true;
        }
    }
    return count;
}

Rcpp::List GEE::get_result() {
    return Rcpp::List::create(Rcpp::Named("beta")=beta,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("phi")=phi,
                              Rcpp::Named("sandwich")=get_sandwich(),
                              Rcpp::Named("gaussian.pseudolikelihood")=gaussian_pseudolikelihood(),
                              Rcpp::Named("geodesic.distance")=geodesic_distance());
}

double GEE::update_beta() {
    mat hessian = zeros<mat>(p, p);
    vec score = zeros<vec>(p);

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec delta_beta = zeros<vec>(p);
    vec err = y - mu;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        mat sub_inverse_var = ((cluster_cor[i] % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        score += sub_D.t() * sub_inverse_var * (err.subvec(start, end));
    }

    delta_beta = solve(hessian, score);
    beta = beta + delta_beta;
    update_intermediate_variable();

    return sum(abs(delta_beta));
}

double GEE::update_beta_penalty(Penalty_Options op) {
    mat hessian = zeros<mat>(p, p);
    vec score = zeros<vec>(p);

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec delta_beta = zeros<vec>(p);
    vec S = y - mu;

    vec E = q_scad(op.lambda) / (abs(beta) + op.eps);
    if (op.pindex.n_elem != 0) {
        E.elem(op.pindex) = zeros<vec>(op.pindex.n_elem);
    }

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        mat sub_inverse_var = ((cluster_cor[i] % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        score += sub_D.t() * sub_inverse_var * (S.subvec(start, end));
    }

    delta_beta = solve(hessian + n * diagmat(E), score - n * (E % beta));
    beta = beta + delta_beta;
    update_intermediate_variable();

    return sum(abs(delta_beta));
}

void GEE::update_phi() {
    vec resid = (y - mu) / sqrt(var);
    phi = (N - p) / sum(resid % resid);
}

void GEE::update_alpha() {
    vec resid = (y - mu) / sqrt(var);

    switch (cor_type) {
        case Independence: {
            break;
        }
        case Exchangable: {
            double numerator = 0;
            double denominator = -p;

            // moment estimation for alpha
            for (int i = 0; i < n; i++) {
                int size = cluster_bound[i] - (i == 0 ? 0:cluster_bound[i-1]);

                vec tmp = resid.subvec(cluster_bound[i] - size, cluster_bound[i] - 1);
                mat R = tmp * tmp.t();
                numerator += accu(R) - sum(diagvec(R));
                denominator += size * (size - 1);
            }
            alpha[0] = numerator / denominator;

            // update the correlation matrix
            for (int i = 0; i < n; i++) {
                cluster_cor[i].fill(alpha[0]);
                cluster_cor[i].diag().fill(1);
            }
            break;
        }
        case AR1: {
            vec next_resid = resid;
            next_resid.head(N - 1) = next_resid.tail(N - 1);
            // set the resid value at cluster bound to zero
            next_resid.elem(cluster_bound - 1) = zeros<vec>(n);

            alpha[0] = sum(resid % next_resid) / (N - n) * phi;

            // update the correlation matrix
            for (int i = 0; i < n; i++) {
                cluster_cor[i].diag(1).fill(alpha[0]);
                cluster_cor[i].diag(-1).fill(alpha[0]);
            }
            break;
        }
    }
}

vec GEE::q_scad(double lambda, double a) {
    vec q = beta;
    q.transform([&](double val) {
        if (val < lambda) {
            val = lambda;
        } else if (val < a * lambda) {
            val = (a * lambda - val) / (a - 1) * lambda;
        }
        return val;
    });
    return q;
}

mat GEE::get_sandwich() {
    mat hessian = zeros<mat>(p, p);
    mat namesand = zeros<mat>(p, p);

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec S = y - mu;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        vec sub_S = S.subvec(start, end);
        mat sub_inverse_var = ((cluster_cor[i] % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        namesand += sub_D.t() * sub_inverse_var * (sub_S * sub_S.t()) * sub_inverse_var * sub_D;
    }

    return solve(hessian, solve(hessian, namesand).t());
}

double GEE::gaussian_pseudolikelihood() {
    vec std_err = sqrt(var);
    vec S = y - mu;
    double result = 0;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0 : cluster_bound[i - 1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        vec sub_S = S.subvec(start, end);
        mat sub_inverse_var = ((cluster_cor[i] % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();
        result += accu(sub_S.t() * sub_inverse_var * sub_S);
        result += det(sub_inverse_var);
    }
    return result / 2;
}

vec GEE::geodesic_distance() {
    vec delta(2);
    mat hessian = zeros<mat>(p, p);
    mat namesand = zeros<mat>(p, p);

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec S = y - mu;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        vec sub_S = S.subvec(start, end);
        mat sub_inverse_var = ((cluster_cor[i] % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();
        mat sub_D = D.rows(start, end);

        hessian += sub_D.t() * sub_inverse_var * sub_D;
        namesand += sub_D.t() * sub_inverse_var * (sub_S * sub_S.t()) * sub_inverse_var * sub_D;
    }

    mat Q = solve(hessian, namesand);
    vec eig = eig_sym(Q);
    delta[0] = sum((eig - 1) % (eig - 1));
    delta[1] = sum(abs(log(eig)));
    return delta;
}