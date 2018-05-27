#include "gee.h"

GEE_Para::GEE_Para(vec y, mat X, vec offset, uvec cluster_sizes,
                   Family family, WorkCor cor_type, GEE_Control ctl,
                   vec beta, vec alpha, double phi):
    y(std::move(y)), X(std::move(X)), offset(std::move(offset)), cluster_bound(cumsum(cluster_sizes)),
    funcs(std::move(family)), cor_type(cor_type), ctl(ctl),
    beta(std::move(beta)), alpha(std::move(alpha)), phi(phi),
    N(this->X.n_rows), p(this->X.n_cols), n(cluster_bound.n_elem) {

    update_intermediate_result();
    this->solved = false;
    this->converged = false;

    for (int i = 0; i < n; i++) {
        cluster_cor.emplace_back(cluster_sizes[i], cluster_sizes[i], fill::eye);
    }
}

int GEE_Para::iterator() {
    bool stop = false;
    int count = 0;

    while(!stop) {
        count++;

        if (!ctl.scale_fix) {
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

int GEE_Para::iterator_penalty(Penalty_Options op) {
    bool stop = false;
    int count = 0;

    while(!stop) {
        count++;
        if (!ctl.scale_fix) {
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

void GEE_Para::update_phi() {
    vec resid = (y - mu) / sqrt(var);
    phi = (N - p) / sum(resid % resid);
}

void GEE_Para::update_alpha() {
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

double GEE_Para::update_beta() {
    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);

    mat hessian = zeros<mat>(p, p);
    vec score = zeros<vec>(p);
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
    update_intermediate_result();

    return sum(abs(delta_beta));
}

double GEE_Para::update_beta_penalty(Penalty_Options op) {
    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);

    mat hessian = zeros<mat>(p, p);
    vec score = zeros<vec>(p);
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
    beta = beta - delta_beta;
    update_intermediate_result();

    return sum(delta_beta);
}

vec GEE_Para::q_scad(double lambda, double a) {
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

void GEE_Para::update_intermediate_result() {
    eta = X * beta;
    mu = funcs.link_inv(eta);
    var = funcs.variance(mu);
    deriv = funcs.mu_eta(eta);
}

vec GEE_Para::get_alpha() {
    return this->alpha;
}

vec GEE_Para::get_beta() {
    return this->beta;
}

double GEE_Para::get_phi() {
    return this->phi;
}

mat GEE_Para::get_sandwich() {
    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec S = y - mu;

    mat hessian = zeros<mat>(p, p);
    mat namesand = zeros<mat>(p, p);

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

double GEE_Para::gaussian_pseudolikelihood() {
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

vec GEE_Para::geodesic_distance() {
    vec delta(2);

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec S = y - mu;

    mat hessian = zeros<mat>(p, p);
    mat namesand = zeros<mat>(p, p);

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