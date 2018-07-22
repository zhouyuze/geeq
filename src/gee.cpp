#include "GEE.h"

GEE::GEE(vec y, mat X, vec offset, vec weight, uvec cluster_sizes,
         Family family, WorkCor cor_type, Control ctl,
         vec beta, vec alpha, double phi, bool fix, mat cor_mat, int Mv):
        Model(std::move(y), std::move(X), std::move(offset),
              std::move(weight), std::move(cluster_sizes), std::move(family),
              cor_type, ctl, std::move(beta)),
        alpha(std::move(alpha)), phi(phi), scale_fix(fix), Mv(Mv),
        H1(p, p, fill::zeros), H2(p, p, fill::zeros), score(p, fill::zeros) {

    switch (cor_type) {
        case Fixed:
            correlation = cor_mat;
            break;
        default:
            correlation = mat(max_cluster, max_cluster, fill::eye);
            break;
    }
}

void GEE::iterator() {
    bool stop = false;

    while(!stop) {
        niter++;
        if (!scale_fix) {
            update_phi();
        }
        update_alpha();
        double diff = update_beta();
        if(diff < ctl.tol) {
            converged = true;
            stop = true;
        } else if (niter >= ctl.maxit){
            stop = true;
        }
    }
    // Calculate H2 after estimation
    calculate_H2();
}

void GEE::iterator_penalty(Penalty_Options op) {
    bool stop = false;

    while(!stop) {
        niter++;
        if (!scale_fix) {
            update_phi();
        }
        update_alpha();
        double diff = update_beta_penalty(op);

        if(diff < ctl.tol) {
            converged = true;
            stop = true;
        } else if (niter >= ctl.maxit){
            stop = true;
        }
    }
    // Calculate H2 after estimation
    calculate_H2();
}

Rcpp::List GEE::get_result() {
    return Rcpp::List::create(Rcpp::Named("beta")=beta,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("phi")=phi,
                              Rcpp::Named("sandwich")=get_sandwich(),
                              Rcpp::Named("niter")=niter,
                              Rcpp::Named("gaussian.pseudolikelihood")=gaussian_pseudolikelihood(),
                              Rcpp::Named("geodesic.distance")=geodesic_distance(),
                              Rcpp::Named("QIC")=QIC());
}

double GEE::update_beta() {
    H1.zeros();
    score.zeros();

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec delta_beta = zeros<vec>(p);
    vec err = y - mu;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        mat cluster_cor = correlation.submat(0, 0, end - start, end - start);
        mat sub_inverse_var = ((cluster_cor % (sub_sqrt_A * sub_sqrt_A.t())) * phi).i();
        mat sub_D = D.rows(start, end);

        mat W = diagmat(weight.subvec(start, end));

        H1 += sub_D.t() * sub_inverse_var * W * sub_D;
        score += sub_D.t() * sub_inverse_var * W * (err.subvec(start, end));
    }

    delta_beta = solve(H1, score);
    beta = beta + delta_beta;
    update_intermediate_variable();

    return sum(abs(delta_beta));
}

double GEE::update_beta_penalty(Penalty_Options op) {
    H1.zeros();
    score.zeros();

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
        mat cluster_cor = correlation.submat(0, 0, end - start, end - start);
        mat sub_inverse_var = ((cluster_cor % (sub_sqrt_A * sub_sqrt_A.t())) * phi).i();
        mat sub_D = D.rows(start, end);

        mat W = diagmat(weight.subvec(start, end));

        H1 += sub_D.t() * sub_inverse_var * W * sub_D;
        score += sub_D.t() * sub_inverse_var * W * (S.subvec(start, end));
    }

    delta_beta = solve(H1 + n * diagmat(E), score - n * (E % beta));
    beta = beta + delta_beta;
    update_intermediate_variable();

    return sum(abs(delta_beta));
}

void GEE::update_phi() {
    vec resid = (y - mu) / sqrt(var);
    phi = sum(resid % resid % weight) / (N - p);
}

void GEE::update_alpha() {
    vec resid = (y - mu) % sqrt(weight) / sqrt(var);

    switch (cor_type) {
        case Independence:
        case Fixed: {
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
            denominator -= p;
            alpha[0] = numerator / (denominator * phi);

            // update the correlation matrix
            correlation.fill(alpha[0]);
            correlation.diag().fill(1);
            break;
        }
        case AR1: {
            vec next_resid = resid;
            next_resid.head(N - 1) = next_resid.tail(N - 1);
            // set the resid value at cluster bound to zero
            next_resid.elem(cluster_bound - 1) = zeros<vec>(n);

            alpha[0] = sum(resid % next_resid) / ((N - n - p) * phi);

            // update the correlation matrix
            correlation.diag(1).fill(alpha[0]);
            correlation.diag(-1).fill(alpha[0]);
            break;
        }
        case Unstructured: {
            mat tmp(max_cluster, max_cluster, fill::zeros);
            // number of clusters for different size
            uvec count(max_cluster, fill::zeros);
            for (int i = 0; i < n; i++) {
                int start = i == 0 ? 0:cluster_bound[i-1];
                int end = cluster_bound[i] - 1;

                count[end - start]++;

                vec sub_resid = resid.subvec(start, end);
                tmp.submat(0, 0, end-start, end-start) += sub_resid * sub_resid.t();
            }

            for (int index = 0, i = 0; i < max_cluster - 1; i++) {
                for (int j = i + 1; j < max_cluster; j++) {
                    double numerator = tmp(i, j);
                    double denominator = sum(count) - sum(count.head(j)) - p;
                    alpha[index] = numerator / (denominator * phi);
                    tmp(i, j) = alpha[index];
                    tmp(j, i) = alpha[index++];
                }
            }
            tmp.diag().fill(1);

            // update the correlation matrix
            correlation = tmp;
            break;
        }
        case M_dependent: {
            vec numerator(Mv, fill::zeros);
            vec denominator(Mv, fill::zeros);

            for (int i = 0; i < n; i++) {
                int start = i == 0 ? 0:cluster_bound[i-1];
                int end = cluster_bound[i] - 1;

                vec sub_resid = resid.subvec(start, end);
                mat tmp = sub_resid * sub_resid.t();
                for (int j = end - start, k = 0; j > 0 && k < Mv; j++, k++) {
                    numerator[k] += sum(tmp.diag(k+1));
                    denominator[k] += j;
                }
            }
            denominator -= p;
            alpha = numerator / (denominator * phi);

            // update the correlation matrix
            mat tmp(max_cluster, max_cluster, fill::zeros);
            for (int i = 0; i < Mv; i++) {
                tmp.diag(i + 1).fill(alpha[i]);
                tmp.diag(-i - 1).fill(alpha[i]);
            }
            tmp.diag().fill(1);
            correlation = tmp;
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
    return solve(H1, solve(H1, H2).t());
}

double GEE::QIC() {
    mat DAD(p, p, fill::zeros);
    mat D = X.each_col() % deriv;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0 : cluster_bound[i - 1];
        int end = cluster_bound[i] - 1;
        vec sub_var = var.subvec(start, end);
        mat sub_D = D.rows(start, end);
        DAD += sub_D.t() * diagmat(1 / sub_var) * sub_D;
    }

    return -2 * sum(funcs.likelyhood(y, mu, weight)) / phi
           + 2 * trace(DAD * get_sandwich());
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
        mat cluster_cor = correlation.submat(0, 0, end - start, end - start);
        mat sub_inverse_var = ((cluster_cor % (sub_sqrt_A * sub_sqrt_A.t())) / phi).i();

        mat W = diagmat(weight.subvec(start, end));

        result += accu(sub_S.t() * sub_inverse_var * W * sub_S);
        result += det(sub_inverse_var);
    }
    return result / 2;
}

vec GEE::geodesic_distance() {
    vec delta(2);
    mat Q = solve(H1, H2);
    vec eig = eig_sym(Q);
    delta[0] = sum((eig - 1) % (eig - 1));
    delta[1] = sum(abs(log(eig)));
    return delta;
}

void GEE::calculate_H2() {
    H2.zeros();

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec S = y - mu;

    for (int i = 0; i < n; i++) {
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        vec sub_S = S.subvec(start, end);
        mat cluster_cor = correlation.submat(0, 0, end - start, end - start);
        mat sub_inverse_var = ((cluster_cor % (sub_sqrt_A * sub_sqrt_A.t())) * phi).i();
        mat sub_D = D.rows(start, end);
        mat W = diagmat(weight.subvec(start, end));

        H2 += sub_D.t() * sub_inverse_var * W * (sub_S * sub_S.t()) * W * sub_inverse_var * sub_D;
    }
}