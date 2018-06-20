#include "QIF.h"

QIF::QIF(vec y, mat X, vec offset, uvec cluster_sizes,
        Family family, WorkCor cor_type, Control ctl, vec beta):
        Model(std::move(y), std::move(X), std::move(offset),
              std::move(cluster_sizes), std::move(family),
              cor_type, std::move(ctl), std::move(beta)),
        Q_first_deriv(p, fill::zeros), Q_second_deriv(p, p, fill::zeros) {
    init_base_mat();
    init_key_mat();
}

void QIF::init_base_mat() {
    base_mat.emplace_back(max_cluster, max_cluster, fill::eye);
    switch (cor_type) {
        case Exchangable:
            base_mat.emplace_back(max_cluster, max_cluster, fill::ones);
            base_mat[1].diag().fill(0);
            m = 2;
            break;

        case AR1:
            base_mat.emplace_back(max_cluster, max_cluster, fill::zeros);
            base_mat[1].diag(1).fill(1);
            base_mat[1].diag(-1).fill(1);
            m = 2;
            break;

        default:
            m = 1;
            break;
    }
}


void QIF::init_key_mat() {
    g = vec(m * p, fill::zeros);
    C = mat(m * p, m * p, fill::zeros);
    dev_g = mat(m * p, p, fill::zeros);
}

double QIF::update_beta() {
    g.zeros();
    C.zeros();
    dev_g.zeros();
    int size = N / n;

    mat D = X.each_col() % deriv;
    vec std_err = sqrt(var);
    vec err = y - mu;

    for (int i = 0; i < n; i++) {
        vec gi(arma::size(g), fill::zeros);
        mat dev_gi(arma::size(dev_g), fill::zeros);
        int start = i == 0 ? 0:cluster_bound[i-1];
        int end = cluster_bound[i] - 1;

        vec sub_sqrt_A = std_err.subvec(start, end);
        mat sub_D = D.rows(start, end);

        for (int j = 0; j < base_mat.size(); j++) {
            mat tmp_mat = base_mat[j].submat(0, 0, end - start, end - start);
            mat sub_inverse_var = tmp_mat % (1 / (sub_sqrt_A * sub_sqrt_A.t()));
            gi.subvec(j*p, (j+1)*p-1) = sub_D.t() * sub_inverse_var * (err.subvec(start, end));
            dev_gi.submat(j*p, 0, (j+1)*p-1, p-1) = -sub_D.t() * sub_inverse_var * sub_D;
        }

        g += gi;
        C += gi * gi.t();
        dev_g += dev_gi;
    }

    Q_first_deriv = dev_g.t() * solve(C, g);
    Q_second_deriv = dev_g.t() * solve(C, dev_g);

    vec delta_beta = solve(Q_second_deriv, Q_first_deriv);
    beta = beta - delta_beta;
    update_intermediate_variable();

    return sum(abs(delta_beta));
}

void QIF::iterator() {
    bool stop = false;

    while(!stop) {
        niter++;
        double diff = update_beta();
        if(diff < ctl.tol) {
            converged = true;
            stop = true;
        } else if (niter >= ctl.maxit){
            stop = true;
        }
    }
}

Rcpp::List QIF::get_result() {
    return Rcpp::List::create(Rcpp::Named("beta")=beta);
}