#include "QIF_Para.h"

void QIF_Para::update_intermediate_result() {
    eta = X * beta;
    mu = funcs.link_inv(eta);
    var = funcs.variance(mu);
    deriv = funcs.mu_eta(eta);
}

QIF_Para::QIF_Para(vec y, mat X, vec offset, uvec cluster_sizes,
                   Family family, WorkCor cor_type, vec beta):
        y(std::move(y)), X(std::move(X)), offset(std::move(offset)), cluster_bound(cumsum(cluster_sizes)),
        funcs(std::move(family)), cor_type(cor_type), beta(std::move(beta)),
        N(this->X.n_rows), p(this->X.n_cols), n(cluster_bound.n_elem) {
    update_intermediate_result();
    this->solved = false;
    this->converged = false;
    switch (cor_type) {
        case Independence:
            g = vec(p, fill::zeros);
            C = mat(p, p, fill::zeros);
            dev_g = mat(p, p, fill::zeros);
            break;
        default:
            g = vec(2 * p, fill::zeros);
            C = mat(2 * p, 2 * p, fill::zeros);
            dev_g = mat(2 * p, p, fill::zeros);
            break;
    }

    int size = N / n;

    base_mat.emplace_back(size, size, fill::eye);
    switch (cor_type) {
        case Exchangable:
            base_mat.emplace_back(size, size, fill::ones);
            base_mat[1].diag().fill(0);
            break;
        case AR1:
            base_mat.emplace_back(size, size, fill::zeros);
            base_mat[1].diag(1).fill(1);
            base_mat[1].diag(-1).fill(1);
            break;
        default:
            break;
    }
}

double QIF_Para::update_beta() {
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
        int start = i * size;
        int end = (i + 1) * size - 1;
        vec sub_sqrt_A = std_err.subvec(start, end);
        mat sub_D = D.rows(start, end);
        for (int j = 0; j < base_mat.size(); j++) {
            mat sub_inverse_var = base_mat[j] % (1 / (sub_sqrt_A * sub_sqrt_A.t()));
            gi.subvec(j*p, (j+1)*p-1) = sub_D.t() * sub_inverse_var * (err.subvec(start, end));
            dev_gi.submat(j*p, 0, (j+1)*p-1, p-1) = -sub_D.t() * sub_inverse_var * sub_D;
        }
        g += gi;
        C += gi * gi.t();
        dev_g += dev_gi;
    }
    mat qif_dev_first = dev_g.t() * solve(C, g);
    mat qif_dev_second = dev_g.t() * solve(C, dev_g);

    vec delta_beta = solve(qif_dev_second, qif_dev_first);
    beta = beta - delta_beta;
    update_intermediate_result();

    return sum(abs(delta_beta));
}

int QIF_Para::iterator() {
    bool stop = false;
    int count = 0;

    while(!stop) {
        count++;
        double diff = update_beta();
        if(diff < 0.001) {
            converged = true;
            stop = true;
        } else if (count >= 30){
            stop = true;
        }
    }
    return count;
}

vec QIF_Para::get_beta() {
    return this->beta;
}