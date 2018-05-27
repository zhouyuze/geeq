#include "QIF.h"

void QIF::init_base_mat() {
    int size = N / n;
    base_mat.emplace_back(size, size, fill::eye);
    switch (cor_type) {
        case Exchangable:
            base_mat.emplace_back(size, size, fill::ones);
            base_mat[1].diag().fill(0);
            m = 2;
            break;

        case AR1:
            base_mat.emplace_back(size, size, fill::zeros);
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

    Q_first_deriv = vec(p, fill::zeros);
    Q_second_deriv = mat(p, p, fill::zeros);
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

    Q_first_deriv = dev_g.t() * solve(C, g);
    Q_second_deriv = dev_g.t() * solve(C, dev_g);

    vec delta_beta = solve(Q_second_deriv, Q_first_deriv);
    beta = beta - delta_beta;
    update_intermediate_variable();

    return sum(abs(delta_beta));
}

int QIF::iterator() {
    bool stop = false;
    int count = 0;

    while(!stop) {
        count++;
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

Rcpp::List QIF::get_result() {
    return Rcpp::List::create(Rcpp::Named("beta")=beta);
}