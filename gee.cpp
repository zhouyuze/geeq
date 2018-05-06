#include <RcppArmadillo.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]

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
    arma::vec link_fun(const arma::vec mu) {
        return Rcpp::as<arma::vec>(link_funr(mu));
    }
    arma::vec link_inv(const arma::vec eta) {
        return Rcpp::as<arma::vec>(link_invr(eta));
    }
    arma::vec variance(const arma::vec mu) {
        return Rcpp::as<arma::vec>(variancer(mu));
    }
    arma::vec mu_eta(const arma::vec eta) {
        return Rcpp::as<arma::vec>(mu_etar(eta));
    }
};

class RO {
public:
    arma::vec beta;
    double alpha;
    double phi;
    RO(arma::vec Beta, double Alpha, double Phi);
};

RO::RO(arma::vec Beta, double Alpha, double Phi) {
    beta = Beta;
    alpha = Alpha;
    phi = Phi;
}

double update_Phi(const arma::vec &resid, int denominator) {
    return denominator / sum(resid % resid);
}

double update_WorkCor(arma::mat &Cor, const arma::vec &resid, const arma::uvec cluster_bound,
                    int p, double phi, WorkCor type) {
    int N = resid.n_elem;
    int n = cluster_bound.n_elem;
    double alpha;
    switch (type) {
        case Independent: {
            break;
        }
        case Exchangable: {
            double Sum = 0;
            double denominator = -p;
            arma::vec r = resid;
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;
                int scope = end - start + 1;

                arma::vec tmp = r.subvec(start, end);
                arma::vec R = tmp * tmp.t();
                Sum += sum(R) - sum(diagvec(R));
                denominator += scope * (scope - 1);

                start = end + 1;
            }

            alpha = Sum / denominator;
            for (int i = 0, start = 0; i < n; i++) {
                int end = cluster_bound[i] - 1;

                Cor.submat(start, start, end, end).fill(alpha);

                start = end + 1;
            }
            Cor.diag().fill(1);
            break;
        }
        case AR1: {
            arma::vec next_resid = resid;
            next_resid.head(N - 1) = next_resid.tail(N - 1);
            next_resid.elem(cluster_bound - 1) = arma::zeros<arma::vec>(n);

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

void update_Beta(arma::vec &Beta, const arma::mat &X, const arma::vec &Y,
                 const arma::mat &Cor, const arma::vec &mu, const arma::vec &deriv,
                 const arma::vec &var, double phi) {
    arma::mat D = diagmat(deriv) * X;
    arma::mat std_err = diagmat(sqrt(var));
    arma::mat invV = (std_err * Cor * std_err / phi).i();

    arma::mat hessian = D.t() * invV * D;
    arma::mat score = D.t() * invV * (Y - mu);
    arma::vec delta_beta = hessian.i() * score;

    Beta = Beta + delta_beta;
}

RO gee_iteration(const arma::mat X, const arma::vec Y, const arma::uvec clusterSizes,
                        Family funcs, WorkCor type, int maxit) {

    arma::uvec cluster_bound = arma::cumsum(clusterSizes);
    arma::mat cor = arma::eye<arma::mat>(sum(clusterSizes), sum(clusterSizes));
    arma::vec beta(X.n_cols, arma::fill::zeros);
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

        arma::vec Beta_old = beta;
        arma::vec eta = X * beta;
        arma::vec mu = funcs.link_inv(eta);
        arma::vec var = funcs.variance(mu);
        arma::vec deriv = funcs.mu_eta(eta);
        arma::vec resid = (Y - mu) / arma::sqrt(var);

        phi = update_Phi(resid, X.n_rows - X.n_cols);
        alpha = update_WorkCor(cor, resid, cluster_bound, X.n_cols, phi, type);
        update_Beta(beta, X, Y, cor, mu, deriv, var, phi);
    }

    RO result(beta, alpha, phi);
    return result;
}

// [[Rcpp::export]]
Rcpp::List gee(SEXP ys, SEXP Xs, SEXP clusterSizes, SEXP family_objs) {
    arma::mat X = Rcpp::as<arma::mat>(Xs);
    arma::vec y = Rcpp::as<arma::vec>(ys);
    arma::uvec cs = Rcpp::as<arma::uvec>(clusterSizes);
    Family family(Rcpp::as<Rcpp::List>(family_objs));

    RO result = gee_iteration(X, y, cs, family, AR1, 20);
    return Rcpp::List::create(Rcpp::Named("beta")=result.beta,
                              Rcpp::Named("alpha")=result.alpha,
                              Rcpp::Named("phi")=result.phi);
}

