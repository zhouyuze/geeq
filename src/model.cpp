#include "model.h"

void Model::update_intermediate_variable() {
    eta = X * beta;
    mu = funcs.link_inv(eta);
    var = funcs.variance(mu);
    deriv = funcs.derivative(eta);
}

Model::Model(vec y, mat X, vec offset, uvec cluster_sizes,
             Family family, WorkCor cor_type, Control ctl, vec init_beta):
        y(std::move(y)), X(std::move(X)), offset(std::move(offset)), cluster_bound(cumsum(cluster_sizes)),
        funcs(std::move(family)), cor_type(cor_type), ctl(ctl), beta(std::move(init_beta)),
        N(this->X.n_rows), p(this->X.n_cols), n(cluster_bound.n_elem), max_cluster(max(cluster_sizes)), converged(false) {
    update_intermediate_variable();
}

