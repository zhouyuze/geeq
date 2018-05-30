#include "utils.h"

vec identity_link(const vec &mu) {
    return mu;
}

vec identity_link_inv(const vec &eta) {
    return eta;
}

vec identity_deriv(const vec &eta) {
    return vec(eta.n_elem, fill::ones);
}

vec log_link(const vec &mu) {
    return log(mu);
}

vec log_link_inv(const vec &eta) {
    return exp(eta);
}

vec log_deriv(const vec &eta) {
    return exp(eta);
}

vec logit_link(const vec &mu) {
    return log(mu/(1-mu));
}

vec logit_link_inv(const vec &eta) {
    return 1 - 1 / (1 + exp(eta));
}

vec logit_deriv(const vec &eta) {
    return exp(eta) / square(1+exp(eta));
}

vec reciprocal_link(const vec &mu) {
    return 1 / mu;
}

vec reciprocal_link_inv(const vec &eta) {
    return 1 / eta;
}

vec reciprocal_deriv(const vec &eta) {
    return -1 / square(eta);
}

vec square_reciprocal_link(const vec &mu) {
    return 1 / square(mu);
}

vec square_reciprocal_link_inv(const vec &eta) {
    return 1 / sqrt(eta);
}

vec square_reciprocal_deriv(const vec &eta) {
    return -1 / (2 * pow(eta, 1.5));
}

vec gaussian_variance(const vec &mu) {
    return vec(mu.n_elem, fill::ones);
}

vec gaussian_likelyhood(const vec &y, const vec &mu) {
    return y % mu - square(mu) / 2;
}

vec poisson_variance(const vec &mu) {
    return mu;
}

vec poisson_likelyhood(const vec &y, const vec &mu) {
    return y % log(mu) - mu;
}

vec binomial_variance(const vec &mu) {
    return mu % (1 - mu);
}

vec binomial_likelyhood(const vec &y, const vec &mu) {
    return y % log(1 / (1 - mu) - 1) + log(1 - mu);
}

vec gamma_variance(const vec &mu) {
    return square(mu);
}

vec gamma_likelyhood(const vec &y, const vec &mu) {
    return -y / mu - log(mu);
}

vec inverse_gaussian_variance(const vec &mu) {
    return pow(mu, 3);
}

vec inverse_gaussian_likelyhood(const vec &y, const vec &mu) {
    return -y / (2 * square(mu)) + 1 / mu;
}