#include "Family.h"

Family::Family(FamilyType family_type):
        family_type(family_type) {
    switch (family_type) {
        case Gaussian:
            link_type = Identity;
            link_fun = identity_link;
            link_inv = identity_link_inv;
            derivative = identity_deriv;
            variance = gaussian_variance;
            likelyhood = gaussian_likelyhood;
            break;
        case Poisson:
            link_type = Log;
            link_fun = log_link;
            link_inv = log_link_inv;
            derivative = log_deriv;
            variance = poisson_variance;
            likelyhood = poisson_likelyhood;
            break;
        case Binomial:
            link_type = Logit;
            link_fun = logit_link;
            link_inv = logit_link_inv;
            derivative = logit_deriv;
            variance = binomial_variance;
            likelyhood = binomial_likelyhood;
            break;
        case Gamma:
            link_type = Reciprocal;
            link_fun = reciprocal_link;
            link_inv = reciprocal_link_inv;
            derivative = reciprocal_deriv;
            variance = gamma_variance;
            likelyhood = gamma_likelyhood;
            break;
        case InvGaussian:
            link_type = SquareReciprocal;
            link_fun = square_reciprocal_link;
            link_inv = square_reciprocal_link_inv;
            derivative = square_reciprocal_deriv;
            variance = inverse_gaussian_variance;
            likelyhood = inverse_gaussian_likelyhood;
            break;
    }
}

Family::Family(FamilyType family_type, LinkType link_type):
        family_type(family_type),
        link_type(link_type) {
    switch (family_type) {
        case Gaussian:
            variance = gaussian_variance;
            likelyhood = gaussian_likelyhood;
            break;
        case Poisson:
            variance = poisson_variance;
            likelyhood = poisson_likelyhood;
            break;
        case Binomial:
            variance = binomial_variance;
            likelyhood = binomial_likelyhood;
            break;
        case Gamma:
            variance = gamma_variance;
            likelyhood = gamma_likelyhood;
            break;
        case InvGaussian:
            variance = inverse_gaussian_variance;
            likelyhood = inverse_gaussian_likelyhood;
            break;
    }

    switch (link_type) {
        case Identity:
            link_fun = identity_link;
            link_inv = identity_link_inv;
            derivative = identity_deriv;
            break;
        case Log:
            link_fun = log_link;
            link_inv = log_link_inv;
            derivative = log_deriv;
            break;
        case Logit:
            link_fun = logit_link;
            link_inv = logit_link_inv;
            derivative = logit_deriv;
            break;
        case Reciprocal:
            link_fun = reciprocal_link;
            link_inv = reciprocal_link_inv;
            derivative = reciprocal_deriv;
            break;
        case SquareReciprocal:
            link_fun = square_reciprocal_link;
            link_inv = square_reciprocal_link_inv;
            derivative = square_reciprocal_deriv;
            break;
    }
}