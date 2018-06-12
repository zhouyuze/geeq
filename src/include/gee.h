#ifndef SRC_GEE_H
#define SRC_GEE_H

#include "model.h"

class Penalty_Options{
public:
    double lambda;
    uvec pindex;
    double eps;
    explicit Penalty_Options(double lambda, uvec pindex = uvec(), double eps = 10^-6):
            lambda(lambda), pindex(std::move(pindex)), eps(eps) {}
};

class GEE : public Model {
protected:
    vec alpha;
    double phi;
    bool scale_fix;
    int Mv;

    // correlation in each cluster
    vector<mat> cluster_cor;

    mat H1;
    mat H2;
    vec score;

    void update_phi();
    void update_alpha();

    double update_beta() override;
    double update_beta_penalty(Penalty_Options op);

    vec q_scad(double lambda, double a = 3.7);

    void calculate_H2();
    mat get_sandwich();
    double gaussian_pseudolikelihood();
    vec geodesic_distance();

public:
    GEE(vec y, mat X, vec offset, uvec cluster_sizes,
        Family family, WorkCor cor_type, Control ctl,
        vec beta, vec alpha, double phi, bool fix, mat cor_mat = zeros(1), int Mv = 0);

    int iterator() override;
    int iterator_penalty(Penalty_Options op);
    Rcpp::List get_result() override;
};



#endif //SRC_GEE_H
