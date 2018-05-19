#ifndef SRC_GEE_H
#define SRC_GEE_H

#include "geeq.h"

class GEE_Control {
public:
    bool scale_fix;
    int maxit;
    double tol;
    bool trace;
    explicit GEE_Control(bool scale_fix = false, int maxit = 30, double tol = 10^-6, bool trace = false):
            scale_fix(scale_fix), maxit(maxit), tol(tol), trace(trace) {}
};

class Penalty_Options{
public:
    double lambda;
    uvec pindex;
    double eps;
    explicit Penalty_Options(double lambda, uvec pindex = uvec(), double eps = 10^-6):
            lambda(lambda), pindex(std::move(pindex)), eps(eps) {}
};
class GEE_Para {
private:
    // ****** passed variable ******
    // model data
    const vec y;
    const mat X;
    const vec offset;
    const uvec cluster_bound;

    // model structure
    Family funcs;
    const WorkCor cor_type;
    const GEE_Control ctl;

    // model parameter
    vec beta;
    vec alpha;
    double phi;

    // ****** calculated variable ******
    bool solved;
    bool converged;

    const int N; // number of all observations
    const int n; // number of clusters
    const int p; // number of beta

    // intermediate result; initialize in Constructor and update when beta is updated
    vec eta;
    vec mu;
    vec var;
    vec deriv;
    // correlation in each cluster
    vector<mat> cluster_cor;

    void update_phi();
    void update_alpha();
    double update_beta();
    double update_beta_penalty(Penalty_Options op);
    void update_intermediate_result();

    vec q_scad(double lambda, double a = 3.7);

public:
    GEE_Para(vec y, mat X, vec offset, uvec cluster_sizes,
             Family family, WorkCor cor_type, GEE_Control ctl,
             vec beta, vec alpha, double phi);
    int iterator();
    int iterator_penalty(Penalty_Options op);
    vec get_beta();
    vec get_alpha();
    double get_phi();
    mat get_sandwich();

    double gaussian_pseudolikelihood();
    vec geodesic_distance();
};



#endif //SRC_GEE_H
