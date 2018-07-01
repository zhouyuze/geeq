#ifndef SRC_MODEL_H
#define SRC_MODEL_H

#include "Family.h"

class Model {
protected:
    // ****** passed variable ******
    // model data
    const vec y;
    const mat X;
    const vec offset;
    const vec weight;
    const uvec cluster_bound;

    // model structure
    Family funcs;
    const WorkCor cor_type;
    const Control ctl;

    // model parameter
    vec beta;

    // ****** calculated variable ******
    const int N; // number of all observations
    const int n; // number of clusters
    const int p; // number of beta
    const int max_cluster; // maximum cluster size

    // intermediate result; initialize in Constructor and update when beta is updated
    vec eta;
    vec mu;
    vec var;
    vec deriv;

    bool converged = false;
    int niter = 0;
    bool balanced; // whether cluster sizes are equal

    void update_intermediate_variable();
    virtual double update_beta() = 0;
public:
    Model() = delete;
    Model(vec y, mat X, vec offset, vec weight, uvec cluster_sizes,
          Family family, WorkCor cor_type, Control ctl, vec init_beta);
    ~Model() = default;

    virtual void iterator() = 0;
    virtual Rcpp::List get_result() = 0;

};
#endif //SRC_MODEL_H
