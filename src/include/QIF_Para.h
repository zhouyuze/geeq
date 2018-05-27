#ifndef SRC_QIF_H
#define SRC_QIF_H

#include "geeq.h"

class QIF_Para {
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

    // model parameter
    vec beta;

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

    vector<mat> base_mat;
    vec g;
    mat C;
    mat dev_g;

    double update_beta();
    void update_intermediate_result();

public:
    QIF_Para(vec y, mat X, vec offset, uvec cluster_sizes,
             Family family, WorkCor cor_type, vec beta);
    int iterator();
    vec get_beta();
};


#endif //SRC_QIF_H
