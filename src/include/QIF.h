#ifndef SRC_QIF_H
#define SRC_QIF_H

#include "model.h"

class QIF : public Model {
protected:
    vector<mat> base_mat;
    int m = 0; // number of base matrix
    vec g;
    mat C;
    mat dev_g;

    double Q = 0;
    vec Q_first_deriv;
    mat Q_second_deriv;

    void init_base_mat();
    void init_key_mat();
    double update_beta() override;
public:
    QIF(vec y, mat X, vec offset, uvec cluster_sizes,
        Family family, WorkCor cor_type, Control ctl, vec beta);
    void iterator() override;
    Rcpp::List get_result() override;
};


#endif //SRC_QIF_H
