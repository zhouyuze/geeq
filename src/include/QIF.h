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
    double phi;

    double Q = 0;
    vec Q_first_deriv;
    mat Q_second_deriv;

    void init_base_mat();
    void init_key_mat();
    double update_beta() override;
    void calculate_phi();
public:
    QIF(vec y, mat X, vec offset, uvec cluster_sizes,
        Family family, WorkCor cor_type, Control ctl, vec beta):
            Model(std::move(y), std::move(X), std::move(offset),
                  std::move(cluster_sizes), std::move(family),
                  cor_type, std::move(ctl), std::move(beta)),
            Q_first_deriv(p, fill::zeros), Q_second_deriv(p, p, fill::zeros) {
        init_base_mat();
        init_key_mat();
    }
    void iterator() override;
    Rcpp::List get_result() override;
};


#endif //SRC_QIF_H
