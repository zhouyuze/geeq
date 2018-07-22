#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

class Control {
public:
    int maxit;
    double tol;
    bool trace;
    explicit Control(int maxit = 30, double tol = 10^-6, bool trace = false):
            maxit(maxit), tol(tol), trace(trace) {}
};

enum WorkCor {
    Independence,
    AR1,
    Exchangable,
    M_dependent,
    Fixed,
    Unstructured
};

const static map<string, WorkCor> workCorMap = {
        {"independence", Independence},
        {"ar1", AR1},
        {"exchangeable", Exchangable},
        {"m-dependent", M_dependent},
        {"fixed", Fixed},
        {"unstructured", Unstructured}
};

#endif //SRC_UTILS_H
