## R package for longitudinal data analysis

This repository is the result of the 2018 Google Summer of Code project "[GEE and QIF for clustered data regression](https://summerofcode.withgoogle.com/projects/#4738290838667264)".

### Introduction

`geeq` stands for  generalized estimating equations(GEE) and quadratic inference functions(QIF). Both methods focus on longitudinal or clustered data analysis, where within-cluster correlation has to be accounted for. Package `geeq` implements both of methods based on `RcppArmadillo` and also provided some new features.

### Why use `geeq`

**Simple:** `geeq` export an unified interface for GEE and QIF, users only need to provide `method` argument to specify which method to use. It's much convenience to compare two methods. 

**Extensible:** Nearly all computing process are written C++ based on package `RcppArmadillo`, which is a popular and easy to learn interface from R to and from Armadillo by utilising the Rcpp R/C++ interface library. Adding new features are feasible for ordinary users.

**Efficiency:** Armadillo is a high quality linear algebra library (matrix maths) for the C++ language, and `geeq` benefits a lot from it. For large data set, `geeq` is much more efficiency than `geeM` and `geepack`.

### How to use `geeq`

Currently, `geeq` provide only one interface for all things. Here is an example of how to use `geeq` to fit ohio data set from `geepack`:

```
> if (require(geepack)) {
  data('ohio', package='geepack')
  fit.geeq <- geeq(resp ~ age + smoke + age:smoke, id = id, data = ohio,
                method = 'gee', family = binomial, corstr = "exchangeable")
  fit.geepack <- geese(resp ~ age + smoke + age:smoke, id = id, data = ohio,
                family = binomial, corstr = "exchangeable")
}

> fit.geeq
            Estimates Model SE     wald       p
(Intercept)  -1.90000  0.11910 -15.9600 0.00000
age          -0.14120  0.05820  -2.4270 0.01524
smoke         0.31380  0.18780   1.6710 0.09478
age:smoke     0.07083  0.08828   0.8024 0.42230

 Correlation Model:  exchangeable
 Correlation Matrix:
       [,1]   [,2]   [,3]   [,4]
[1,] 1.0000 0.3544 0.3544 0.3544
[2,] 0.3544 1.0000 0.3544 0.3544
[3,] 0.3544 0.3544 1.0000 0.3544
[4,] 0.3544 0.3544 0.3544 1.0000

 Scale Parameter:  1.001

 Null deviance:  1829 on  2147  degrees of freedom
 Residual deviance:  1819 on  2144  degrees of freedom
 QIC:  1828.022

 Converged:  TRUE
 Number of iterations:  2
 Number of clusters:  537
 Maximum cluster size:  4

> fit.geepack
Call:
geese(formula = resp ~ age + smoke + age:smoke, id = id, data = ohio,
    family = binomial, corstr = "exchangeable")

Mean Model:
 Mean Link:                 logit
 Variance to Mean Relation: binomial

 Coefficients:
(Intercept)         age       smoke   age:smoke
-1.90049518 -0.14123591  0.31382579  0.07083184

Scale Model:
 Scale Link:                identity

 Estimated Scale Parameters:
(Intercept)
  0.9994078

Correlation Model:
 Correlation Structure:     exchangeable
 Correlation Link:          identity

 Estimated Correlation Parameters:
   alpha
0.354605

Returned Error Value:  0
Number of clusters:   537   Maximum cluster size: 4

```

There are some other examples in the `example` folder.

### Features in the future

* **Penalized Generalized Estimating Equations:** PGEE is important for variable selection. We've implemented solving penalized generalized estimating equations but not tuning paramter estimation. It's also essential to estimate the tuning paramter from users' perspective, so we decide not to export the interface now. The tuning parameter estimation part is in developing.

[Reference](https://github.com/zhouyuze/geeq/wiki/Reference)
