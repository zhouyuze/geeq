## R package for longitudinal data analysis

### Introduction

`geeq` stands for  generalized estimating equations(GEE) and quadratic inference functions(QIF). Both methods focus on longitudinal or clustered data analysis, where within-cluster correlation has to be accounted for. Package `geeq` implements both of methods based on `RcppArmadillo` and also provided some new features.

### Why use `geeq`

**Simple:** `geeq` export an unified interface for GEE and QIF, users only need to provide `method` argument to specify which method to use. It's much convenience to compare two methods. 

**Extensible:** Nearly all computing process are written C++ based on package `RcppArmadillo`, which is a popular and easy to learn interface from R to and from Armadillo by utilising the Rcpp R/C++ interface library. Adding new features are feasible for ordinary users.

**Efficiency:** Armadillo is a high quality linear algebra library (matrix maths) for the C++ language, and `geeq` benefits a lot from it. For large data set, `geeq` is much more efficiency than `geeM` and `geepack`.

### How to use `geeq`

Currently, `geeq` provide only one interface for all things. Here is an example of how to use `geeq` to fit ohio data set from `geepack`:

```
if (require(geepack)) {
  data('ohio', package='geepack')
  fit <- geeq(resp ~ age + smoke + age:smoke, id=id, data=ohio,
              method='gee', family=binomial(), corstr="exchangeable")
}
```

There are some other examples in the `example` folder.

[Reference](https://github.com/zhouyuze/geeq/wiki/Reference)