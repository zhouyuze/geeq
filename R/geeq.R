#' Solve Generalized Estimating Equations and Quadratic Inference Function
#'
#' @description
#' Calculate coefficients of generalized estimating equations. Both generalised estimating
#' equations and quadratic inference function method are provided. Produces an object of
#' class 'geeq' which is a fit of the model.
#'
#' @param formula a formula expression similar to that for \code{\link[stats]{glm}}, of the
#' form \code{response ~ predictors}. A term of the form \code{offset(expression)} is allowed.
#'
#' @param id a vector which identifies the clusters. The length of \code{id} should be the same
#' as the number of observations. Observations with same \code{id} belong to the same cluster.
#'
#' @param data an optional data frame containing the variables and \code{id}. \code{weights} and
#' \code{waves} should also be included if necessary.
#'
#' @param family \code{family} argument determine the details of the model. The argument can be
#' one of three options: a \code{\link[stats]{family}} object, a character string corresponding
#' to one of the family objects, or a list of user specified functions.
#'
#' User specified functions must contain following components: \code{linkfun}, \code{linkinv},
#' \code{variance}, \code{mu.eta} and \code{dev.resids}. See \code{\link[stats]{family}} for more
#' information.
#'
#' @param method a character string specifying the method to estimate. Options could be 'gee' of 'qif'.
#'
#' @param weights an optional vector for each observation. Observations with weight 0 are excluded
#' from the calculation process. All observations will be assigned 1 by default. If an observation
#' misses any variable, it will be assigned 0.
#'
#' @param waves an optional vector identifying time ordering in a cluster.
#'
#' @param maxit integer giving the maximal number of iteration.
#'
#' @param tol positive convergence tolerance. The iterations converge when the absolute value
#' of the difference in parameter estimate is below \code{tol}
#'
#' @param corstr a character string specifying the correlation structure. Options could be:
#' 'independence', 'fixed', 'ar1', 'exchangeable', 'm-dependent' and 'unstructured'.
#'
#' @param Mv for 'm-dependent', the value for \code{m}.
#'
#' @param corr.mat the correlations matrix for 'fixed'. Matrix should be symmetric with dimensions
#' >= the maximum cluster size.
#'
#' @param init.alpha an optional vector with initial values of alpha. The length is different for
#' different correlation structure.
#'
#' @param init.beta an optional vector with initial values of beta. If not specificed, then it will
#' be set to the result of \code{\link[stats]{glm}}
#'
#' @param init.phi an optional initial overdispersion parameter.
#'
#' @param scale.fix if set to \code{TRUE}, then the scale parameter is fixed at the value of \code{
#' init.phi}.
#'
#' @return
#' An object of class 'geeq' representing the fit. Some of the attributes are defined as follows:
#' \item{\code{beta}}{a named vector of estimated coefficients.}
#' \item{\code{variance}}{estimated variance of \code{beta}.}
#' \item{\code{phi}}{estimated scale parameter.}
#' \item{\code{fitted.values}}{the fitted mean values, obtained by transforming the linear predictors
#' by the inverse of the link function.}
#' \item{\code{call}}{the matched call.}
#' \item{\code{model}}{the model frame.}
#' \item{\code{offset}}{the priori known component in the linear predictor.}
#' \item{\code{df.null}}{the residual degrees of freedom for the null model.}
#' \item{\code{df.residual}}{the residual degrees of freedom.}
#' \item{\code{null.deviance}}{the deviance for the null model, comparable with deviance. The null model
#' will include the offset, and an intercept if there is one in the model.}
#' \item{\code{deviance}}{up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant
#' is chosen so that a saturated model has deviance zero.}
#' \item{\code{niter}}{the number of iterations to solve equations.}
#' \item{\code{converged}}{logical. Whether the algorithm have converged.}
#' \item{\code{cluster.size}}{a vector of size for each cluster}
#' \item{\code{ncluster}}{number of clusters}
#' \item{\code{nobs}}{number of observations from all clusters}
#'
#' For the result of 'gee' method, the following attributes are provided:
#' \item{\code{alpha}}{a vector of estimated correlation parameters.}
#' \item{\code{correlation}}{a matrix of estimated correlation.}
#' \item{\code{QIC}}{a criteria to select correlation structure.}
#'
#' As for 'qif' method, the result has following attributes:
#' \item{\code{Q}}{value of quadratic inference function}
#' \item{\code{statistics}}{some statistics to evaluate the result}
#'
#' @author
#' Yuze Zhou
#'
#' @seealso
#' \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#'
#' @examples
#' if (require(geepack)) {
#'   data('ohio', package='geepack')
#'   fit <- geeq(resp ~ age + smoke + age:smoke, id=id, data=ohio,
#'               method='gee', family=binomial(), corstr="exchangeable")
#' }
#'
#' @export
geeq <- function(formula, id = NULL, data = parent.frame(), family = gaussian, method = 'gee',
                weights = NULL, waves = NULL, maxit = 30, tol = 10^-6,
                corstr = 'independence', Mv = 1, cor.mat = matrix(),
                init.alpha = NULL, init.beta = NULL, init.phi = 1, scale.fix = FALSE) {
  call <- match.call()

  # data
  dat <- model.frame(formula, data, na.action=na.pass)

  if (class(data) == 'data.frame') {
    if ('id' %in% names(call)) {
        subj.col <- which(colnames(data) == call$id)
        id <- data[, subj.col]
    }

    if ('weights' %in% names(call)) {
      subj.col <- which(colnames(data) == call$weights)
      weights <- data[, subj.col]
    }

    if ('waves' %in% names(call)) {
      subj.col <- which(colnames(data) == call$waves)
      weights <-data[, subj.col]
    }
  }

  if (is.null(id)) {
    id = seq(1, dim(dat)[1])
  }

  if (is.null(weights)) {
    weights = rep(1, dim(dat)[1])
  }

  na.inds <- NULL

  if (any(is.na(dat))) {
    na.inds <- which(is.na(dat), arr.ind = T)
  }

  if (!is.null(waves)) {
    dat <- dat[order(id, waves),]
  } else {
    dat <- dat[order(id),]
  }

  if (!is.null(na.inds)) {
    for (i in 1:dim(na.inds)[1]) {
      if (is.factor(dat[na.inds[i, 1], na.inds[i, 2]])) {
        dat[na.inds[i, 1], na.inds[i, 2]] <- levels(dat[, na.inds[i, 2]])[1]
      } else {
        dat[na.inds[i, 1], na.inds[i, 2]] <- median(dat[, na.inds[i, 2]], na.rm=T)
      }
    }
    weights[unique(na.inds[,1])] <- 0
  }

  X <- model.matrix(formula, dat)
  Y <- model.response(dat)
  offset <- model.offset(dat)
  if (is.null(offset)) {
    offset <- rep(0, length(Y))
  }
  intercept <- attr(attr(dat, 'terms'), 'intercept') > 0L

  cluster.size <- as.numeric(summary(split(id, id, drop=T))[,1])
  max.cluster <- max(cluster.size)

  # structure
  if (is.character(family)) {
    family <- get(family, mode = 'function', envir = parent.frame(2))
  } else if (is.function(family)) {
    family <- family()
  }
  if (sum(c('linkfun', 'linkinv', 'variance', 'mu.eta', 'dev.resids') %in% names(family)) != 5) {
    stop('Problem with family parameter: should contains four functions')
  }

  cor.vec <- c('independence', 'fixed', 'ar1', 'exchangeable', 'm-dependent', 'unstructured')
  cor.match <- charmatch(corstr, cor.vec)
  if (is.na(cor.match)) {
    stop('Unsupported correlation structure')
  } else if (cor.match == 5 & is.null(Mv)) {
    stop('Parameter Mv is not appropriate for m-dependent correlation')
  } else if (cor.match == 2 & ((is.null(cor.mat)) |
            (dim(cor.mat)[1] != dim(cor.mat)[2]) |
            (dim(cor.mat)[1] > max.cluster))) {
    stop('Parameter cor.mat is not appropriate')
  }

  if (method == 'qif' & (cor.match == 2 | cor.match == 5)) {
    stop('Unsupported correlation structure for qif')
  }

  # initialize
  if (is.null(init.beta)) {
    fit0 <- glm.fit(X, Y, offset=offset, family=family)
    init.beta <- fit0$coef
  } else if (length(init.beta) != dim(X)[2]) {
    stop('Length of init.beta is not correct.')
  }

  if (method == 'gee') {
    if (cor.match < 3) {
      alpha.length <- 0
    } else if (cor.match < 5) {
      alpha.length <- 1
    } else if (cor.match == 5) {
      alpha.length <- Mv
    } else {
      alpha.length <- sum(1:(max.cluster - 1))
    }

    if (is.null(init.alpha)) {
      init.alpha <- rep(0.2, alpha.length)
    } else if (length(init.alpha) != alpha.length) {
      stop('Length of init.alpha is not correct.')
    }
  }

  if (method == 'gee') {
    result <- gee_c(Y, X, offset, weights, cluster.size, family, corstr,
                    init.beta, init.alpha, init.phi, scale.fix, maxit, tol, cor.mat, Mv)
    result$alpha <- as.numeric(result$alpha)
  } else if (method == 'qif') {
    result <- qif_c(Y, X, offset, weights, cluster.size, family, corstr, init.beta, maxit, tol)
  }

  result$call <- call
  result$method <- method
  result$formula <- formula
  result$corr.type <- corstr
  result$beta <- as.numeric(result$beta)
  result$family <- family
  names(result$beta) <- colnames(X)
  colnames(result$variance) <- colnames(X)
  rownames(result$variance) <- colnames(X)

  result$model <- dat
  result$offset <- offset
  result$id <- id
  result$weights <- weights
  result$fitted.values <- as.vector(family$linkinv(X %*% result$beta + offset))
  result$cluster.size <- cluster.size
  result$ncluster <- length(cluster.size)
  result$nobs <- dim(dat)[1]

  if (method == 'qif') {
      np = length(result$beta)
      pvalue = 1 - pchisq(result$Q, np)

      if (corstr == 'independence') {
        AIC <- result$Q
        BIC <- result$Q
      } else {
        AIC <- result$Q + 2 * np
        BIC <- result$Q + np * log(result$nobs)
      }
      result$statistics <- c(result$Q, np, pvalue, AIC, BIC)
      names(result$statistics) <- c('Q', 'D.F.', 'pvalue', 'AIC', 'BIC')
  }

  n.ok <- result$nobs - sum(weights==0)
  result$df.null <- n.ok - as.integer(intercept)
  result$df.residual <- n.ok - dim(X)[2]

  wtdmu <- if (intercept) sum(weights * Y) / sum(weights) else family$linkinv(offset)
  result$null.deviance <- sum(family$dev.resids(Y, wtdmu, weights))
  result$deviance <- sum(family$dev.resids(Y, result$fitted.values, weights))

  class(result) <- 'geeq'
  return(result)
}
