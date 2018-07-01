gee <- function(formula, data, id, family = gaussian(), weight = NULL,
                corstr = "independence", Mv = 0, cor.mat = matrix(),
                init.alpha = NULL, init.beta = NULL, init.phi = NULL, scale.fix = FALSE,
                penalty = FALSE, lambda = 10^-3, pindex = vector(), eps = 10^-6,
                maxit = 30, tol = 10^-6) {
  call <- match.call()

  # data
  # TODO data cleaning
  dat <- model.frame(formula, data)
  X <- model.matrix(formula, dat)
  Y <- model.response(dat)
  offset <- model.offset(dat)
  if (is.null(offset)) {
    offset <- rep(0, length(Y))
  }

  subj.col <- which(colnames(data) == call$id)
  id <- data[,subj.col]

  cluster.size <- as.numeric(summary(split(id, id, drop=T))[,1])
  max.cluster <- max(cluster.size)

  # structure
  if (sum(c("linkfun", "linkinv", "variance", "mu.eta") %in% names(family)) != 4) {
    stop("Problem with family parameter: should contains four functions")
  }

  cor.vec <- c("independence", "fixed", "ar1", "exchangeable", "m-dependent", "unstructured")
  cor.match <- charmatch(corstr, cor.vec)
  if (is.na(cor.match)) {
    stop("Unsupported correlation structure")
  } else if (cor.match == 5 & is.null(Mv)) {
    stop("Parameter Mv is not appropriate for m-dependent correlation")
  } else if (cor.match == 2 & ((is.null(cor.mat)) | 
            (dim(cor.mat)[1] != dim(cor.mat)[2]) |
            (dim(cor.mat)[1] > max.cluster))) {
    stop("Parameter cor.mat is not appropriate")
  }

  # initialize
  intercept.col <- apply(X==1, 2, all)
  if (is.null(init.beta)) {
    fit0 <- glm.fit(X, Y, offset=offset, family=family)
    init.beta <- fit0$coef
  } else if (length(init.beta) != dim(X)[2]) {
    stop("Length of init.beta is not correct.")
  }

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
    stop("Length of init.alpha is not correct.")
  }

  if (scale.fix & is.null(init.phi)) {
    stop("If scale.fix = TRUE, then init.phi must be supplied")
  }
  if (is.null(init.phi)) {
    init.phi = 1
  }

  if (is.null(weight)) {
    weight = rep(1, length(Y))
  }

  if (penalty) {
    result <- pgee_c(Y, X, offset, weight, cluster.size, family, corstr,
                    init.beta, init.alpha, init.phi, scale.fix, lambda, pindex, eps, maxit, tol, cor.mat, Mv)
  } else {
    result <- gee_c(Y, X, offset, weight, cluster.size, family, corstr,
                    init.beta, init.alpha, init.phi, scale.fix, maxit, tol, cor.mat, Mv)
  }

  result$call <- call
  result$corr.type <- corstr
  result$family.type <- family$family
  result$link.type <- family$link

  result$coefnames <- colnames(X)
  result$formula <- formula
  result$X <- X
  result$offset <- offset
  result$Y <- Y
  result$cluster.size <- cluster.size

  class(result) <- "geeq"
  return(result)
}