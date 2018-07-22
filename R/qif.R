qif <- function(formula, data, id, family = gaussian(), weight = NULL,
                corstr = "independence", Mv = 0, cor.mat = matrix(),
                init.beta = NULL, maxit = 30, tol = 10^-6) {
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
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame(2))
  } else if (is.function(family)) {
    family <- family()
  }
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
  if (is.null(init.beta)) {
    fit0 <- glm.fit(X, Y, offset=offset, family=family)
    init.beta <- fit0$coef
  } else if (length(init.beta) != dim(X)[2]) {
    stop("Length of init.beta is not correct.")
  }

  if (is.null(weight)) {
    weight = rep(1, length(Y))
  }

  result <- qif_c(Y, X, offset, weight, cluster.size, family, corstr, init.beta, maxit, tol)

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

  class(result) <- "qif"
  return(result)
}