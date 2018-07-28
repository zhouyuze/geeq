#' Solve Quadratic Inference Functions
#'
#' @param formula
#'
#' @param data
#'
#' @export
qif <- function(formula, id, data = parent.frame(), family = gaussian, weights = NULL,
                waves = NULL, corstr = "independence", Mv = 0, cor.mat = matrix(),
                init.beta = NULL, maxit = 30, tol = 10^-6) {
  call <- match.call()

  # data
  dat <- model.frame(formula, data, na.action=na.pass)

  if (class(data) == "data.frame") {
    subj.col <- which(colnames(data) == call$id)
    id <- data[, subj.col]

    if ("weights" %in% names(call)) {
      subj.col <- which(colnames(data) == call$weights)
      weights <- data[, subj.col]
    }

    if ("waves" %in% names(call)) {
      subj.col <- which(colnames(data) == call$waves)
      weights <-data[, subj.col]
    }
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

  result <- qif_c(Y, X, offset, weights, cluster.size, family, corstr, init.beta, maxit, tol)

  result$call <- call
  result$corr.type <- corstr
  names(result$beta) <- colnames(X)
  result$beta <- as.numeric(result$beta)
  result$family <- family


  result$formula <- formula
  result$weights <- weights
  result$X <- X
  result$offset <- offset
  result$Y <- Y
  result$mu <- result$X %*% result$beta
  result$cluster.size <- cluster.size

  class(result) <- "qif"
  return(result)
}
