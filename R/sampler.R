########################################################################
## Persistent, updatable BASS sampler (one block of a larger Gibbs sampler)
########################################################################

#' @title Updatable BASS Sampler
#'
#' @description Sets up the state of a BASS Markov chain (continuous inputs,
#'   scalar response, no tempering) so that it can be advanced a few iterations
#'   at a time, with a response that may change between calls. Use it when BASS
#'   is one block of a larger Gibbs sampler whose response is itself a latent
#'   variable that is redrawn every outer iteration.
#'
#' @param xx a numeric matrix (or data frame of numeric columns) of predictors.
#' @param y the starting response vector.
#' @param maxInt,maxBasis,npart,g1,g2,h1,h2,a.tau,b.tau,w1,w2 as in \link{bass}
#'   (\code{degree} is 1).
#' @param intercept logical; if \code{TRUE} (default) the model has a constant
#'   basis function.
#' @param center.basis logical; if \code{TRUE}, every basis function is centred
#'   at its mean over the training inputs. The default is \code{!intercept}:
#'   without an intercept, every draw of \eqn{f} then averages exactly zero over
#'   the training inputs.
#' @param s2.fixed if not \code{NULL}, the residual variance is fixed at this
#'   value and not sampled.
#'
#' @details The chain state (the basis functions, their cross-products and the
#'   current coefficients and variances) is kept in an environment. Only
#'   \eqn{B'y}, \eqn{y'y} and the quantities derived from them depend on the
#'   response, so \code{$set_y()} does not rebuild the basis.
#'
#'   Each iteration is the BASS update: a birth, death or change of a basis
#'   function given \eqn{(\sigma^2, \tau, \lambda)}, with the coefficients
#'   integrated out, then Gibbs draws of \eqn{\sigma^2} (coefficients
#'   integrated out), the coefficients, \eqn{\lambda} and \eqn{\tau}. That
#'   kernel leaves the BASS posterior given \eqn{y} invariant, so interleaving
#'   it with updates of \eqn{y} is a valid Gibbs sampler for the joint model.
#'   Unlike \link{bass}, the Gibbs draws run even when the proposal is invalid,
#'   so the coefficients always reflect the current response.
#'
#'   The returned environment has these methods:
#'   \describe{
#'     \item{\code{set_y(y)}}{replace the response.}
#'     \item{\code{iterate()}}{one RJMCMC move and the Gibbs draws.}
#'     \item{\code{step(k = 1)}}{\code{k} calls of \code{iterate()}.}
#'     \item{\code{fitted()}}{the current mean function at the training inputs.}
#'     \item{\code{s2()}}{the current residual variance.}
#'     \item{\code{save()}}{append the current state to the stored draws.}
#'     \item{\code{n_saved()}}{number of stored draws.}
#'     \item{\code{accept_rate()}}{fraction of accepted RJMCMC proposals.}
#'     \item{\code{as_bass()}}{the stored draws as a \code{"bass"} object,
#'       usable with \link{predict.bass} and \link{sobol}.}
#'   }
#'   A centred basis function \eqn{B_m(x) - \bar B_m} is stored as \eqn{B_m}
#'   with its constant folded into the intercept, so \code{as_bass()} predicts
#'   with the usual BASS basis.
#' @return An environment of class \code{"bass_sampler"}.
#' @seealso \link{bass}.
#' @export
#' @examples
#' set.seed(1)
#' n <- 100
#' X <- matrix(runif(n * 2), n, 2)
#' s <- bass_sampler(X, X[, 1]^2 + rnorm(n, sd = 0.1))
#' s$step(200)
#' for(it in 1:50){
#'   s$set_y(X[, 1]^2 + rnorm(n, sd = 0.1)) # a new response every iteration
#'   s$step(1)
#'   s$save()
#' }
#' fit <- s$as_bass()
#' dim(predict(fit, X))
#'
bass_sampler <- function(xx, y, maxInt = 3, maxBasis = 1000, npart = NULL, g1 = 0, g2 = 0, h1 = 10,
                         h2 = 10, a.tau = .5, b.tau = NULL, w1 = 5, w2 = 5, intercept = TRUE,
                         center.basis = !intercept, s2.fixed = NULL){
  self <- environment()

  xx <- as.matrix(xx)
  if(!is.numeric(xx))
    stop('xx must be numeric (bass_sampler has no categorical inputs)')
  y <- as.numeric(y)
  n <- length(y)
  if(nrow(xx) != n)
    stop('nrow(xx) must equal length(y)')
  if(anyNA(xx) || anyNA(y))
    stop('xx and y must not contain missing values')
  if(!posInt(maxInt) || !posInt(maxBasis))
    stop('invalid maxInt or maxBasis')
  if(any(c(g1, g2) < 0))
    stop('g1 and g2 must be greater than or equal to 0')
  if(!is.null(b.tau) && b.tau <= 0)
    stop('b.tau must be greater than 0')
  if(any(c(h1, h2, a.tau, w1, w2) <= 0))
    stop('h1,h2,a.tau,w1,w2 must be greater than 0')
  if(!is.logical(intercept) || length(intercept) != 1 || is.na(intercept))
    stop('intercept must be TRUE or FALSE')
  if(!is.null(s2.fixed) && !(is.numeric(s2.fixed) && length(s2.fixed) == 1 && s2.fixed > 0))
    stop('s2.fixed must be NULL or a positive number')

  p <- ncol(xx)
  range.des <- apply(xx, 2, range)
  xx.des <- apply(xx, 2, scale_range)
  if(is.null(dim(xx.des)))
    xx.des <- matrix(xx.des, n, p)
  xxt <- t(xx.des)
  unique.ind <- lapply(1:p, function(j) which(!duplicated(xxt[j, ])))
  vars.len <- sapply(unique.ind, length)

  if(is.null(npart))
    npart <- min(20, .1 * n)
  maxInt <- min(maxInt, p)
  maxBasis <- min(maxBasis, n)
  if(is.null(b.tau))
    b.tau <- 2 / n
  prior <- list(miC = 0) # minInt = 1, as for bass() with continuous inputs only
  i0 <- as.integer(intercept)

  ## chain state
  ssy <- sum(y^2)
  nbasis <- 0L
  nc <- i0
  basis <- matrix(1, n, i0)          # columns: [intercept,] basis functions (centred if center.basis)
  bmean <- numeric(0)                # training mean of each (uncentred) basis function
  vars <- signs <- knotInd <- matrix(NA_integer_, 0, maxInt)
  knots <- matrix(NA_real_, 0, maxInt)
  n.int <- integer(0)
  XtX <- crossprod(basis)
  Xty <- c(crossprod(basis, y))
  I.star <- rep(w1, maxInt)
  I.vec <- I.star / sum(I.star)
  z.star <- rep(w2, p)
  z.vec <- z.star / sum(z.star)
  s2_cur <- if(is.null(s2.fixed)) 1 else s2.fixed
  lam <- 1
  tau <- 1
  beta <- numeric(nc)
  n_prop <- 0
  n_acc <- 0
  cmod <- TRUE                       # basis changed since the last save()

  qf_info <- function(XtX, Xty){
    if(length(Xty) == 0)
      return(list(R = matrix(0, 0, 0), bhat = numeric(0), qf = 0))
    out <- getQf(XtX, Xty)
    if(is.null(out))
      return(NULL)
    out$qf <- c(out$qf)
    out$bhat <- c(out$bhat)
    out
  }
  cur <- qf_info(XtX, Xty)
  R <- cur$R
  bhat <- cur$bhat
  qf <- cur$qf

  accept <- function(q){
    R <<- q$R
    bhat <<- q$bhat
    qf <<- q$qf
    n_acc <<- n_acc + 1
    cmod <<- TRUE
  }

  prep <- function(b){
    m <- mean(b)
    if(center.basis) b <- b - m
    list(b = b, m = m)
  }

  birth <- function(){
    cand <- genCandBasis(minInt = 1, maxInt = maxInt, I.vec = I.vec, z.vec = z.vec, p = p, xxt = xxt,
                         q = 1, xx.unique.ind = unique.ind, vars.len = vars.len, prior = prior)
    if(sum(cand$basis != 0) < npart)
      return(invisible(NULL))
    pb <- prep(cand$basis)
    Xta <- c(crossprod(basis, pb$b))
    XtX.cand <- rbind(cbind(XtX, Xta), c(Xta, sum(pb$b^2)))
    Xty.cand <- c(Xty, sum(pb$b * y))
    q <- qf_info(XtX.cand, Xty.cand)
    if(is.null(q))
      return(invisible(NULL))
    alpha <- .5 / s2_cur * (q$qf - qf) / (1 + tau) + log(lam) - log(nbasis + 1) - cand$lbmcmp +
      .5 * log(tau) - .5 * log(1 + tau)
    if(log(runif(1)) < alpha){
      accept(q)
      XtX <<- unname(XtX.cand)
      Xty <<- Xty.cand
      basis <<- cbind(basis, pb$b, deparse.level = 0)
      bmean <<- c(bmean, pb$m)
      fill <- rep(NA, maxInt - cand$n.int)
      vars <<- rbind(vars, c(cand$vars, fill))
      signs <<- rbind(signs, c(cand$signs, fill))
      knotInd <<- rbind(knotInd, c(cand$knotInd, fill))
      knots <<- rbind(knots, c(cand$knots, fill))
      n.int <<- c(n.int, cand$n.int)
      nbasis <<- nbasis + 1L
      nc <<- nbasis + i0
      I.star[cand$n.int] <<- I.star[cand$n.int] + 1
      I.vec <<- I.star / sum(I.star)
      z.star[cand$vars] <<- z.star[cand$vars] + 1
      z.vec <<- z.star / sum(z.star)
    }
  }

  death <- function(){
    k <- sample.int(nbasis, 1)
    ind <- (1:nc)[-(k + i0)]
    q <- qf_info(XtX[ind, ind, drop = FALSE], Xty[ind])
    if(is.null(q))
      return(invisible(NULL))
    vk <- vars[k, 1:n.int[k]]
    I.star1 <- I.star
    I.star1[n.int[k]] <- I.star1[n.int[k]] - 1
    I.vec1 <- I.star1 / sum(I.star1)
    z.star1 <- z.star
    z.star1[vk] <- z.star1[vk] - 1
    z.vec1 <- z.star1 / sum(z.star1)
    lpbmcmp <- logProbChangeMod(n.int[k], vk, I.vec1, z.vec1, p, vars.len, maxInt, prior$miC)
    alpha <- .5 / s2_cur * (q$qf - qf) / (1 + tau) - log(lam) + log(nbasis) + lpbmcmp -
      .5 * log(tau) + .5 * log(1 + tau)
    if(log(runif(1)) < alpha){
      accept(q)
      XtX <<- XtX[ind, ind, drop = FALSE]
      Xty <<- Xty[ind]
      basis <<- basis[, ind, drop = FALSE]
      bmean <<- bmean[-k]
      vars <<- vars[-k, , drop = FALSE]
      signs <<- signs[-k, , drop = FALSE]
      knotInd <<- knotInd[-k, , drop = FALSE]
      knots <<- knots[-k, , drop = FALSE]
      n.int <<- n.int[-k]
      nbasis <<- nbasis - 1L
      nc <<- nbasis + i0
      I.star <<- I.star1
      I.vec <<- I.vec1
      z.star <<- z.star1
      z.vec <<- z.vec1
    }
  }

  change <- function(){
    k <- sample.int(nbasis, 1)
    use <- 1:n.int[k]
    int.change <- use[sample.int(length(use), 1)]
    cand <- genBasisChange(NULL, k, int.change, xxt, 1, knots = knots[k, use], knotInd = knotInd[k, use],
                           signs = signs[k, use], vars = vars[k, use], xx.unique.ind = unique.ind)
    if(sum(cand$basis != 0) < npart)
      return(invisible(NULL))
    pb <- prep(cand$basis)
    j <- k + i0
    XtX.cand <- XtX
    XtX.cand[j, ] <- XtX.cand[, j] <- c(crossprod(basis, pb$b))
    XtX.cand[j, j] <- sum(pb$b^2)
    Xty.cand <- Xty
    Xty.cand[j] <- sum(pb$b * y)
    q <- qf_info(XtX.cand, Xty.cand)
    if(is.null(q))
      return(invisible(NULL))
    alpha <- .5 / s2_cur * (q$qf - qf) / (1 + tau)
    if(log(runif(1)) < alpha){
      accept(q)
      XtX <<- XtX.cand
      Xty <<- Xty.cand
      basis[, j] <<- pb$b
      bmean[k] <<- pb$m
      signs[k, use] <<- cand$signs
      knotInd[k, use] <<- cand$knotInd
      knots[k, use] <<- cand$knots
    }
  }

  gibbs <- function(){
    if(is.null(s2.fixed)){
      rate <- g2 + .5 * (ssy - sum(bhat * Xty) / (1 + tau))
      rate <- max(rate, .Machine$double.eps)
      s2_cur <<- 1 / rgamma(1, g1 + n / 2, rate)
    }
    if(nc > 0){
      beta <<- bhat / (1 + tau) + c(backsolve(R, rnorm(nc))) * sqrt(s2_cur / (1 + tau))
    } else {
      beta <<- numeric(0)
    }
    lam <<- rgamma(1, h1 + nbasis, h2 + 1)
    qf2 <- if(nc > 0) sum((R %*% beta)^2) else 0
    tau <<- rgamma(1, a.tau + nc / 2, b.tau + .5 * qf2 / s2_cur)
    invisible(self)
  }

  set_y <- function(y_new){
    y_new <- as.numeric(y_new)
    if(length(y_new) != n)
      stop(sprintf('y must have length %d', n))
    if(anyNA(y_new))
      stop('y must not contain missing values')
    y <<- y_new
    ssy <<- sum(y^2)
    if(nc > 0){
      Xty <<- c(crossprod(basis, y))
      bhat <<- c(backsolve(R, forwardsolve(t(R), Xty)))
      qf <<- sum(bhat * Xty)
    }
    invisible(self)
  }

  iterate <- function(){
    n_prop <<- n_prop + 1
    if(nbasis == 0){
      birth()
    } else {
      u <- sample.int(3, 1)
      if(nbasis >= maxBasis) u <- sample(2:3, 1)
      if(u == 1) birth() else if(u == 2) death() else change()
    }
    gibbs()
  }

  step <- function(k = 1L){
    for(i in seq_len(k)) iterate()
    invisible(self)
  }

  fitted <- function() if(nc > 0) c(basis %*% beta) else numeric(n)
  s2 <- function() s2_cur
  accept_rate <- function() n_acc / max(n_prop, 1)

  ## stored draws
  keep <- list(beta = list(), s2 = numeric(0), lam = numeric(0), tau = numeric(0),
               nbasis = integer(0), model.lookup = integer(0), models = list())

  save <- function(){
    # draw in the usual BASS form: an intercept, then the uncentred basis functions
    b0 <- if(intercept) beta[1] else 0
    bb <- beta[i0 + seq_len(nbasis)]
    if(center.basis) b0 <- b0 - sum(bb * bmean)
    if(cmod){
      keep$models[[length(keep$models) + 1L]] <<- list(vars = vars, signs = signs, knotInd = knotInd,
                                                       n.int = n.int)
      cmod <<- FALSE
    }
    k <- length(keep$s2) + 1L
    keep$beta[[k]] <<- c(b0, bb)
    keep$s2[k] <<- s2_cur
    keep$lam[k] <<- lam
    keep$tau[k] <<- tau
    keep$nbasis[k] <<- nbasis
    keep$model.lookup[k] <<- length(keep$models)
    invisible(self)
  }
  n_saved <- function() length(keep$s2)

  as_bass <- function(){
    ns <- n_saved()
    if(ns == 0)
      stop('no draws have been saved')
    nm <- length(keep$models)
    mb <- max(1L, keep$nbasis)
    beta.mat <- matrix(NA_real_, ns, mb + 1)
    for(k in 1:ns) beta.mat[k, seq_along(keep$beta[[k]])] <- keep$beta[[k]]
    arr <- function(f) {
      a <- array(NA_integer_, c(nm, mb, maxInt))
      for(i in 1:nm){
        m <- keep$models[[i]][[f]]
        if(nrow(m) > 0) a[i, 1:nrow(m), ] <- m
      }
      a
    }
    n.int.mat <- matrix(NA_integer_, nm, mb)
    for(i in 1:nm){
      ni <- keep$models[[i]]$n.int
      if(length(ni) > 0) n.int.mat[i, seq_along(ni)] <- ni
    }
    out <- list(beta = beta.mat, s2 = keep$s2, lam = keep$lam, nbasis = keep$nbasis, degree = 1,
                nmcmc = ns, nburn = 0, thin = 1, p = p, beta.prec = keep$tau, y = y,
                n.models = nm, model.lookup = keep$model.lookup, des = TRUE, func = FALSE,
                cat = FALSE, type = '_des', cx = rep('numeric', p),
                knotInd.des = arr('knotInd'), signs.des = arr('signs'), vars.des = arr('vars'),
                n.int.des = n.int.mat, maxInt.des = maxInt, pdes = p, xx.des = xx.des,
                range.des = range.des, unique.ind.des = unique.ind)
    class(out) <- 'bass'
    out
  }

  gibbs() # so the coefficients exist before the first step
  class(self) <- 'bass_sampler'
  self
}

#' @export
print.bass_sampler <- function(x, ...){
  cat(sprintf('BASS sampler: n = %d, p = %d, %d basis functions, %d saved draws, acceptance %.3f\n',
              x$n, x$p, as.integer(x$nbasis), x$n_saved(), x$accept_rate()))
  invisible(x)
}
