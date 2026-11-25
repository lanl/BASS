context("bass_sampler")

sampler_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(n * 3), n, 3)
  f <- sin(2 * pi * X[, 1]) + 2 * (X[, 2] - .5)^2
  list(X = X, f = f, y = f + rnorm(n, sd = .1))
}

test_that("bass_sampler fits, and as_bass() predicts its draws", {
  d <- sampler_data()
  s <- bass_sampler(d$X, d$y)
  s$step(2000)
  for (i in 1:100) { s$set_y(d$y); s$step(1); s$save() }
  expect_equal(s$n_saved(), 100L)
  fit <- s$as_bass()
  expect_s3_class(fit, "bass")
  P <- predict(fit, d$X)
  expect_equal(dim(P), c(100L, 200L))
  expect_equal(P[100, ], s$fitted(), tolerance = 1e-10)
  expect_lt(sqrt(mean((colMeans(P) - d$f)^2)), .1)
  expect_true(s$accept_rate() > 0 && s$accept_rate() < 1)
  expect_output(print(s), "BASS sampler")
})

test_that("set_y keeps the basis and recomputes the least-squares fit", {
  d <- sampler_data()
  s <- bass_sampler(d$X, d$y)
  s$step(200)
  nb <- s$nbasis
  s$set_y(-d$y)
  expect_identical(s$nbasis, nb)
  B <- s$basis
  bh <- c(solve(crossprod(B), crossprod(B, -d$y)))
  expect_equal(s$bhat, bh, tolerance = 1e-8)
  expect_error(s$set_y(d$y[-1]), "length")
})

test_that("without an intercept, centred draws average zero over the training x", {
  d <- sampler_data()
  s <- bass_sampler(d$X, d$y - mean(d$y), intercept = FALSE)
  expect_true(s$center.basis)
  s$save()                                  # an empty basis is a valid draw
  s$step(500)
  for (i in 1:19) { s$step(1); s$save() }
  P <- predict(s$as_bass(), d$X)
  expect_equal(unname(P[1, ]), rep(0, 200))
  expect_lt(max(abs(rowMeans(P))), 1e-10)
  expect_equal(P[20, ], s$fitted(), tolerance = 1e-10)
})

test_that("s2.fixed holds the residual variance", {
  d <- sampler_data()
  s <- bass_sampler(d$X, d$y, s2.fixed = .05)
  s$step(50); s$save()
  expect_identical(s$s2(), .05)
  expect_identical(s$as_bass()$s2, .05)
  expect_error(bass_sampler(d$X, d$y, s2.fixed = -1), "s2.fixed")
})
