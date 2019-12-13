test_that("f5 sobol (3 priors) & prediction matches previous validated values", {
  cat('Ishigami function test (sobol & prediction) multiple priors \n')

  eps<-1e-12

  mod<-readRDS('../f5_mod.rda') # previous model

  ss<-readRDS('../f5_sob1.rda') # previous validated sobol
  prior<-readRDS('../f5_sobPrior1.rda')
  sens<-sobol(mod,prior,verbose = F) # new sobol
  diff<-max(abs(range(ss$var.tot-sens$var.tot))) # difference between total variance
  expect_that(diff,is_less_than(eps))
  diff<-max(abs(range(ss$S-sens$S))) # difference between sobol indices
  expect_that(diff,is_less_than(eps))

  ss<-readRDS('../f5_sob2.rda') # previous validated sobol
  prior<-readRDS('../f5_sobPrior2.rda')
  sens<-sobol(mod,prior,verbose = F) # new sobol
  diff<-max(abs(range(ss$var.tot-sens$var.tot))) # difference between total variance
  expect_that(diff,is_less_than(eps))
  diff<-max(abs(range(ss$S-sens$S))) # difference between sobol indices
  expect_that(diff,is_less_than(eps))

  ss<-readRDS('../f5_sob3.rda') # previous validated sobol
  prior<-readRDS('../f5_sobPrior3.rda')
  sens<-sobol(mod,prior,verbose = F) # new sobol
  diff<-max(abs(range(ss$var.tot-sens$var.tot))) # difference between total variance
  expect_that(diff,is_less_than(eps))
  diff<-max(abs(range(ss$S-sens$S))) # difference between sobol indices
  expect_that(diff,is_less_than(eps))

  xtest<-readRDS('../f5_testX.rda') # x values
  pp<-predict(mod,xtest) # new predictions at x values
  pred<-readRDS('../f5_testPred.rda') # old predictions at x values
  diff<-max(abs(range(pred-pp))) # difference between new and old predictions
  expect_that(diff,is_less_than(eps))
})

