test_that("f1 sobol & prediction matches previous validated values", {
  cat('McKay function test (sobol & prediction) \n')

  eps<-1e-12

  mod<-readRDS('../f1_mod.rda') # previous model
  ss<-readRDS('../f1_sob.rda') # previous validated sobol
  sens<-sobol(mod,verbose = F) # new sobol
  diff<-max(abs(range(ss$var.tot-sens$var.tot))) # difference between total variance
  expect_that(diff,is_less_than(eps))

  diff<-max(abs(range(ss$S-sens$S))) # difference between sobol indices
  expect_that(diff,is_less_than(eps))

  xftest<-readRDS('../f1_testX.rda') # x values
  pp<-predict(mod,xftest) # new predictions at x values
  pred<-readRDS('../f1_testPred.rda') # old predictions at x values
  diff<-max(abs(range(pred-pp))) # difference between new and old predictions
  expect_that(diff,is_less_than(eps))
})

