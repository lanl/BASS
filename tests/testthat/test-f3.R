test_that("f3 sobol & prediction matches previous validated values", {
  cat('Modified Friedman function (functional output) test (sobol & prediction) \n')

  eps<-1e-12

  mod<-readRDS('../f3_mod.rda') # previous model
  ss<-readRDS('../f3_sob.rda') # previous validated sobol
  sens<-sobol(mod,verbose = F) # new sobol
  diff<-max(abs(range(ss$var.tot-sens$var.tot))) # difference between total variance
  expect_that(diff,is_less_than(eps))

  diff<-max(abs(range(ss$S-sens$S))) # difference between sobol indices
  expect_that(diff,is_less_than(eps))

  ss.func<-readRDS('../f3_sobFunc.rda') # previous validated sobol - functional
  sens.func<-sobol(mod,func.var=1,verbose = F) # new sobol
  diff<-max(abs(range(ss.func$S.var-sens.func$S.var))) # difference between total variance
  expect_that(diff,is_less_than(eps))

  xtest<-readRDS('../f3_testX.rda') # x values
  pp<-predict(mod,xtest) # new predictions at x values
  pred<-readRDS('../f3_testPred.rda') # old predictions at x values
  diff<-max(abs(range(pred-pp))) # difference between new and old predictions
  expect_that(diff,is_less_than(eps))
})

