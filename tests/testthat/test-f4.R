test_that("f4 sobol & prediction matches previous validated values", {
  cat('Modified Friedman function (functional output, PCA approach) test (sobol & prediction) \n')

  eps<-1e-12

  mod<-readRDS('../f4_mod.rda') # previous model
  ss<-readRDS('../f4_sob.rda') # previous validated sobol
  sens<-sobolBasis(mod,int.order = 2,mcmc.use = 1,verbose=F) # new sobol

  diff<-max(abs(range(ss$S.var-sens$S.var))) # difference between sobol indices
  expect_that(diff,is_less_than(eps))

  xtest<-readRDS('../f4_testX.rda') # x values
  pp<-predict(mod,xtest,trunc.error = F,nugget=F) # new predictions at x values
  pred<-readRDS('../f4_testPred.rda') # old predictions at x values
  diff<-max(abs(range(pred-pp))) # difference between new and old predictions
  expect_that(diff,is_less_than(eps))
})

