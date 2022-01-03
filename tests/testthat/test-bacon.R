## test bacon interface

library(hamstrbacon)

skip_on_cran()
skip_if_not_installed("rbacon", "2.5.2")

test_that("basic hamstr_bacon", {
  
  hb1 <- hamstr_bacon(id = "sdf", 
                      depth = MSB2K$depth,
                      obs_age = MSB2K$age,
                      obs_err = MSB2K$error, 
                      ssize = 10, burnin = 1)
  
  expect_equal(class(hb1), c("hamstr_bacon_fit", "list"))
  
  p <- plot(hb1)
  
  expect_equal(class(p), c("gg", "ggplot"))

  })



test_that("hamstr_bacon boundaries", {
  
  
  hb1 <- hamstr_bacon(id = "sdf", 
                      depth = MSB2K$depth,
                      obs_age = MSB2K$age,
                      obs_err = MSB2K$error, 
                      boundary = c(40),
                      ssize = 10, burnin = 1)
  
  expect_equal(class(hb1), c("hamstr_bacon_fit", "list"))
  
  expect_equal(length(hb1$pars$acc.mean), 2)
  expect_equal(length(hb1$info$acc.mean), 2)
  
  
  
})


test_that("hamstr_bacon hiatuses", {
  
  hb1 <- hamstr_bacon(id = "sdf", 
                      depth = MSB2K$depth,
                      obs_age = MSB2K$age,
                      obs_err = MSB2K$error, 
                      hiatus.depths = c(40),
                      ssize = 10, burnin = 1)
  
  expect_equal(class(hb1), c("hamstr_bacon_fit", "list"))
  
  expect_equal(length(hb1$pars$acc.mean), 2)
  expect_equal(length(hb1$info$acc.mean), 2)
  
  
})


test_that("hamstr_bacon marine20 and postbomb", {
  
  MSB2K.2 <- MSB2K
  
  MSB2K.2$age <- MSB2K.2$age - min(MSB2K.2$age) - 30
  
  hb1 <- hamstr_bacon(id = "sdf", 
                      depth = MSB2K.2$depth,
                      obs_age = MSB2K.2$age,
                      obs_err = MSB2K.2$error,
                      cc = 2, postbomb = 1,
                      ssize = 10, burnin = 1
                      )
  
  expect_equal(class(hb1), c("hamstr_bacon_fit", "list"))
  
  expect_equal(hb1$pars$postbomb, 1)
  
})


test_that("other options work: delta.R, delta.STD,
          t.a, t.b,
          d.min, d.max, th0", {
  
  hb2 <- hamstr_bacon(id = "sdf", 
                      depth = MSB2K$depth,
                      obs_age = MSB2K$age,
                      obs_err = MSB2K$error, 
                      th0 = 3000,
                      d.min = 9, d.max = 125,
                      thick = 5,
                      t.a = 33, t.b = 34,
                      delta.R = 100, delta.STD = 40,
                      ssize = 10, burnin = 1)
  
 
  expect_equal(class(hb2), c("hamstr_bacon_fit", "list"))
  
  expect_equal(hb2$pars$t.a, 33)
  expect_equal(hb2$info$t.a, 33)
  
  expect_equal(hb2$info$delta.R, 100)
  expect_equal(hb2$info$delta.STD, 40)
  
  fttd <- predict(hb2)
  
  expect_equal(min(fttd$depth), 9)
  expect_gt(max(fttd$depth), hb2$pars$d.max)
  
})
