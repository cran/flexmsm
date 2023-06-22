
test_that("Model setup (without model fitting on simulated data).", {


  IDM_sim <- simulateIDM(N = 100, seed = 42, og.12 = TRUE)

  formula = list(years ~ s(years, bs = 'cr', k = 10), # 1-2
                 years ~ s(years, bs = 'cr', k = 10), # 1-3
                 0,                                   # 2-1
                 years ~ s(years, bs = 'cr', k = 10), # 2-3
                 0,                                   # 3-1
                 0                                    # 3-2
  )


  msm.out = fmsm(formula = formula, data = IDM_sim, id = PTNUM, state = state, death = TRUE,
                fit = FALSE)

  expect_equal(round(sum(msm.out$suStf$full.S), 3), round(1.778928, 3))
  expect_equal(round(sum(msm.out$suStf$full.X), 3), round(2898,3))

})




test_that("Model setup (without model fitting).", {


  IDM_cav <- read.csv(test_path("testing_data", "IDM_cav.csv"))

  formula = list(years ~ s(years, bs = 'cr', k = 10) + s(dage, bs = 'cr', k = 10) + pdiag, # 1-2
                 years ~ s(years, bs = 'cr', k = 10) + dage + pdiag, # 1-3
                 0,                                                  # 2-1
                 years ~ s(years, bs = 'cr', k = 10) + s(dage, bs = 'cr', k = 10) + ti(years, dage, bs = 'cr') + pdiag, # 2-3
                 0,                                  # 3-1
                 0                                   # 3-2
  )

  msm.out = fmsm(formula = formula, data = IDM_cav, id = PTNUM, state = state, death = TRUE,
                fit = FALSE)

  expect_equal(round(sum(msm.out$suStf$full.S), 3), round(17.15939, 3))
  expect_equal(round(sum(msm.out$suStf$full.X), 3), round(6687.129, 3))

})




test_that("Q matrix setup check.", {


  IDM_cav <- read.csv(test_path("testing_data", "IDM_cav.csv"))

  formula = list(years ~ s(years, bs = 'cr', k = 10) + s(dage, bs = 'cr', k = 10) + pdiag, # 1-2
                 years ~ s(years, bs = 'cr', k = 10) + dage + pdiag, # 1-3
                 0,                                                  # 2-1
                 years ~ s(years, bs = 'cr', k = 10) + s(dage, bs = 'cr', k = 10) + ti(years, dage, bs = 'cr') + pdiag, # 2-3
                 0,                                  # 3-1
                 0                                   # 3-2
  )

  msm.out = fmsm(formula = formula, data = IDM_cav, id = PTNUM, state = state, death = TRUE,
                fit = FALSE)

  Qmatr = Q.matr.setup.general(params = msm.out$suStf$params, nstates = msm.out$suStf$nstates,
                               full.X = msm.out$suStf$full.X, start.pos.par = msm.out$suStf$start.pos.par,
                               l.short.formula = msm.out$suStf$l.short.formula, whereQ = msm.out$suStf$whereQ,
                               firstD = FALSE, secondD = FALSE, bound.eta = FALSE,
                               pos.optparams = msm.out$suStf$pos.optparams, pos.optparams2 = msm.out$suStf$pos.optparams2)$Qmatr


  expect_equal(round(sum(Qmatr), 3), round(-1.769418e-15, 3))




})



