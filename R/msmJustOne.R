




msmJustOne = function(params.0, sp, pen.matr.S.lambda, pmethod, suStf, death,
                   Q.diagnostics = FALSE, parallel, no_cores, justComp = NULL){

  Qmatr.diagnostics.list = gradient = hessian = NULL

  # Unpack necessary variables -------------------------------------------------------------------
  data = suStf$data # HERE WE USED TO HAVE $data BUT I THINK THAT IN ALL OTHER TESTS I NORMALLY DO data AND data.long COINCIDE (not in no intercept test though)
  nstates = suStf$nstates
  start.pos.par = suStf$start.pos.par

  pos.optparams = suStf$pos.optparams
  pos.optparams2 = suStf$pos.optparams2

  l.short.formula = suStf$l.short.formula
  whereQ = suStf$whereQ
  full.X = suStf$full.X

  Sl.sf = suStf$Sl.sf # what is the difference between Sl.sf and S.list ? (check this - it seems former is only for EFS)


  start.pos.par.only.smooth = suStf$start.pos.par.only.smooth
  start.pos.par.only.smooth.FPC = suStf$start.pos.par.only.smooth.FPC
  start.pos.par.detailed = suStf$start.pos.par.detailed

  MM = list(start.pos.par = start.pos.par,
            pos.optparams = pos.optparams,
            pos.optparams2 = pos.optparams2,
            l.short.formula = l.short.formula,
            whereQ = whereQ,
            nstates = nstates,
            l.params = length(params.0),
            cens.state = suStf$cens.state)


  # ******************** #
  # Penalty matrix setup #
  # ******************** #
  # Setup full penalty matrix to be used for penalized likelihood estimation
  pen.matr.S.lambda = penalty.setup(sp = sp, suStf = suStf)
  pen.matr.S.lambda = pen.matr.S.lambda[1:max(pos.optparams2) == pos.optparams2, 1:max(pos.optparams2) == pos.optparams2]
  # ----------------------------------------------------------------------------------------------


  do.grad.inner = ifelse('grad' %in% justComp, TRUE, FALSE)
  do.hess.inner = ifelse('hess' %in% justComp, TRUE, FALSE)

  singleComp <- LikGradHess.CM(params.0, data = data, full.X = full.X, MM = MM, pen.matr.S.lambda = pen.matr.S.lambda, # *****
                                    aggregated.provided = FALSE, do.gradient = do.grad.inner, do.hessian = do.hess.inner,
                                    pmethod = pmethod, death = death,
                                    Qmatr.diagnostics.list = Qmatr.diagnostics.list,
                                    parallel = parallel, no_cores = no_cores,
                               P.save.all = TRUE, CM.comp = TRUE)



  if(do.grad.inner) gradient = singleComp$gradient
  if(do.hess.inner) hessian = singleComp$hessian

  list(value = singleComp$value, gradient = gradient, hessian = hessian,
       argument = params.0,
       S.h = singleComp$S.h, S.h1 = singleComp$S.h1, S.h2 = singleComp$S.h2,
       apprHessian = singleComp$apprHessian,
       l = singleComp$l,
       Qmatr.diagnostics.list = singleComp$Qmatr.diagnostics.list,
       sp = sp,
       P.hist = singleComp$P.hist,
       dP.hist = singleComp$dP.hist,
       d2P.hist = singleComp$d2P.hist,
       ind.CM = singleComp$ind.CM,
       P.CM.contr = singleComp$P.CM.contr)




}
