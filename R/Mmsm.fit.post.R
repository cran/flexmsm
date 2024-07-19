


Mmsm.fit.post <- function(Mmsm.fit.object, l.params1, l.params2){


  suStf1 = Mmsm.fit.object$fit$proc1$suStf
  mod.list1 = suStf1$mod.list

  suStf2 = Mmsm.fit.object$fit$proc2$suStf
  mod.list2 = suStf2$mod.list


  pVbres <- Mmsm.postVb(Mmsm.fit.object)
  He <- pVbres$He
  Vb <- pVbres$Vb
  HeSh <- pVbres$HeSh
  F <- pVbres$F
  F1 <- pVbres$F1
  R <- pVbres$R
  Ve <- pVbres$Ve
  t.edf <- pVbres$t.edf
  Mmsm.fit.object <- pVbres$Mmsm.fit.object


  Mmsm.post.object = list()


  # PROCESS 1
  ii = 1
  t.edf1 = 0
  for(i in 1:sum(suStf1$whereQ != 0)){

    Mmsm.post.object[[ii]] = mod.list1[[i]]
    Mmsm.post.object[[ii]]$coefficients <- Mmsm.fit.object$fit$argument[suStf1$pos.optparams][suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)] # first indexing is to augment back to the full set of values (i.e. repeating the constrained parameters, same for below); the second is to make the transition specific selection
    Mmsm.post.object[[ii]]$coef.names = attr(suStf1$full.X, 'dimnames')[[2]][suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)]
    Mmsm.post.object[[ii]]$Vp <- Vb[suStf1$pos.optparams,suStf1$pos.optparams][suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1), suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)]
    Mmsm.post.object[[ii]]$Ve <- Ve[suStf1$pos.optparams,suStf1$pos.optparams][suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1), suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)]
    Mmsm.post.object[[ii]]$edf <- diag(F[suStf1$pos.optparams,suStf1$pos.optparams])[suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)]

    t.edf1 = t.edf1 + sum(diag(F[suStf1$pos.optparams,suStf1$pos.optparams])[suStf1$start.pos.par[i]:(suStf1$start.pos.par[i+1]-1)])

    ii = ii + 1
  }

  # PROCESS 2 -- MAKE INDEXING ADJUSTMENTS FOR PROCESS 2
  dbg.chk = 0
  t.edf2 = 0
  for(i in 1:sum(suStf2$whereQ != 0)){

    Mmsm.post.object[[ii]] = mod.list2[[i]]
    Mmsm.post.object[[ii]]$coefficients <- Mmsm.fit.object$fit$argument[l.params1+suStf2$pos.optparams][suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)] # first indexing is to augment back to the full set of values (i.e. repeating the constrained parameters, same for below); the second is to make the transition specific selection
    Mmsm.post.object[[ii]]$coef.names = attr(suStf2$full.X, 'dimnames')[[2]][suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)]
    Mmsm.post.object[[ii]]$Vp <- Vb[l.params1+suStf2$pos.optparams,l.params1+suStf2$pos.optparams][suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1), suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)]
    Mmsm.post.object[[ii]]$Ve <- Ve[l.params1+suStf2$pos.optparams,l.params1+suStf2$pos.optparams][suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1), suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)]
    Mmsm.post.object[[ii]]$edf <- diag(F[l.params1+suStf2$pos.optparams,l.params1+suStf2$pos.optparams])[suStf2$start.pos.par[i]:(suStf2$start.pos.par[i+1]-1)]

    dbg.chk = dbg.chk + length(Mmsm.post.object[[ii]]$coefficients)

    t.edf2 = t.edf2 + sum(Mmsm.post.object[[ii]]$edf)

    ii = ii + 1
  }

  if(dbg.chk != l.params2) stop('There is something wrong with the post-fit object setup.') # only for debug purposes - remove later


  list(Mmsm.post.object = Mmsm.post.object,
       He = He, Vb = Vb, HeSh = HeSh, F = F,
       F1 = F1, R = R, Ve = Ve, t.edf = t.edf,
       t.edf1 = t.edf1, t.edf2 = t.edf2,
       logLik = -Mmsm.fit.object$fit$l)




}
