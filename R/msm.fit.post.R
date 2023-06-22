


msm.fit.post = function(msm.fit.object, mod.list, suStf){


  pVbres = msm.postVb(msm.fit.object)
  He <- pVbres$He
  Vb <- pVbres$Vb
  HeSh <- pVbres$HeSh
  F <- pVbres$F
  F1 <- pVbres$F1
  R <- pVbres$R
  Ve <- pVbres$Ve
  t.edf <- pVbres$t.edf
  msm.fit.object <- pVbres$msm.fit.object


  msm.post.object = list()


  for(i in 1:sum(suStf$whereQ != 0)){

    msm.post.object[[i]] = mod.list[[i]]
    msm.post.object[[i]]$coefficients <- msm.fit.object$fit$argument[suStf$pos.optparams][suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)] # first indexing is to augment back to the full set of values (i.e. repeating the constrained parameters, same for below); the second is to make the transition specific selection
    msm.post.object[[i]]$coef.names = attr(suStf$full.X, 'dimnames')[[2]][suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)]
    msm.post.object[[i]]$Vp <- Vb[suStf$pos.optparams,suStf$pos.optparams][suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1), suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)]
    msm.post.object[[i]]$Ve <- Ve[suStf$pos.optparams,suStf$pos.optparams][suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1), suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)]
    msm.post.object[[i]]$edf <- diag(F[suStf$pos.optparams,suStf$pos.optparams])[suStf$start.pos.par[i]:(suStf$start.pos.par[i+1]-1)]

  }





  list(msm.post.object = msm.post.object,
       He = He, Vb = Vb, HeSh = HeSh, F = F,
       F1 = F1, R = R, Ve = Ve, t.edf = t.edf,
       logLik = -msm.fit.object$fit$l)

}
