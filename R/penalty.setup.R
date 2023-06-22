


penalty.setup = function(sp, suStf){

  # Quantities needed
  S.list = suStf$S.list
  start.pos.par.detailed = suStf$start.pos.par.detailed
  start.pos.par.only.smooth = suStf$start.pos.par.only.smooth
  start.pos.par.only.smooth.FPC = suStf$start.pos.par.only.smooth.FPC

  block.dims = diff(start.pos.par.detailed)
  pen.matr.S.lambda = matrix(0, ncol = block.dims[1], nrow = block.dims[1]) # this is for the non-smooth related portions, so lambda is just set to 0 (intercept in particular)

  sp.mapping = match(start.pos.par.detailed, start.pos.par.only.smooth)

  sp.mapping.FPC = list()
  j = 1
  for(idx in start.pos.par.detailed) {if(idx %in% start.pos.par.only.smooth.FPC) {sp.mapping.FPC[[j]] = which(idx == start.pos.par.only.smooth.FPC); j = j + 1}}
  rm(j)


  for(i in 2:length(block.dims)){

    if(is.na(sp.mapping[i])){
      pen.matr.S.lambda = cbind( rbind(pen.matr.S.lambda,
                                   matrix(0, ncol = ncol(pen.matr.S.lambda), nrow = block.dims[i])),
                             rbind(matrix(0, ncol = block.dims[i], nrow = nrow(pen.matr.S.lambda)),
                                   matrix(0, ncol = block.dims[i], nrow = block.dims[i])))
    } else {
      ind = sp.mapping.FPC[[sp.mapping[i]]]
      S.k.lambda = mapply('*', sp[ind] , S.list[ind], SIMPLIFY = FALSE)
      S.k.lambda = Reduce("+", S.k.lambda)
      pen.matr.S.lambda = cbind( rbind(pen.matr.S.lambda,
                                   matrix(0, ncol = ncol(pen.matr.S.lambda), nrow = block.dims[i])),
                             rbind(matrix(0, ncol = block.dims[i], nrow = nrow(pen.matr.S.lambda)),
                                   matrix(S.k.lambda, ncol = block.dims[i], nrow = block.dims[i])))
    }
  }



  pen.matr.S.lambda
}
