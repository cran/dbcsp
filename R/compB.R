calculateDistance <- function(X, type){

  distances_TSdist <- c("infnorm", "ccor", "sts", "lb.keogh", "edr", "erp", "lcss", "fourier", "tquest", "dissim", "acf",
                        "pacf", "ar.lpc.ceps", "ar.mah", "ar.mah.statistic", "ar.mah.pvalue", "ar.pic", "cdm", "cid", "cor",
                        "cort", "int.per", "per", "mindist.sax", "ncd", "pred",  "spec.glk", "spec.isd", "spec.llr", "pdc",
                        "frechet", "tam")

  distances_parallelDist <- c("bhjattacharyya", "bray", "canberra", "chord", "divergence", "dtw", "euclidean", "fJaccard",
                              "geodesic", "hellinger", "kullback", "mahalanobis", "manhattan", "maximum", "minkowski",
                              "podani", "soergel", "wave", "whittaker") #LAS DE BINARY INPUT VARIABLES NO LAS HE PUESTO

  if(type==toupper(type)) type <- tolower(type) #si está en mayus pasar a minúscula

  if(type %in% distances_parallelDist) dist_type <- 'parDist'
  else if(type %in% distances_TSdist) dist_type <- 'TSdist'
  else stop('Invalid distance type. See ?dbcsp for more information.')

  Ds <- switch(
    dist_type,
    'TSdist' = Ds <- plyr::llply(X, TSdist::TSDatabaseDistances, distance=type),
    'parDist' = Ds <- plyr::llply(X, parallelDist::parDist, method=type)
  )
  return(Ds)
}

compB <- function(X, mixture, type, w, eig.tol = 1e-06, getWarning=TRUE)
{
  # Put together all de Bi-s
  #
  # EUCL
  if(type=='EUCL' || type=='euclidean') Bs <- plyr::llply(X, function(x){x%*%t(x)})
  else # DB and MIX
  {
    if(is.list(X)) n <- length(X)
    else n <- 1
    # DB
    Ds <- calculateDistance(X,type)
    # MIX
    if(mixture){
      Ds1 <- calculateDistance(X,'euclidean')
      Ds2 <- Ds
      Ds <- vector("list", n)
      for (r in 1:n)
      {
        aux <- list(Ds1[[r]], Ds2[[r]])
        Ds[[r]] <- Mixture.dist(aux, w=w)
      }
    }
    Bs <- vector("list", n)
    for (r in 1:n) Bs[[r]] <- Bd(D=Ds[[r]], X=X[[r]])
  }

  # pooled B matrix
  Bs <- plyr::laply(Bs, as.matrix)
  B <- plyr::aaply(Bs, c(2, 3), mean)

  # DB and MIX
  if(type!='EUCL'){
    # Check whether it is a positive definite matrix
    vp <- eigen(B)$values
    if (min(vp) < -eig.tol)
    {
      B <- Matrix::nearPD(B)$mat
      B <- as.matrix(B)
      if(getWarning) warning('Distance matrix was converted to be definite positive',immediate. = TRUE)
      #warning("B was converted to be definite positive")
    }
  }

  return(B)
}
