estandarizar <- function(X, kte=10, centrarfila=FALSE)
{
  if (centrarfila)
  {
    m <- matrix(apply(X, 1, mean), ncol=1)
    uno <- rep(1, dim(X)[2])
    X <- X - m%*%t(uno)
  }
  aux <- sqrt(sum(diag(X%*%t(X))))
  X <- X/aux*sqrt(kte)
  return(X)
}

norm1 <- function(W)
{
  # Input: matrix of eigenvectors
  # Output: eigenvector with norm 1
  apply(W, 2, function(x){x/sqrt(sum(x^2))})
}

doublecenter <- function(D)
{
  # Computes the doble centering of a distance matrix
  D <- as.matrix(D)
  n <- dim(D)[1]
  uno <- matrix(1, nrow=n, ncol=1)
  H <- diag(1, n) - 1/n*uno%*%t(uno)
  A <- -1/2*D^2
  B <- H%*%A%*%H
  return(B)
}

Bd <- function(D, X){
  # Distance Based B is computed.
  # If B not definite positive, then it is converted
  D <- as.matrix(D)
  ch <- dim(D)[1]
  m <- apply(X, 2, mean)
  m <- as.matrix(m, ncol=1)
  uno <- matrix(rep(1, ch), ncol=1)
  Bi <- doublecenter(D) + X%*%m%*%t(uno) + uno%*%t(m)%*%t(X) - uno%*%t(m)%*%m%*%t(uno)
  return(Bi)
}

calcVarianzas <- function(X, W)
{
  # Input: X, i-esimo EEG
  # Output: varianza de la curva proyectada en las direcciones W
  diag(t(W)%*%X%*%t(X)%*%W)
}

Mixture.dist <- function(d, w, eps=10^-10)
{
  #----- To calculate Related MEtric Scaling
  # Input:
  #     d: a list with k=2 distance matrices
  #     w: weight for mixture
  #     eps: accuracy for numerical 0. Values in (-eps, eps) are considered as 0
  # Output:
  #     dmixt: Mixture distance --> w D1 + (1-w)D2/max{D2}*max{D1}
  #            Tomamos como referencia la distancia euclidea, para luego, al rehacer la B, podamos añadir la X sin problema de unidades
  #-----------------------------------------

  k <- length(d)
  if(k != 2) stop("Mixture distance must be calculated on 2 distances\n")

  M <- unlist(lapply(d, max))
  dmixt <- w*d[[1]] + (1-w)*d[[2]]/M[2]*M[1]
  return(dmixt)
}
