
#conditions

lre_auto_bk <- function(A,x0){
  Asch <- Matrix::Schur(A)
  Asch2 <- QZ::qz.dtrsen(Asch$T, Asch$Q, abs(Asch$EValues) <= 1)
  n<-length(x0)
  Q1S <- Asch2$Q[1:n,1:n]
  Q2S <- Asch2$Q[(n+1):(nrow(Asch2$Q)),1:n]

  list(
    g = function(x0) Q2S %*% solve(Q1S) %*% x0,
    h = function(x0) A[1:n, 1:n] %*% x0
    + A[1:n, (n+1):ncol(A)] %*%
    (Q2S %*% solve(Q1S) %*% x0)
    )}

lre_auto_klein <- function(A, E, x0) {
  ret <- QZ::qz(A, E)
  ord <- abs(ret$ALPHA / ret$BETA) <= 1
  ret2 <- QZ::qz.dtgsen(ret$S, ret$T, ret$Q, ret$Z,select = ord)
  n <- length(x0)
  Z1S <- ret2$Z[1:n,1:n]
  Z2S <- ret2$Z[(n+1):(nrow(ret2$Z)),1:n]
  SSS <- ret2$T[1:n, 1:n]
  TSS <- ret2$S[1:n, 1:n]

  list(
    g = function(x0) Z2S %*% solve(Z1S) %*% x0,
    h = function(x0) Z1S %*% solve(SSS) %*% TSS %*% solve(Z1S) %*% x0
  )}

lre_auto <- function(A, E = NULL, x0) {
  if (is.null(E)) {
    lre_auto_bk(A, x0)
  } else {
    lre_auto_klein(A, E, x0)
  }
}


