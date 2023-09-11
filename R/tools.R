# Various auxiliary functions

matrix.angle <- function(U, V){
  acos(sum(diag(t(U) %*% V))/(sqrt(sum(diag(t(U) %*% U))) * sqrt(sum(diag(t(V) %*% V)))))
}
