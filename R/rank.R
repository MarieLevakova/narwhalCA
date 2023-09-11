rank.johansen <- function(Y, Z, Z2, conf.level = 0.05, type = c("trace", "max")){
  # A function for estimating the cointegration rank by Johansen procedure,
  # namely the trace and the maximum eigenvalue test of the rank
  # Input: Y - matrix with differences
  #        Z - matrix with lagged values
  #        Z2 - matrix with lagged differences and other predictors
  #        conf.level - confidence level for the tests
  #        type - type of the test (trace test, maximum eigenvalue test)
  # Output: list with elements
  #         $r - one- or two-element vector with estimated rank(s)
  #         $lambdas - eigenvalues of the eigenvalue problem solved in the Johansen procedure
  #         $lrts - test statistics for testing all possible ranks

  alpha.levels <- c(0.1, 0.05, 0.01)
  P <- dim(Y)[2]  # no. of columns

  if (P > 10) print("Dimension exceed the allowed maximum of 10.") else {
    N <- dim(Y)[1]  # no. of rows

    cvals <- array(c(6.5, 12.91, 18.9, 24.78, 30.84, 36.25,
                     42.06, 48.43, 54.01, 59, 65.07, 8.18, 14.9, 21.07, 27.14,
                     33.32, 39.43, 44.91, 51.07, 57, 62.42, 68.27, 11.65,
                     19.19, 25.75, 32.14, 38.78, 44.59, 51.3, 57.07, 63.37,
                     68.61, 74.36, 6.5, 15.66, 28.71, 45.23, 66.49, 85.18,
                     118.99, 151.38, 186.54, 226.34, 269.53, 8.18, 17.95,
                     31.52, 48.28, 70.6, 90.39, 124.25, 157.11, 192.84, 232.49,
                     277.39, 11.65, 23.52, 37.22, 55.43, 78.87, 104.2, 136.06,
                     168.92, 204.79, 246.27, 292.65), c(11, 3, 2))

    lambdas <- vecm(Y, Z, Z2, r = 1, dt = 1, intercept = T)[P+5,]
    lrts1 <- -N*rev(cumsum(log(1-rev(lambdas))))
    lrts2 <- -N*log(1-lambdas)

    r <- c(NA, NA)
    names(r) <- c("trace", "max")
    r["trace"] <- which(c(lrts1, 0) <= c(cvals[0:P, conf.level == alpha.levels, 2], 1))[1] - 1
    r["max"] <- which(c(lrts2, 0) <= c(cvals[0:P, conf.level == alpha.levels, 1], 1))[1] - 1
    return(list(r = r[type],
                lambdas = lambdas,
                lrts = data.frame(trace = lrts1, max = lrts2,
                                  row.names = paste("r =", 0:(P-1)))))
  }
}
