// Perform Johansen procedure for alpha/beta-restricted VECM models.
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional... for now.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>     // Cpp linear algebra library
#include <stdio.h>
#include <math.h>
#include "johaTools.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat vecm(arma::mat Z0, arma::mat Z1, arma::mat Z2,
               int r, double dt, bool intercept){
  int N = Z0.n_rows;   // Number of observations
  int p = Z0.n_cols;   // Dimension of the system

  Z0 = Z0.t();  // transpose   // Easier to work with column vectors
  Z1 = Z1.t();

  mat R0,R1;
  mat M00,M11,M22,M01,M02,M12,M22_1;
  cx_mat N12;

  if(intercept){ // Check if Z2 has at least one column
    Z2 = Z2.t();

    // Find matrices M_{ij}
    M00 = M(Z0,Z0);  // crossproduct divided by N, operator defined in johaTools.h
    M11 = M(Z1,Z1);
    M22 = M(Z2,Z2);
    M01 = M(Z0,Z1);
    M02 = M(Z0,Z2);
    M12 = M(Z1,Z2);

    M22_1 = inv(M22);

    // Find residuals R0,R1
    R0 = Z0-M02*M22_1*Z2;
    R1 = Z1-M12*M22_1*Z2;
  } else {
    R0 = Z0;
    R1 = Z1;
  }

  // Find matrices S_{ij}
  mat S00,S11,S01,S10,S00_1; // Function S does the same operation as the function M
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);

  mat solveMat;
  mat cholS;

  // Solve for the eigenvalues (and vectors)

  N12 = powmat(S11, -0.5);
  solveMat = conv_to<mat>::from(N12)*S10*S00_1*S01*conv_to<mat>::from(N12);
  // solveMat = inv(S11)*S10*S00_1*S01;
  cholS = chol(S11,"lower");

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  // Normalize eigenvectors
  // for(int i=0;i<r;i++){
  //   double nf = as_scalar(sqrt(eigvec.col(i).t()*(S11)*eigvec.col(i)));
  //   eigvec.col(i) = eigvec.col(i)/sqrt(nf);
  // }

  vec testStat = -N*cumsum(log(1-eigval(sort_index(eigval,"ascend")))); // This is the standard test of H(r) versus H(p) in an unrestricted model.

  // To use cumsum for the teststats, the eigenvalues must be sorted "in reverse", this is ascending...
  testStat = sort(testStat,"descend");

  mat b_hat = conv_to<mat>::from(N12) * eigvec.cols(0,r-1); // Normalize to have beta' S_11 beta = I_r

  // Normalize beta using c {p x r} = (I{r x r},0{()p-r) x r})
  // mat Ir = eye<mat>(r,r);
  // mat c = join_cols(Ir,zeros<mat>(p-r,r));
  // b_hat = b_hat*inv(c.t()*b_hat);

  // Estimation of alpha and other parameters
  mat BS11B_1 = inv(b_hat.t()*S11*b_hat);

  mat a_hat = S01*b_hat*BS11B_1;
  mat Psi_hat;
  if(intercept){
    Psi_hat = (M02*M22_1-a_hat*b_hat.t()*M12*M22_1); // coefficients of lags and deterministic terms, but in fact only a constant is allowed
  } else {
    Psi_hat = zeros<mat>(p,1);
  }

  // Calculate residuals
  mat Pi_hat = a_hat*b_hat.t();
  mat res = zeros<mat>(p,N);
  mat Omega_hat = zeros<mat>(p,p);
  if(intercept){
    for(int n=0;n<N;n++){
      res.col(n) = Z0.col(n)-Pi_hat*(Z1.col(n))-Psi_hat*(Z2.col(n));
      Omega_hat = Omega_hat + (res.col(n)*(res.col(n)).t())/N;
    }
  } else {
    for(int n=0;n<N;n++){
      res.col(n) = Z0.col(n)-Pi_hat*(Z1.col(n));
      Omega_hat = Omega_hat + (res.col(n)*(res.col(n)).t())/N;
    }
  }

  int outRows = a_hat.n_cols+b_hat.n_cols+Psi_hat.n_cols+p+2+N;
  int outCols = p;
  mat out = zeros<mat>(outRows,outCols);

  // Insert alpha estimate in output
  for(int i=0; i<a_hat.n_cols;i++){
    out.row(i) = a_hat.col(i).t()/dt;
  }

  // Insert beta estimate in output
  int j = a_hat.n_cols;
  for(int i=0;i<b_hat.n_cols;i++){
    out.row(j+i) = b_hat.col(i).t();
  }
  // Insert Psi estimate in output
  j = a_hat.n_cols+b_hat.n_cols;
  for(int i=0;i<Psi_hat.n_cols;i++){
    out.row(j+i) = Psi_hat.col(i).t()/dt;
  }

  // Insert Omega estimate in output
  j = a_hat.n_cols+b_hat.n_cols+Psi_hat.n_cols;
  for(int i=0;i<p;i++){
    out.row(j+i) = Omega_hat.col(i).t()/dt;
  }

  // Insert test statistic in output
  j = a_hat.n_cols+b_hat.n_cols+Psi_hat.n_cols+p;
  out.row(j) = testStat.t();
  j = a_hat.n_cols+b_hat.n_cols+Psi_hat.n_cols+p+1;
  out.row(j) = eigval.t();

  // Insert residuals in output
  for(int i=0;i<N;i++){
    j = a_hat.n_cols+b_hat.n_cols+Psi_hat.n_cols+p+2;
    out.row(j+i) = res.col(i).t();
  }

  return out;
}


arma::mat var(arma::mat Z0, arma::mat Z1, arma::mat Z2, double dt, bool intercept){
 // Z2 is not incorporated in the code yet

  int N = Z0.n_rows;
  int p = Z0.n_cols;

  Z0 = Z0.t();
  Z1 = Z1.t();
  Z2 = Z2.t();

  mat Psi_hat = zeros<mat>(p,1);

  // Estimation using Least Squares
  // Formula: Psi_hat = T^-1*sum_{t=1}^T dX_t (this is 3x1 dimensional)
  if(intercept){
    for(int n=0;n<N;n++){
      Psi_hat += Z0.col(n);
    }
    Psi_hat = Psi_hat/N;
  }

  // Calculate residuals and covariance estimator (see Lütkepohls book, p. 75)
  mat res = zeros<mat>(p,N);
  mat Omega = zeros<mat>(p,p);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Psi_hat;
    Omega += res.col(n)*res.col(n).t();
  }
  Omega = (Omega/N)*(N/(N-1));


  // Calculate r-hypotheses statistics
  int r_tmp = 1;
  mat tmp = vecm(Z0.t(),Z1.t(), Z2.t(), r_tmp, dt, intercept);
  mat test  = tmp.rows(2*r_tmp+p+1,2*r_tmp+p+1);
  mat eigs  = tmp.rows(2*r_tmp+p+2,2*r_tmp+p+2);
  mat joha = join_cols(test,eigs);

  // Model estimates
  mat est = join_cols(Psi_hat.t()/dt, Omega/dt);

  // Add test statistics for r- hypotheses
  mat est_test = join_cols(est,joha);

  // Output estimates and residuals
  mat out  = join_cols(est_test,res.t());

  return out;
}


// Johansen estimation procedure, returns eigenvalues, eigenvectors
// [[Rcpp::export]]
arma::mat johansenCpp(arma::mat Z0, arma::mat Z1, arma::mat Z2, int r,
                                double dt, bool intercept){
  int N = Z0.n_rows;   // Number of observations
  int p = Z0.n_cols;     // Dimension of the system

  mat out = zeros<mat>(N,p);
  if(r > 0){
    out = vecm(Z0, Z1, Z2, r, dt, intercept);
  } else {
    out = var(Z0, Z1, Z2, dt, intercept);
  }
  return out;
}
