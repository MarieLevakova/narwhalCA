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
arma::mat glsCpp(arma::mat Y, arma::mat Z, arma::mat Z2,
               arma::mat H, arma::mat h,
               arma::mat A0, arma::mat B0, arma::mat C0,
               arma::cube Omega0, arma::vec break_pts,
               int r, int n_iter){

  int N = Y.n_rows;   // Number of observations
  int p = Y.n_cols;
  int p2 = Z2.n_cols;
  int m = break_pts.n_elem + 1;

  vec breaks_extended = zeros<vec>(m+1);
  breaks_extended.subvec(1,m-1) = break_pts;
  breaks_extended(m) = N+1;

  vec Lmax = zeros<vec>(n_iter);
  mat res_all = zeros<mat>(p, N);

  for(int it=0; it<n_iter; it++){
    // Update A and C
    mat sum_mat1 = zeros<mat>(p*(m*r+p2), p*(m*r+p2));
    mat sum_mat2 = zeros<mat>(p*(m*r+p2), 1);
    for(int mm=1; mm<m+1; mm++){
      for(int j=breaks_extended(mm-1); j<breaks_extended(mm)-1; j++){
        mat blok11 = B0.t() * (Z.row(j)).t() * Z.row(j) * B0;
        mat blok21 = (Z2.row(j)).t() * Z.row(j) * B0;
        mat blok22 = (Z2.row(j)).t() * Z2.row(j);

        mat aux_mat1 = join_cols(join_rows(blok11, blok21.t()),
                                 join_rows(blok21, blok22));
        sum_mat1 = sum_mat1 + kron(aux_mat1, (Omega0.slice(mm-1)).i());
        sum_mat2 = sum_mat2 +
          vectorise((Omega0.slice(mm-1)).i() * (Y.row(j)).t() * join_rows(Z.row(j)*B0, Z2.row(j)));
      }
    }

    mat vec_AC = sum_mat1.i() * sum_mat2;
    mat AC = reshape(vec_AC, 3, m*r+p2);

    A0 = AC.cols(0,r*m-1);
    C0 = AC.cols(r*m,r*m+p2-1);

    // Update B
    mat sum_mat3 = zeros<mat>(m*m*r*p, m*m*r*p);
    mat sum_mat4 = zeros<mat>(m*m*r*p, 1);
    for(int mm=1; mm<m+1; mm++){
      for(int j=breaks_extended(mm-1); j<breaks_extended(mm)-1; j++){
        sum_mat3 = sum_mat3 +
          kron(A0.t() * (Omega0.slice(mm-1)).i() * A0, (Z.row(j)).t() * Z.row(j));
        sum_mat4 = sum_mat4 +
          vectorise((Z.row(j)).t() * ((Y.row(j)).t() - C0 * (Z2.row(j)).t()).t() * (Omega0.slice(mm-1)).i() * A0);
      }
    }

    mat vec_B = H * (H.t() * sum_mat3 * H).i() * H.t() * (sum_mat4 - sum_mat3 * h) + h;
    B0 = reshape(vec_B, m*p, m*r);

    // Update Omega

    cube Omega_aux = zeros<cube>(p,p,m);
    for(int mm=1; mm<m+1; mm++){
      // Omega0.slice(mm-1) = Omega0.slice(mm-1)-Omega0.slice(mm-1);
      mat res = zeros<mat>(p, breaks_extended(mm)-1-breaks_extended(mm-1));
      for(int j=breaks_extended(mm-1); j<breaks_extended(mm)-1; j++){
        res_all.col(j) = (Y.row(j)).t() - A0 * B0.t() * (Z.row(j)).t() - C0 * (Z2.row(j)).t();
        Omega_aux.slice(mm-1) = Omega_aux.slice(mm-1) + (res_all.col(j) * (res_all.col(j)).t())/(breaks_extended(mm)-1-breaks_extended(mm-1));
      }
    }

    Omega0 = Omega_aux;

    Lmax(it) = -N*3/2*log(2*datum::pi) - N*3/2;
    for(int mm=1; mm<m+1; mm++){
      double lOmega = log_det_sympd(Omega0.slice(mm-1));
      Lmax(it) = Lmax(it) -
        (breaks_extended(mm)-1-breaks_extended(mm-1))*lOmega/2;
    }
  }

  mat out = zeros<mat>(m*(2*r+p)+p2+N+n_iter,p);

  // Insert A estimate in output
  for(int i=0; i<A0.n_cols;i++){
    out.row(i) = A0.col(i).t();
  }

  // Insert B estimate in output
  int i = A0.n_cols;
  for(int mm=0; mm<m; mm++){
    mat B0_mm = B0.submat(mm*p, mm*r, (mm+1)*p-1, (mm+1)*r-1);
    for(int rr=0; rr<r; rr++){
      out.row(i) = B0_mm.col(rr).t();
      i = i+1;
    }
  }

  // Insert C estimate in output
  for(int j=0; j<p2;j++){
    out.row(i+j) = C0.col(j).t();
  }
  i = i + p2;

  // Insert Omega estimate in output
  for(int mm=0; mm<m; mm++){
    mat Omega_mm = Omega0.slice(mm);
    for(int rr=0; rr<p;rr++){
      out.row(i) = Omega_mm.col(rr).t();
      i = i+1;
    }
  }

  // Insert the residuals in output
  for(int j=0; j<N;j++){
    out.row(i+j) = res_all.col(j).t();
  }

  // Insert the values of the log-likelihood
  i = i+N;
  for(int j=0; j<n_iter; j++){
    out(i+j,0) = Lmax(j);
  }

  return out;

}

// [[Rcpp::export]]
arma::mat glsConstCpp(arma::mat Y, arma::mat Z, arma::mat Z2,
                 arma::mat H, arma::mat h,
                 arma::mat A0, arma::mat B0, arma::mat C0,
                 arma::mat Omega0, arma::vec break_pts,
                 int r, int n_iter){

  int N = Y.n_rows;   // Number of observations
  int p = Y.n_cols;
  int p2 = Z2.n_cols;
  int m = break_pts.n_elem + 1;
  int MM = B0.n_cols;
  int MMM = B0.n_rows;
  int p1 = MMM/MM;
  int K = H.n_cols/(p-1); // Number of cointegration relationships

  // Calculate moment matrices
  mat M00 = M(Y.t(), Y.t());
  mat M01 = M(Y.t(), Z.t());
  mat M02 = M(Y.t(), Z2.t());
  mat M10 = M(Z.t(), Y.t());
  mat M11 = M(Z.t(), Z.t());
  mat M12 = M(Z.t(), Z2.t());
  mat M20 = M(Z2.t(), Y.t());
  mat M21 = M(Z2.t(), Z.t());
  mat M22 = M(Z2.t(), Z2.t());

  mat M22_1 = inv(M22);

  // Find residuals R0,R1
  mat R0 = Y-Z2*M22_1*M20;
  mat R1 = Z-Z2*M22_1*M21;

  // Find moment matrices of the residuals
  mat S00 = S(R0.t(), R0.t());
  mat S01 = S(R0.t(), R1.t());
  mat S10 = S(R1.t(), R0.t());
  mat S11 = S(R1.t(), R1.t());

  // vec breaks_extended = zeros<vec>(m+1);
  // breaks_extended.subvec(1,m-1) = break_pts;
  // breaks_extended(m) = N+1;
  //
  vec Lmax = zeros<vec>(n_iter);
  mat res_all = zeros<mat>(p, N);
  //
  for(int it=0; it<n_iter; it++){

    // Update A

    A0 = S01*B0*(B0.t()*S11*B0).i();

    // Update C
    C0 = M02*M22_1-A0*B0.t()*M12*M22_1;

    // Update B

    mat vec_B = H * (H.t() * kron(A0.t()*(inv(Omega0))*A0, S11) * H).i() * H.t() * (kron(A0.t(),S10)*vectorise(inv(Omega0)) -
      kron(A0.t()*inv(Omega0)*A0,S11) * h) + h;

    B0 = reshape(vec_B, MMM, MM);

    // Update Omega

    Omega0 = S00 - S01*B0*(B0.t()*S11*B0).i()*B0.t()*S10;
    for(int j=1; j<N; j++){
      res_all.col(j) = (Y.row(j)).t() - A0 * B0.t() * (Z.row(j)).t() - C0 * (Z2.row(j)).t();
    }

    Lmax(it) = -N*p*(1+log(2*datum::pi))/2 -N*log_det_sympd(Omega0)/2;
  }

  mat out = zeros<mat>(MM+MM+p+p2+N+n_iter,p1);

  int i;

  // Insert A estimate in output
  for(int i=0; i<A0.n_cols;i++){
    out.submat(i,0,i,p-1) = A0.col(i).t();
  }

  // Insert B estimate in output
  i = A0.n_cols;

  for(int mm=0; mm<MM; mm++){
    mat B0_mm = B0.col(mm);
    // Pick up nonzero values
    mat B0_mm_nonzero = zeros<mat>(p1,1);
    int i_B0 = 0;
    for(int iii=0; iii<m*p1; iii++){
      if(B0_mm(iii,0)!=0){
        B0_mm_nonzero(i_B0,0) = B0_mm(iii,0);
        i_B0 = i_B0+1;
      }
    }
    out.row(i) = B0_mm_nonzero.t();
    i = i+1;
  }

  // Insert C estimate in output
  for(int j=0; j<p2;j++){
    out.submat(i+j,0,i+j,p-1) = (C0.col(j)).t();
  }
  i = i + p2;

  // Insert Omega estimate in output
  for(int rr=0; rr<p;rr++){
    out.submat(i,0,i,p-1) = Omega0.col(rr).t();
    i = i+1;
  }

  // Insert the residuals in output
  for(int j=0; j<N;j++){
    out.submat(i+j,0,i+j,p-1) = res_all.col(j).t();
  }

  // Insert the values of the log-likelihood
  i = i+N;
  for(int j=0; j<n_iter; j++){
    out(i+j,0) = Lmax(j);
  }

  return out;

}
