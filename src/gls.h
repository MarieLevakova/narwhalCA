// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef gls_H
#define gls_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::mat glsCpp(arma::mat Y, arma::mat Z, arma::mat Z2, arma::mat H, arma::mat h,
                 arma::mat A0, arma::mat B0, arma::mat C0,
                 arma::cube Omega0, arma::vec break_pts,
                 int r, int n_iter);


// This is the end of the header guard
#endif
