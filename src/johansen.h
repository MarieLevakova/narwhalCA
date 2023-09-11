// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef johansen_H
#define johansen_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::mat johansenCpp(arma::mat Z0, arma::mat Z1, int r, double dt, bool intercept);


// This is the end of the header guard
#endif
