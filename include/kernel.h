#ifndef __kernel__
#define __kernel__

#include "numerics.h"



long double Compute_square_amplitude_extended( complex<double> h1, complex<double> h2, complex<double> fA, complex<double> fV, complex<double> he1, complex<double> he2, complex<double> feA, complex<double> feV, double MK, double fk,double rl,double xk, double xq, double A, double B, double y, string MODE);

long double Compute_square_amplitude_different_lepton(complex<double> h1, complex<double> h2, complex<double> fA, complex<double> fV, double MK, double fk,double rl, double rll, double xk, double xq, double a, double b, double Y, string MODE );


#endif
