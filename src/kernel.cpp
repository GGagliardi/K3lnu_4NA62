#include "../include/kernel.h"
using namespace std;
const double LeviCivitaSign=-1;

//variables:
//x is normalized square root of photon invariant mass
//y is normalized square root of W invariant mass
//rl and rl2 are the two normalized lepton masses
//mk is the meson mass







long double Compute_square_amplitude_different_lepton(complex<double> h1, complex<double> h2, complex<double> fA, complex<double> fV, double MK, double fk,double rl, double rll, double xk, double xq, double A, double B, double y, string MODE ) {


   
  h1 /= pow(xk,2);
  h2 /= pow(xk,2)*(1.0-pow(xq,2));
  fA /= pow(xk,2);
  fV /= LeviCivitaSign*pow(xk,2);


  double pt = (8*pow(MK,4)*rl*(-2*rll*pow(xq,2) + pow(xk,2)*pow(xq,2) - 4*B*pow(xk,2)*pow(xq,2) + 2*pow(B,2)*pow(xk,2)*pow(xq,2) - pow(xk,4)*pow(xq,2) + 
       5*B*pow(xk,4)*pow(xq,2) - 2*pow(B,2)*pow(xk,4)*pow(xq,2) - 2*B*pow(xk,6)*pow(xq,2) - 2*pow(xq,4) + 4*B*pow(xq,4) - 
       2*pow(B,2)*pow(xq,4) + 4*rll*pow(xq,4) + 3*pow(xk,2)*pow(xq,4) - 5*B*pow(xk,2)*pow(xq,4) + 2*pow(B,2)*pow(xk,2)*pow(xq,4) - 
       pow(xk,4)*pow(xq,4) + 3*B*pow(xk,4)*pow(xq,4) - 2*rll*pow(xq,6) - B*pow(xk,2)*pow(xq,6) + 2*rll*y - pow(xk,2)*y + 4*B*pow(xk,2)*y - 
       2*pow(B,2)*pow(xk,2)*y - 2*rll*pow(xk,2)*y + 2*pow(xk,4)*y - 3*B*pow(xk,4)*y - pow(xk,6)*y + 4*pow(xq,2)*y - 8*B*pow(xq,2)*y + 
       4*pow(B,2)*pow(xq,2)*y - 2*rll*pow(xq,2)*y - 9*pow(xk,2)*pow(xq,2)*y + 12*B*pow(xk,2)*pow(xq,2)*y - 2*pow(B,2)*pow(xk,2)*pow(xq,2)*y + 
       4*rll*pow(xk,2)*pow(xq,2)*y + 5*pow(xk,4)*pow(xq,2)*y - 5*B*pow(xk,4)*pow(xq,2)*y - pow(xk,6)*pow(xq,2)*y + 4*pow(xq,4)*y - 
       4*B*pow(xq,4)*y - 2*rll*pow(xq,4)*y - 4*pow(xk,2)*pow(xq,4)*y + 4*B*pow(xk,2)*pow(xq,4)*y - 2*rll*pow(xk,2)*pow(xq,4)*y + 
       pow(xk,4)*pow(xq,4)*y + 2*rll*pow(xq,6)*y - 2*pow(y,2) + 4*B*pow(y,2) - 2*pow(B,2)*pow(y,2) - 2*rll*pow(y,2) + 6*pow(xk,2)*pow(y,2) - 
       7*B*pow(xk,2)*pow(y,2) - 4*pow(xk,4)*pow(y,2) - 8*pow(xq,2)*pow(y,2) + 8*B*pow(xq,2)*pow(y,2) + 4*rll*pow(xq,2)*pow(y,2) + 
       9*pow(xk,2)*pow(xq,2)*pow(y,2) - 3*B*pow(xk,2)*pow(xq,2)*pow(y,2) - 2*pow(xk,4)*pow(xq,2)*pow(y,2) - 2*pow(xq,4)*pow(y,2) - 
       2*rll*pow(xq,4)*pow(y,2) + pow(xk,2)*pow(xq,4)*pow(y,2) + 4*pow(y,3) - 4*B*pow(y,3) - 5*pow(xk,2)*pow(y,3) + 4*pow(xq,2)*pow(y,3) - 
       pow(xk,2)*pow(xq,2)*pow(y,3) - 2*pow(y,4) + 2*pow(rl,2)*pow(xk,2)*(-1 + pow(xq,2))*(-1 + pow(xk,2) + y) + 
       2*pow(A,2)*(-1 + rl)*pow(-1 + pow(xk,2) + y,2) + A*(pow(xk,6)*(1 - 2*rl + pow(xq,2)) + 
          pow(xk,4)*(B*(2 - 4*rl + 2*pow(xq,2)) + (-3 + 4*rl - pow(xq,2))*(1 + pow(xq,2) - 2*y)) + 4*(-1 + rl)*(pow(xq,2) - y)*(-1 + y)*(-1 + B + y) + 
          pow(xk,2)*(-2*B*(2 + pow(xq,4) - 2*rl*(1 + pow(xq,2) - 2*y) - pow(xq,2)*(-1 + y) - 3*y) + 
             (-1 + y)*(-2 - pow(xq,4) + 2*rl*(1 + 4*pow(xq,2) - 5*y) + pow(xq,2)*(-7 + y) + 9*y))) + 
       rl*(2*pow(pow(xq,2) - y,2)*pow(-1 + B + y,2) + pow(xk,6)*(-1 + 2*B + pow(xq,2) + 2*y) + 
          pow(xk,4)*(3 + 2*pow(B,2) - 5*y + 6*pow(y,2) - pow(xq,2)*(1 + 3*y) + B*(-2 - 6*pow(xq,2) + 8*y)) + 
          2*pow(xk,2)*(-1 - pow(xq,4) + rll*pow(-1 + pow(xq,2),2) - 2*pow(B,2)*(pow(xq,2) - y) + 2*y + 4*pow(xq,2)*y + pow(xq,4)*y - 4*pow(y,2) - 
             4*pow(xq,2)*pow(y,2) + 3*pow(y,3) + B*(2*pow(xq,4) + pow(xq,2)*(3 - 7*y) + y*(-3 + 5*y))))))/
    (pow(xk,4)*pow(-1 + xq,2)*pow(1 + xq,2)*pow(pow(xk,2) - pow(xq,2) + y,2));


  
  double inth1 = (16*fk*pow(MK,5)*rl*(pow(xk,2) + 2*B*pow(xk,2) - 2*pow(B,2)*pow(xk,2) - rl*pow(xk,2) - pow(xk,4) - 2*B*pow(xk,4) + rl*pow(xk,4) + pow(xq,2) - 
       3*B*pow(xq,2) + 2*pow(B,2)*pow(xq,2) - 2*pow(xk,2)*pow(xq,2) + 3*B*pow(xk,2)*pow(xq,2) + pow(xk,4)*pow(xq,2) - y + 3*B*y - 2*pow(B,2)*y + 
       pow(xk,2)*y - 5*B*pow(xk,2)*y + rl*pow(xk,2)*y - pow(xk,4)*y - 2*pow(xq,2)*y + 3*B*pow(xq,2)*y + 2*pow(xk,2)*pow(xq,2)*y + 2*pow(y,2) - 
       3*B*pow(y,2) - 2*pow(xk,2)*pow(y,2) + pow(xq,2)*pow(y,2) - pow(y,3) + 2*rll*(-1 + pow(xq,2))*(-1 + pow(xk,2) + y) + 
					 A*(-1 + pow(xk,2) + y)*(-1 + 2*B + pow(xk,2) + y)))/((-1 + xq)*(1 + xq)*(-pow(xk,2) + pow(xq,2) - y));




  double inth2 = (4*fk*pow(MK,5)*rl*(pow(xq,2) - 3*B*pow(xq,2) + 2*pow(B,2)*pow(xq,2) - pow(xk,2)*pow(xq,2) + 7*B*pow(xk,2)*pow(xq,2) - 
       4*pow(B,2)*pow(xk,2)*pow(xq,2) - 4*B*pow(xk,4)*pow(xq,2) + 2*pow(xq,4) - 4*B*pow(xq,4) + 2*pow(B,2)*pow(xq,4) - 3*pow(xk,2)*pow(xq,4) + 
       5*B*pow(xk,2)*pow(xq,4) + pow(xq,6) - B*pow(xq,6) + 2*pow(rl,2)*pow(xk,2)*(-1 + pow(xq,2)) + 
       A*B*(pow(xk,2)*(2 + 6*pow(xq,2)) + rl*(4 - 8*pow(xk,2) + 4*pow(xq,2) - 8*y) - 2*(1 + pow(xq,2))*(1 + pow(xq,2) - 2*y)) + 
       2*rl*(-2*pow(B,2) + B*(3 + pow(xq,2) - 4*y) + (1 + pow(xq,2) - 2*y)*(-1 + y))*(pow(xq,2) - y) - y + 3*B*y - 2*pow(B,2)*y + 2*pow(xk,2)*y - 
       3*B*pow(xk,2)*y - pow(xk,4)*y - 5*pow(xq,2)*y + 8*B*pow(xq,2)*y - 2*pow(B,2)*pow(xq,2)*y + 6*pow(xk,2)*pow(xq,2)*y - 
       9*B*pow(xk,2)*pow(xq,2)*y - 3*pow(xk,4)*pow(xq,2)*y - 5*pow(xq,4)*y + 5*B*pow(xq,4)*y + 4*pow(xk,2)*pow(xq,4)*y - pow(xq,6)*y + 
       3*pow(y,2) - 4*B*pow(y,2) - 3*pow(xk,2)*pow(y,2) + 6*pow(xq,2)*pow(y,2) - 4*B*pow(xq,2)*pow(y,2) - 5*pow(xk,2)*pow(xq,2)*pow(y,2) + 
       3*pow(xq,4)*pow(y,2) - 2*pow(y,3) - 2*pow(xq,2)*pow(y,3) + 2*pow(A,2)*(-1 + 2*rl - pow(xq,2))*(-1 + pow(xk,2) + y) + 
       A*(pow(xk,2)*(1 + 3*pow(xq,2)) + rl*(2 - 4*pow(xk,2) + 6*pow(xq,2) - 8*y) - (1 + pow(xq,2))*(1 + 3*pow(xq,2) - 4*y))*(-1 + pow(xk,2) + y) + 
       rl*pow(xk,4)*(-1 + 4*B + pow(xq,2) + 4*y) + rl*pow(xk,2)*
        (1 + 4*pow(B,2) - pow(xq,4) + pow(xq,2)*(4 - 6*y) - 4*B*(1 + 2*pow(xq,2) - 3*y) - 6*y + 8*pow(y,2))))/
    ((-1 + xq)*(1 + xq)*(-pow(xk,2) + pow(xq,2) - y));

  double intfA = (8*fk*pow(MK,5)*rl*(-pow(xk,2) + 2*B*pow(xk,2) - 2*pow(B,2)*pow(xk,2) + 2*rl*pow(xk,2) + 2*pow(xk,4) - 2*B*pow(xk,4) - 2*rl*pow(xk,4) - 
       pow(xk,6) + rl*pow(xk,6) + 2*B*pow(xk,2)*pow(xq,2) - 2*rl*pow(xk,2)*pow(xq,2) - 2*pow(xk,4)*pow(xq,2) - B*pow(xk,4)*pow(xq,2) + 
       rl*pow(xk,4)*pow(xq,2) + pow(xk,6)*pow(xq,2) + 2*pow(xq,4) - 4*B*pow(xq,4) + 2*pow(B,2)*pow(xq,4) - 2*pow(xk,2)*pow(xq,4) + 
       3*B*pow(xk,2)*pow(xq,4) + pow(xk,4)*pow(xq,4) + 3*pow(xk,2)*y - 6*B*pow(xk,2)*y + 2*pow(B,2)*pow(xk,2)*y - 2*rl*pow(xk,2)*y - 
       4*pow(xk,4)*y + 3*B*pow(xk,4)*y + rl*pow(xk,4)*y + pow(xk,6)*y - 4*pow(xq,2)*y + 8*B*pow(xq,2)*y - 4*pow(B,2)*pow(xq,2)*y + 
       6*pow(xk,2)*pow(xq,2)*y - 8*B*pow(xk,2)*pow(xq,2)*y + 2*rl*pow(xk,2)*pow(xq,2)*y - pow(xk,4)*pow(xq,2)*y - 4*pow(xq,4)*y + 
       4*B*pow(xq,4)*y + 2*pow(xk,2)*pow(xq,4)*y + 2*pow(y,2) - 4*B*pow(y,2) + 2*pow(B,2)*pow(y,2) - 7*pow(xk,2)*pow(y,2) + 
       7*B*pow(xk,2)*pow(y,2) + 4*pow(xk,4)*pow(y,2) + 8*pow(xq,2)*pow(y,2) - 8*B*pow(xq,2)*pow(y,2) - 6*pow(xk,2)*pow(xq,2)*pow(y,2) + 
       2*pow(xq,4)*pow(y,2) - 4*pow(y,3) + 4*B*pow(y,3) + 5*pow(xk,2)*pow(y,3) - 4*pow(xq,2)*pow(y,3) + 2*pow(y,4) + 
       2*rll*(-1 + pow(xq,2))*(-1 + pow(xk,2) + pow(xq,2))*(-1 + pow(xk,2) + y) + 2*pow(A,2)*pow(-1 + pow(xk,2) + y,2) - 
       A*(pow(xk,6) + pow(xk,2)*(-2*B*(2 + pow(xq,2) - 3*y) - (2 + 7*pow(xq,2) - 9*y)*(-1 + y)) - 4*(pow(xq,2) - y)*(-1 + y)*(-1 + B + y) + 
          pow(xk,4)*(-3 + 2*B - 3*pow(xq,2) + 6*y))))/(pow(xk,2)*(-1 + xq)*(1 + xq)*(pow(xk,2) - pow(xq,2) + y));

  double intfV = (8*fk*pow(MK,5)*rl*(2*rll*(pow(xk,4) + (-1 + pow(xq,2))*(-1 + y) + pow(xk,2)*(-2 + 2*rl - pow(xq,2) + y)) + 
       pow(xk,2)*(1 + 2*pow(B,2) - 2*pow(xk,2) + rl*pow(xk,2) + pow(xk,4) - pow(xk,2)*pow(xq,2) - 2*y + 2*pow(xk,2)*y + pow(y,2) - 
		    A*(-1 + 2*B + pow(xk,2) + y) + B*(-2 + 2*pow(xk,2) - pow(xq,2) + 3*y))))/(pow(xk,2)*(pow(xk,2) - pow(xq,2) + y));

  double sdh1 = 8*pow(MK,6)*pow(xk,4)*(-2*rl*rll - pow(xq,2) + B*pow(xq,2) + 2*rll*pow(xq,2) + pow(xk,2)*pow(xq,2) + y - B*y - pow(xk,2)*y + pow(xq,2)*y - 
					   pow(y,2) + A*(-1 + 2*B + pow(xk,2) + y));


  double sdh2 = 4*pow(MK,6)*rl*pow(xk,4)*(rl - pow(xq,2))*(pow(A,2) + pow(B,2) + pow(xq,2) + A*(1 - 2*B - pow(xk,2) + pow(xq,2) - 2*y) - y + pow(xk,2)*y - 
     pow(xq,2)*y + pow(y,2) + B*(-1 + pow(xk,2) - pow(xq,2) + 2*y)) ;


  double sdfA = 2*pow(MK,6)*(2*rll*pow(xq,2) + pow(xk,2)*pow(xq,2) - 2*B*pow(xk,2)*pow(xq,2) + 2*pow(B,2)*pow(xk,2)*pow(xq,2) - 4*rll*pow(xk,2)*pow(xq,2) - 
     2*pow(xk,4)*pow(xq,2) + 2*B*pow(xk,4)*pow(xq,2) + 2*rll*pow(xk,4)*pow(xq,2) + pow(xk,6)*pow(xq,2) + 2*pow(xq,4) - 4*B*pow(xq,4) + 
     2*pow(B,2)*pow(xq,4) - 4*rll*pow(xq,4) - 2*pow(xk,2)*pow(xq,4) + 2*B*pow(xk,2)*pow(xq,4) + 4*rll*pow(xk,2)*pow(xq,4) + 
     2*pow(xk,4)*pow(xq,4) + 2*rll*pow(xq,6) + pow(xk,2)*pow(xq,6) - 2*rl*rll*pow(-1 + pow(xk,2) + pow(xq,2),2) - 4*pow(xq,2)*y + 
     8*B*pow(xq,2)*y - 4*pow(B,2)*pow(xq,2)*y + 4*pow(xk,2)*pow(xq,2)*y - 4*B*pow(xk,2)*pow(xq,2)*y - 4*pow(xq,4)*y + 4*B*pow(xq,4)*y + 
     2*pow(y,2) - 4*B*pow(y,2) + 2*pow(B,2)*pow(y,2) - 4*pow(xk,2)*pow(y,2) + 4*B*pow(xk,2)*pow(y,2) + 2*pow(xk,4)*pow(y,2) + 
     8*pow(xq,2)*pow(y,2) - 8*B*pow(xq,2)*pow(y,2) - 4*pow(xk,2)*pow(xq,2)*pow(y,2) + 2*pow(xq,4)*pow(y,2) - 4*pow(y,3) + 4*B*pow(y,3) + 
     4*pow(xk,2)*pow(y,3) - 4*pow(xq,2)*pow(y,3) + 2*pow(y,4) + 
     rl*pow(xk,2)*(1 - 2*pow(B,2) - 2*pow(xk,2) + pow(xk,4) - 2*pow(xq,2) - pow(xq,4) + 4*pow(xq,2)*y - 2*pow(y,2) - 
        2*B*(-1 + pow(xk,2) - pow(xq,2) + 2*y)) + 2*pow(A,2)*(pow(xk,4) + pow(-1 + y,2) + pow(xk,2)*(-2 - rl + pow(xq,2) + 2*y)) + 
     2*A*(pow(xk,2)*pow(xq,4) - 2*y*(-1 + pow(xk,2) + y)*(-1 + B + pow(xk,2) + y) + rl*pow(xk,2)*(-1 + 2*B + pow(xk,2) - pow(xq,2) + 2*y) + 
        pow(xq,2)*(pow(xk,4) + 2*B*(-1 + y) + 2*pow(-1 + y,2) + pow(xk,2)*(-3 + 2*y)))) ;


  double sdfV =  -2*pow(MK,6)*(-6*rl*rll - 8*pow(rl,2)*rll*pow(xk,2) - (2*rll + (1 + 2*(-1 + B)*B - 4*rll)*pow(xk,2) + 2*(-1 + B + rll)*pow(xk,4) + pow(xk,6))*
      pow(xq,2) + 2*(1 - 2*B + pow(B,2) - 2*rll + (-1 + B + 2*rll)*pow(xk,2) + pow(xk,4))*pow(xq,4) - (2*rll + pow(xk,2))*pow(xq,6) - 
     4*pow(xq,2)*(pow(B,2) + (-1 + 2*rll + pow(xk,2))*(-1 + pow(xk,2) - pow(xq,2)) + B*(-2 + 2*pow(xk,2) - pow(xq,2)))*y + 
     2*(pow(-1 + B + pow(xk,2),2) - 4*(-1 + B + rll + pow(xk,2))*pow(xq,2) + pow(xq,4))*pow(y,2) + 4*(-1 + B + pow(xk,2) - pow(xq,2))*pow(y,3) + 
     2*pow(y,4) + 2*pow(A,2)*(pow(xk,4) + pow(-1 + y,2) + pow(xk,2)*(-2 + rl - pow(xq,2) + 2*y)) - 
     rl*(pow(xk,6) + pow(xk,4)*(-2 + 6*B + 6*rll + 4*y) - 2*rll*(pow(xq,4) + pow(xq,2)*(2 - 4*y) + 4*y) + 
        pow(xk,2)*(1 + 6*pow(B,2) - pow(xq,2)*(-2 + pow(xq,2)) - 2*B*(3 + pow(xq,2) - 4*y) - 4*rll*(3 + pow(xq,2) - 2*y) + 2*(-2 + y)*y)) + 
     2*A*(-(pow(xk,2)*pow(xq,4)) + rl*pow(xk,2)*(-1 + 2*B + pow(xk,2) + pow(xq,2)) - 2*y*(-1 + pow(xk,2) + y)*(-1 + B + pow(xk,2) + y) + 
        pow(xq,2)*(pow(xk,4) + 2*B*(-1 + y) + 2*pow(-1 + y,2) + pow(xk,2)*(-3 + 4*y)))) ;


  if(sdfV < 0) crash ("sdfv < 0");




  

  double sdh1h2 = 8*pow(MK,6)*rl*pow(xk,4)*(-2*pow(B,2) + rl*pow(xk,2) - pow(xq,2) + B*(2 - 2*pow(xk,2) + pow(xq,2) - 3*y) + y - pow(xk,2)*y + pow(xq,2)*y - 
     pow(y,2) + A*(-1 + 2*B + pow(xk,2) + y)) ;

  double sdh1fA = -4*pow(MK,6)*pow(xk,2)*(pow(xq,2) - 3*B*pow(xq,2) + 2*pow(B,2)*pow(xq,2) - 4*rll*pow(xq,2) - 3*pow(xk,2)*pow(xq,2) + 3*B*pow(xk,2)*pow(xq,2) + 
     4*rll*pow(xk,2)*pow(xq,2) + 2*pow(xk,4)*pow(xq,2) - pow(xq,4) + B*pow(xq,4) + 4*rll*pow(xq,4) + 2*pow(xk,2)*pow(xq,4) + 
     2*A*B*(-1 + pow(xk,2) + pow(xq,2)) - 4*rl*rll*(-1 + pow(xk,2) + pow(xq,2)) - y + 3*B*y - 2*pow(B,2)*y + 2*pow(xk,2)*y - 3*B*pow(xk,2)*y - 
     pow(xk,4)*y + B*pow(xq,2)*y + pow(xq,4)*y + pow(y,2) - 2*B*pow(y,2) - pow(xk,2)*pow(y,2) - pow(xq,2)*pow(y,2) + 
     2*pow(A,2)*(-1 + pow(xk,2) + y) + A*(-1 + pow(xk,2) + 3*pow(xq,2) - 2*y)*(-1 + pow(xk,2) + y) + rl*pow(xk,2)*(-1 + pow(xk,2) - pow(xq,2) + 2*y)) ;

  double sdh1fV = -4*pow(MK,6)*pow(xk,2)*(pow(xq,2) - 3*B*pow(xq,2) + 2*pow(B,2)*pow(xq,2) - 4*rll*pow(xq,2) - 3*pow(xk,2)*pow(xq,2) + 3*B*pow(xk,2)*pow(xq,2) + 
     4*rll*pow(xk,2)*pow(xq,2) + 2*pow(xk,4)*pow(xq,2) + pow(xq,4) - B*pow(xq,4) - 4*rll*pow(xq,4) - 2*pow(xk,2)*pow(xq,4) + 
     rl*(4*rll + pow(xk,2))*(-1 + pow(xk,2) + pow(xq,2)) - y + 3*B*y - 2*pow(B,2)*y + 2*pow(xk,2)*y - 3*B*pow(xk,2)*y - pow(xk,4)*y - 
     4*pow(xq,2)*y + 5*B*pow(xq,2)*y + 8*rll*pow(xq,2)*y + 6*pow(xk,2)*pow(xq,2)*y - pow(xq,4)*y + 3*pow(y,2) - 4*B*pow(y,2) - 
     3*pow(xk,2)*pow(y,2) + 3*pow(xq,2)*pow(y,2) - 2*pow(y,3) - 2*pow(A,2)*(-1 + pow(xk,2) + y) + 2*A*B*(-1 + pow(xk,2) - pow(xq,2) + 2*y) + 
     A*(-1 + pow(xk,2) + y)*(-1 + pow(xk,2) - 3*pow(xq,2) + 4*y));

  double sdh2fA = -4*pow(MK,6)*rl*pow(xk,2)*(-pow(xq,2) + 3*B*pow(xq,2) - 2*pow(B,2)*pow(xq,2) + pow(xk,2)*pow(xq,2) - 3*B*pow(xk,2)*pow(xq,2) - pow(xq,4) + 
     B*pow(xq,4) + rl*pow(xk,2)*(-1 + pow(xk,2) + pow(xq,2)) + y - 3*B*y + 2*pow(B,2)*y - 2*pow(xk,2)*y + 3*B*pow(xk,2)*y + pow(xk,4)*y + 
     4*pow(xq,2)*y - 5*B*pow(xq,2)*y - 2*pow(xk,2)*pow(xq,2)*y + pow(xq,4)*y - 3*pow(y,2) + 4*B*pow(y,2) + 3*pow(xk,2)*pow(y,2) - 
     3*pow(xq,2)*pow(y,2) + 2*pow(y,3) + 2*pow(A,2)*(-1 + pow(xk,2) + y) - 
     A*(2*B*(-1 + pow(xk,2) - pow(xq,2) + 2*y) + (-1 + pow(xk,2) + y)*(-1 + pow(xk,2) - 3*pow(xq,2) + 4*y))) ;


  double sdfAfV = 4*pow(MK,6)*(pow(xq,2)*(2*rll*(-1 + pow(xk,2) + pow(xq,2))*(-1 + pow(xk,2) - pow(xq,2) + 2*y) + 
        pow(xk,2)*(-2*pow(A,2) + 2*pow(B,2) - 2*A*(pow(xq,2) - y) + 2*B*(-1 + pow(xk,2) + y) + 
           (-1 + pow(xk,2) + pow(xq,2))*(-1 + pow(xk,2) - pow(xq,2) + 2*y))) + 
     rl*(2*rll*pow(-1 + pow(xk,2) + pow(xq,2),2) + pow(xk,2)*
         (1 + 2*pow(A,2) + 2*pow(B,2) - 2*pow(xk,2) + pow(xk,4) + pow(xq,4) - 2*y + 2*pow(xk,2)*y - 2*pow(xq,2)*y + 2*pow(y,2) + 
           2*B*(-1 + pow(xk,2) - pow(xq,2) + 2*y) - 2*A*(-1 + 2*B + pow(xk,2) - pow(xq,2) + 2*y)))) ;



 //pt
 pt *= pow(fk,2);


 //int
 inth1 *= real(h1);
 inth2 *= real(h2);
 intfA *= real(fA);
 intfV *= real(fV);



 //sd

 sdh1 *= norm(h1);
 sdh2 *= norm(h2);
 sdfA *= norm(fA);
 sdfV *= norm(fV);
 sdh1h2 *= real(h1*conj(h2));
 sdh1fA *= real(h1*conj(fA));
 sdh1fV *= real(h1*conj(fV));
 sdh2fA *= real(h2*conj(fA));
 sdfAfV *= real(fA*conj(fV));



 //Kahan sum the terms

 vector<double> Int_Kernel_To_Sum({inth1,inth2,intfA,intfV});
 vector<double> SD_Kernel_To_Sum({sdh1,sdh2,sdfA,sdfV,sdh1h2,sdh1fA,sdh1fV, sdh2fA, sdfAfV});


 double SD_Kernel_To_Sum_naive= sdh1+sdh2+sdfA+sdfV+sdh1h2+sdh1fA+sdh1fV+ sdh2fA+sdfAfV;
 double Int_Kernel_To_Sum_naive = inth1 + inth2 + intfA + intfV;


 double Int_Kernel= Kahan_sum(Int_Kernel_To_Sum);

 double SD_Kernel = Kahan_sum(SD_Kernel_To_Sum);

 long double result;
 if(MODE=="PT") result= pt;
 else if(MODE=="INTERFERENCE") result= Int_Kernel;
 else if(MODE=="QUADRATIC") result= SD_Kernel;
 else if(MODE=="SD") result = SD_Kernel+Int_Kernel;
 else {
    vector<double> summ_vector({pt, Int_Kernel, SD_Kernel});
    result=Kahan_sum(summ_vector);
 }



 
   

 if(result <0 && ( MODE != "INTERFERENCE" && MODE != "SD") ) {
   cout<<"Kernel unequal leptons is negative!"<<endl;
    cout<<"Printing info...."<<endl;
    cout<<"tot: "<<result<<endl;
    cout<<"point-like: "<<pt<<" "<<endl;
    cout<<"interference: "<<Int_Kernel<<endl;
    cout<<"Quadratic: "<<SD_Kernel<<"   Naive: "<<SD_Kernel_To_Sum_naive<<endl;
    cout<<"xk: "<<xk<<" xq: "<<xq<<endl;
    cout<<"#############"<<endl<<flush;
    cout<<sdh1<<" "<<sdh2<<" "<<sdfA<<" "<<sdfV<<" "<<sdh1h2<<" "<<sdh1fA<<" "<<sdh1fV<<" "<<sdh2fA<<" "<<sdfAfV<<" "<<endl;

    cout<<"###########"<<endl<<flush;
   
    crash("");
 }

 return result;

}



long double Compute_square_amplitude_extended( complex<double> h1, complex<double> h2, complex<double> fA, complex<double> fV, complex<double> he1, complex<double> he2, complex<double> feA, complex<double> feV, double MK, double fk,double rl,  double xk, double xq, double A, double B, double y, string MODE) {


  double xk_prime= sqrt(2*rl + A);
  double xq_prime= sqrt(rl +B);
  h1 /= pow(xk,2);
  h2 /= pow(xk,2)*(1.0-pow(xq,2));
  fA /= pow(xk,2);
  fV /= LeviCivitaSign*pow(xk,2);
  he1 /= pow(xk_prime,2);
  he2 /= pow(xk_prime,2)*(1.0-pow(xq_prime,2));
  feA /= pow(xk_prime,2);
  feV /= LeviCivitaSign*pow(xk_prime,2);

  
  double pt =(8*pow(MK,4)*rl*(2*pow(A,4)*(-1 + rl)*pow(-1 + B + rl,2)*pow(-1 + pow(xk,2) + y,2) - pow(A,3)*(-1 + B + rl)*(-8*pow(rl,3)*pow(-1 + pow(xk,2) + y,2) + 2*pow(rl,2)*(-1 + pow(xk,2) + y)*(pow(xk,4) - 2*(-4 + pow(xq,2) - y)*(-1 + y) - 2*B*(-2 + pow(xk,2) + pow(xq,2) + y) +
                                      pow(xk,2)*(8 - 3*pow(xq,2) + 3*y)) - 2*pow(B,2)*
                                    (pow(xk,4)*(1 + pow(xq,2)) - 2*(pow(xq,2) - y)*(-1 + y) - pow(xk,2)*(2 + pow(xq,2) + pow(xq,4) - 3*y - pow(xq,2)*y)) -
                                   2*B*(pow(xk,6)*pow(xq,2) - pow(xk,4)*(1 + pow(xq,4) - 2*pow(xq,2)*(-2 + y) - 2*y) - 2*(pow(xq,2) - y)*(2 - 3*y + pow(y,2)) +
                                      pow(xk,2)*(2 - pow(xq,4)*(-2 + y) - 7*y + 4*pow(y,2) + pow(xq,2)*(5 - 6*y + pow(y,2)))) -
                                   (-1 + pow(xk,2) + y)*(pow(xk,6)*(-1 + pow(xq,2)) + 4*(pow(xq,2) - y)*(-1 + y) - pow(xk,4)*(-1 + pow(xq,4) - 2*pow(xq,2)*(-1 + y) + 2*y) -
                                      pow(xk,2)*(pow(xq,4)*(-3 + y) + y*(3 + y) - pow(xq,2)*(3 - 2*y + pow(y,2)))) +
                                   2*rl*(2*pow(B,2)*(-1 + pow(xk,2) + y)*(pow(xk,2) - pow(xq,2) + y) -
                                      (-1 + pow(xk,2) + y)*(pow(xk,4)*(1 + pow(xq,2)) + pow(xk,2)*(4 - pow(xq,4) + pow(xq,2)*(-5 + y) + 5*y) +
                                         4*(-1 + pow(xq,2) - pow(xq,2)*y + pow(y,2))) +
                                      B*(pow(xk,6) - 4*pow(xk,4)*(pow(xq,2) - y) - 2*(-1 + y)*(2 + pow(xq,2)*(-3 + y) + y - pow(y,2)) +
                                         pow(xk,2)*(-4 + pow(xq,4) + pow(xq,2)*(8 - 6*y) - 4*y + 5*pow(y,2))))) +
                                pow(A,2)*(-1 + B + rl)*(pow(xk,8)*(-1 + pow(xq,2))*(-1 - B + pow(xq,2) - 2*y) - 2*(-1 + B)*pow(pow(xq,2) - y,2)*pow(-1 + B + y,2) +
                                   8*pow(rl,4)*pow(-1 + pow(xk,2) + y,2) - 2*pow(rl,3)*
                                    (4*pow(xk,6) - 8*(-1 + pow(xq,2) - y)*pow(-1 + y,2) + 4*B*(-1 + pow(xk,2) + y)*(1 + pow(xk,2) - 2*pow(xq,2) + y) +
                                      pow(xk,4)*(7 - 15*pow(xq,2) + 16*y) - pow(xk,2)*(20 - 25*pow(xq,2) + pow(xq,4) + y + 23*pow(xq,2)*y - 20*pow(y,2))) +
                                   pow(xk,6)*(-2 - pow(xq,6) + pow(B,2)*(1 - 3*pow(xq,2)) - 3*y + 6*pow(y,2) + pow(xq,4)*(-3 + 5*y) +
                                      B*(-2 + pow(xq,2) + 3*pow(xq,4) + 5*y - 7*pow(xq,2)*y) - 6*pow(xq,2)*(-1 + pow(y,2))) +
                                   pow(xk,4)*(1 - 2*pow(B,3)*pow(xq,2) + pow(B,2)*(-3 + 5*pow(xq,4) - 8*pow(xq,2)*(-1 + y)) + pow(xq,6)*(5 - 2*y) + 4*y - 7*pow(y,2) +
                                      6*pow(y,3) + pow(xq,4)*(1 - 12*y + 7*pow(y,2)) + pow(xq,2)*(-5 + 2*y + 6*pow(y,2) - 6*pow(y,3)) +
                                      B*(1 - 2*pow(xq,6) - 6*y + 5*pow(y,2) + 2*pow(xq,4)*(-6 + 5*y) + pow(xq,2)*(3 + 12*y - 11*pow(y,2)))) +
                                   pow(rl,2)*(pow(xk,6)*(3 + 13*pow(xq,2) + 2*y) + pow(xk,4)*(19 - 7*pow(xq,4) + 38*y + 6*pow(y,2) + 6*pow(xq,2)*(-11 + 3*y)) -
                                      2*pow(B,2)*(7*pow(xk,4) - pow(xq,4) - 2*pow(xk,2)*(4 + 3*pow(xq,2) - 7*y) + pow(xq,2)*(8 - 6*y) + y*(-8 + 7*y)) +
                                      2*(-1 + y)*(-4 - 3*pow(xq,4) + pow(xq,6) - 13*y + 15*pow(y,2) + pow(y,3) + pow(xq,2)*(17 - 12*y - 2*pow(y,2))) +
                                      pow(xk,2)*(-30 + pow(xq,4)*(5 - 9*y) - 43*y + 63*pow(y,2) + 6*pow(y,3) + pow(xq,2)*(85 - 78*y + pow(y,2))) -
                                      2*B*(3*pow(xk,6) - 3*pow(xk,4)*(3 + 5*pow(xq,2) - 4*y) - 2*(-1 + y)*(2 + pow(xq,4) + 2*pow(xq,2)*(-6 + y) + 10*y - 3*pow(y,2)) +
                                         pow(xk,2)*(4 + pow(xq,2)*(35 - 19*y) - 35*y + 15*pow(y,2)))) -
                                   pow(xk,2)*(pow(xq,6)*(2 + 2*pow(B,2) + 2*B*(-2 + y) - 3*y + pow(y,2)) +
                                      y*(3 + 2*pow(B,3) - 2*y + pow(y,2) - 2*pow(y,3) + pow(B,2)*(-3 + 5*y) + B*(-1 - 4*y + pow(y,2))) -
                                      pow(xq,4)*(-5 + 2*pow(B,3) + 11*y - 9*pow(y,2) + 3*pow(y,3) + pow(B,2)*(-9 + 7*y) + B*(11 - 18*y + 7*pow(y,2))) +
                                      pow(xq,2)*(-3 + 2*pow(B,3)*(-1 + y) - 3*y + 8*pow(y,2) - 4*pow(y,3) + 2*pow(y,4) + pow(B,2)*(3 - 14*y + 5*pow(y,2)) +
                                         B*(1 + 15*y - 15*pow(y,2) + 5*pow(y,3)))) + rl*
                                    (5*pow(xk,8)*(-1 + pow(xq,2)) + pow(xk,6)*(14 + 2*pow(B,2) - 3*pow(xq,4) - 15*y + 2*B*(-2 + 4*pow(xq,2) + y) + pow(xq,2)*(-19 + 11*y)) +
                                      2*(pow(xq,2) - y)*(pow(B,3)*(pow(xq,2) - y) + B*(-1 + y)*(17 + pow(xq,4) + pow(xq,2)*(-7 + y) - 3*y - pow(y,2)) +
                                         2*pow(B,2)*(4 + pow(xq,2)*(-2 + y) - 2*y - pow(y,2)) + (-1 + y)*(-9 - pow(xq,4) - 2*pow(xq,2)*(-2 + y) + 6*y + 2*pow(y,2))) +
                                      pow(xk,4)*(-17 + 2*pow(B,3) - 2*pow(xq,6) + 6*B*pow(xq,2)*(-5 + y) - 4*pow(xq,4)*(-3 + y) + 6*y - 19*pow(y,2) + 2*B*y*(1 + 3*y) +
                                         pow(B,2)*(6 - 2*pow(xq,2) + 8*y) + pow(xq,2)*(35 - 18*y + 7*pow(y,2))) +
                                      pow(xk,2)*(6 - 4*pow(B,3)*(pow(xq,2) - y) + 25*y - 16*pow(y,2) - 13*pow(y,3) - pow(xq,4)*(3 - 6*y + pow(y,2)) +
                                         2*pow(B,2)*(-8 + pow(xq,2) + 7*y - 5*pow(xq,2)*y + 5*pow(y,2)) + pow(xq,2)*(-35 + 21*y + 9*pow(y,2) + pow(y,3)) -
                                         2*B*(-5 + pow(xq,6) + 23*y - 5*pow(y,2) - 3*pow(y,3) + pow(xq,2)*(-20 + 6*y + 3*pow(y,2)))))) +
                                2*(4*pow(rl,6)*pow(xk,2)*(-1 + pow(xq,2))*(-2 + pow(xk,2) + pow(xq,2) + y) -
                                   pow(xk,4)*pow(-1 + pow(xq,2),2)*pow(pow(xq,2)*(-1 + B + y) - y*(-1 + B + pow(xk,2) + y),2) +
                                   pow(rl,5)*(pow(xk,6)*(1 + 4*B - pow(xq,2) + 4*y) +
                                      pow(xk,4)*(20 + 4*pow(B,2) + 17*pow(xq,4) - 4*B*(1 + 3*pow(xq,2) - 4*y) + y + 12*pow(y,2) - pow(xq,2)*(33 + 17*y)) +
                                      4*(pow(xq,2) - y)*(-1 + pow(xq,4)*(-1 + y) - pow(B,2)*y - 2*B*(-1 + y)*y + 2*pow(y,2) - pow(y,3) +
                                         pow(xq,2)*(3 + pow(B,2) + 2*B*(-1 + y) - 4*y + pow(y,2))) -
                                      4*pow(xk,2)*(5 + pow(xq,4)*(5 - 2*y) + 2*pow(B,2)*(pow(xq,2) - y) - y + 2*pow(y,2) - 3*pow(y,3) + pow(xq,2)*(-8 - 4*y + 6*pow(y,2)) +
                                         B*(-4 - 6*pow(xq,4) + 3*y - 5*pow(y,2) + pow(xq,2)*(5 + 7*y)))) +
                                   pow(rl,4)*(-2*pow(xk,8)*(-1 + pow(xq,2)) + pow(xk,6)*(-3 + 8*pow(B,2) + 13*pow(xq,4) - 2*pow(xq,2)*(5 + 6*y) + B*(1 - 13*pow(xq,2) + 8*y)) +
                                      4*(pow(xq,2) - y)*(2*pow(B,3)*(pow(xq,2) - y) + pow(B,2)*(pow(xq,2) - y)*(-7 + 4*y) +
                                         2*B*(-1 + y)*(1 + pow(xq,4) + pow(xq,2)*(-6 + y) + 4*y - pow(y,2)) + (-1 + y)*(-2 - 2*pow(xq,4) + pow(xq,2)*(7 - 3*y) - 3*y + 3*pow(y,2)))\
                                       + pow(xk,4)*(-28 + 8*pow(B,3) - 11*pow(xq,6) - 4*pow(B,2)*(1 + 10*pow(xq,2) - 8*y) - 18*y - 19*pow(y,2) + pow(xq,4)*(-39 + 20*y) +
                                         pow(xq,2)*(66 + 46*y - 17*pow(y,2)) + 4*B*(4 + 13*pow(xq,4) - 6*y + 6*pow(y,2) - pow(xq,2)*(3 + 14*y))) +
                                      pow(xk,2)*(16 - 16*pow(B,3)*(pow(xq,2) - y) + pow(xq,6)*(-1 + y) + 16*y + 9*pow(y,2) - 29*pow(y,3) +
                                         8*pow(B,2)*(1 + 5*pow(xq,4) + pow(xq,2)*(2 - 9*y) - 4*y + 5*pow(y,2)) + pow(xq,4)*(49 - 43*y + 6*pow(y,2)) -
                                         pow(xq,2)*(40 + 58*y - 81*pow(y,2) + 7*pow(y,3)) +
                                         B*(-24 - 7*pow(xq,6) + 8*y - 65*pow(y,2) + 24*pow(y,3) + pow(xq,4)*(-77 + 34*y) + pow(xq,2)*(32 + 134*y - 59*pow(y,2))))) -
                                   pow(rl,2)*(4*pow(-1 + B,2)*pow(pow(xq,2) - y,2)*pow(-1 + B + y,2) +
                                      pow(xk,8)*(-1 + pow(xq,2))*(4 + 5*pow(B,2) + pow(xq,2)*(-2 + y) - 4*y - pow(y,2) + B*(-8 + pow(xq,2) + 2*y)) +
                                      pow(xk,6)*(4 + pow(B,3)*(-5 + 9*pow(xq,2)) - 14*y + 4*pow(y,2) + 3*pow(y,3) + 2*pow(xq,6)*(1 + y) + pow(xq,4)*(-14 - 7*y + 2*pow(y,2)) +
                                         pow(xq,2)*(8 + 23*y - 6*pow(y,2) - 3*pow(y,3)) + pow(B,2)*(19 - 10*pow(xq,4) - 11*y + pow(xq,2)*(-17 + 15*y)) +
                                         B*(pow(xq,4)*(22 + 3*y) - 2*(8 - 12*y + pow(y,2)) + pow(xq,2)*(-2 - 35*y + 2*pow(y,2)))) +
                                      pow(xk,4)*(2 + 4*pow(B,4)*pow(xq,2) - 4*pow(B,3)*(-2 + 3*pow(xq,4) + pow(xq,2)*(5 - 4*y)) - 3*pow(xq,8)*(-1 + y) + 16*y - 4*pow(y,2) -
                                         2*pow(y,3) + 3*pow(y,4) + pow(xq,6)*(1 + 5*y + 3*pow(y,2)) + pow(xq,4)*(14 - 12*y - 12*pow(y,2) + 3*pow(y,3)) -
                                         pow(xq,2)*(16 + 22*y - 25*pow(y,2) + pow(y,3) + 3*pow(y,4)) +
                                         2*pow(B,2)*(-6 + 3*pow(xq,6) + pow(xq,4)*(22 - 7*y) + 9*y + pow(xq,2)*(1 - 26*y + 6*pow(y,2))) +
                                         B*(-9*pow(xq,6) - 3*pow(xq,8) + 7*pow(xq,4)*(-6 + 4*y + pow(y,2)) + pow(xq,2)*(28 + 52*y - 37*pow(y,2) - 2*pow(y,3)) +
                                            2*(1 - 16*y + 3*pow(y,2) + pow(y,3)))) - pow(xk,2)*(pow(xq,2) - y)*
                                       (4*pow(B,4)*(1 + pow(xq,2)) + pow(B,3)*(-8 - 3*pow(xq,4) + 13*y + pow(xq,2)*(-25 + 7*y)) +
                                         pow(B,2)*(-4 + pow(xq,6) - 25*y + 10*pow(y,2) + pow(xq,4)*(6 + y) + pow(xq,2)*(53 - 36*y + 2*pow(y,2))) +
                                         (-1 + y)*(10 + pow(xq,6) + 6*y - pow(y,2) + pow(y,3) + pow(xq,4)*(1 - 4*y + pow(y,2)) - pow(xq,2)*(20 - 10*y + pow(y,3))) -
                                         B*(pow(xq,6)*y + pow(xq,4)*(2 + 5*y - 5*pow(y,2)) - 2*(9 + 4*y - 8*pow(y,2) + pow(y,3)) +
                                            pow(xq,2)*(52 - 58*y + 13*pow(y,2) + 2*pow(y,3))))) +
                                   rl*pow(xk,2)*(-1 + pow(xq,2))*(pow(xk,6)*(2 - y - 2*pow(y,2) + pow(B,2)*(2 - 2*pow(xq,2) + y) + B*(-2 - pow(xq,2)*(-2 + y) + pow(y,2)) +
                                         pow(xq,2)*(-2 + y + pow(y,2))) + (-1 + B)*pow(pow(xq,2) - y,2)*
                                       (-2 + pow(B,3) + 3*pow(B,2)*(-1 + y) + 4*y - 3*pow(y,2) + pow(y,3) + B*(4 - 6*y + 3*pow(y,2))) -
                                      pow(xk,4)*(2 - 3*y - 5*pow(y,2) + 5*pow(y,3) + pow(xq,4)*(-3 - 4*y + 2*pow(y,2)) + pow(xq,2)*(1 + 10*y - 2*pow(y,2) - 2*pow(y,3)) +
                                         pow(B,3)*(3*pow(xq,2) - 2*(1 + y)) + pow(B,2)*(4 - 3*pow(xq,4) - 5*pow(y,2) + pow(xq,2)*(-3 + 8*y)) +
                                         B*(-4 + y + 8*pow(y,2) - 3*pow(y,3) + 2*pow(xq,4)*(1 + y) + pow(xq,2)*(3 - 12*y + 2*pow(y,2)))) +
                                      pow(xk,2)*(pow(B,4)*(-pow(xq,2) + y) + pow(y,2)*(-6 + 9*y - 4*pow(y,2)) + pow(xq,6)*(2 - 5*y + pow(y,2)) -
                                         2*pow(xq,4)*(2 - 2*y - 3*pow(y,2) + pow(y,3)) + pow(xq,2)*y*(10 - 15*y + 3*pow(y,2) + pow(y,3)) +
                                         pow(B,3)*(3*pow(xq,4) + pow(xq,2)*(2 - 8*y) + y*(-2 + 5*y)) +
                                         pow(B,2)*(pow(xq,4)*(-6 + 4*y) + pow(xq,2)*(-3 + 18*y - 11*pow(y,2)) + y*(3 - 12*y + 7*pow(y,2))) +
                                         B*(pow(xq,6)*(-4 + 3*y) + pow(xq,4)*(9 - 2*y - 3*pow(y,2)) + pow(xq,2)*(2 - 24*y + 20*pow(y,2) - 3*pow(y,3)) +
                                            y*(-2 + 15*y - 14*pow(y,2) + 3*pow(y,3))))) +
                                   pow(rl,3)*(-(pow(xk,8)*(-1 + pow(xq,2))*(-4 + 7*B + 3*y)) +
                                      pow(xk,6)*(4 + 4*pow(B,3) - pow(xq,6) - 15*y + 7*pow(y,2) + pow(xq,4)*(-23 + 3*y) + pow(B,2)*(5 - 21*pow(xq,2) + 4*y) +
                                         B*(-18 + 26*pow(xq,4) + pow(xq,2)*(4 - 29*y) + 13*y) + pow(xq,2)*(20 + 24*y - 7*pow(y,2))) +
                                      4*(-1 + B)*(pow(xq,2) - y)*(pow(B,3)*(pow(xq,2) - y) + pow(B,2)*(pow(xq,2) - y)*(-5 + 2*y) +
                                         B*(-1 + y)*(1 + pow(xq,4) + pow(xq,2)*(-9 + y) + 7*y - pow(y,2)) - (-1 + y)*(1 + pow(xq,4) + 3*y - 3*pow(y,2) + pow(xq,2)*(-5 + 3*y))) +
                                      pow(xk,4)*(14 + 4*pow(B,4) + 2*pow(xq,8) - 4*pow(B,3)*(1 + 7*pow(xq,2) - 4*y) - 2*pow(xq,6)*(-8 + y) + 33*y - 2*pow(y,2) + 5*pow(y,3) +
                                         3*pow(xq,4)*(11 - 12*y + pow(y,2)) - pow(xq,2)*(53 + 43*y - 35*pow(y,2) + 5*pow(y,3)) +
                                         4*pow(B,2)*(-2 + 12*pow(xq,4) - 6*y + 3*pow(y,2) - 2*pow(xq,2)*(-5 + 7*y)) -
                                         2*B*(4 + 10*pow(xq,6) + pow(xq,4)*(43 - 21*y) + 13*y + 7*pow(y,2) + pow(xq,2)*(-25 - 48*y + 17*pow(y,2)))) +
                                      pow(xk,2)*(-4 - 8*pow(B,4)*(pow(xq,2) - y) - pow(xq,8)*(-1 + y) - 26*y + 3*pow(xq,6)*(-1 + y)*y + 5*pow(y,2) + 20*pow(y,3) + pow(y,4) -
                                         pow(xq,4)*(51 - 64*y + 16*pow(y,2) + pow(y,3)) + pow(xq,2)*(30 + 50*y - 88*pow(y,2) + 17*pow(y,3) - pow(y,4)) +
                                         4*pow(B,3)*(6*pow(xq,4) + pow(xq,2)*(7 - 11*y) + y*(-7 + 5*y)) -
                                         pow(B,2)*(4 + 11*pow(xq,6) + pow(xq,4)*(81 - 38*y) - 20*y + 69*pow(y,2) - 12*pow(y,3) + pow(xq,2)*(16 - 154*y + 43*pow(y,2))) +
                                         B*(8 + pow(xq,8) + pow(xq,6)*(11 - 2*y) + 24*y + 48*pow(y,2) - 36*pow(y,3) + pow(xq,4)*(104 - 94*y + 13*pow(y,2)) -
                                            pow(xq,2)*(32 + 160*y - 127*pow(y,2) + 12*pow(y,3)))))) +
                                A*(-8*pow(rl,5)*(pow(xk,6) - 2*(pow(xq,2) - y)*(-1 + y)*(-1 + B + y) + pow(xk,4)*(1 + 2*B - 5*pow(xq,2) + 4*y) -
                                      pow(xk,2)*(3 - 9*pow(xq,2) + pow(xq,4) + 2*B*(1 + pow(xq,2) - 2*y) + 3*y + 7*pow(xq,2)*y - 5*pow(y,2))) +
                                   pow(rl,4)*(pow(xk,6)*(7 - 8*B + 17*pow(xq,2) + 8*y) +
                                      pow(xk,4)*(50 - 24*pow(B,2) + 7*pow(xq,4) + 4*B*(7 + 13*pow(xq,2) - 8*y) + 61*y + 24*pow(y,2) + pow(xq,2)*(-145 + 3*y)) +
                                      8*(pow(xq,2) - y)*(2*B*(-1 + y)*(-5 + pow(xq,2) + y) + pow(B,2)*(-4 + pow(xq,2) + 3*y) +
                                         (-1 + y)*(7 + pow(xq,4) + pow(xq,2)*(-3 + y) - 5*y - pow(y,2))) +
                                      2*pow(xk,2)*(-32 + 8*pow(B,2)*(2 + pow(xq,2) - 3*y) - 35*y + 43*pow(y,2) + 12*pow(y,3) - pow(xq,4)*(11 + y) +
                                         pow(xq,2)*(95 - 56*y - 15*pow(y,2)) + 2*B*(9*pow(xq,4) + (31 - 10*y)*y + pow(xq,2)*(-43 + 13*y)))) +
                                   pow(rl,3)*(4*pow(xk,8)*(-1 + pow(xq,2)) + pow(xk,6)*
                                       (16 + 8*pow(B,2) + 8*pow(xq,4) - 22*y - 2*pow(xq,2)*(24 + y) + B*(-5 + 13*pow(xq,2) + 16*y)) +
                                      8*(pow(xq,2) - y)*(2*pow(B,3)*(-1 + pow(xq,2)) + 2*B*(-1 + y)*(8 + pow(xq,4) + pow(xq,2)*(-6 + y) - pow(y,2)) +
                                         (-1 + y)*(-8 - 2*pow(xq,4) + pow(xq,2)*(7 - 3*y) + 3*y + 3*pow(y,2)) + pow(B,2)*(10 - 5*y - 2*pow(y,2) + pow(xq,2)*(-7 + 4*y))) +
                                      pow(xk,4)*(-76 - 8*pow(xq,6) - 41*y - 57*pow(y,2) + pow(xq,4)*(-27 + 13*y) + pow(xq,2)*(183 + 28*y - 15*pow(y,2)) -
                                         8*pow(B,2)*(3*pow(xq,2) - 4*(1 + y)) + 2*B*(13 + 23*pow(xq,4) + 4*y + 24*pow(y,2) - 4*pow(xq,2)*(17 + 5*y))) -
                                      pow(xk,2)*(-56 + 16*pow(B,3)*(-1 + pow(xq,2)) - pow(xq,6)*(-1 + y) - 84*y + 53*pow(y,2) + 63*pow(y,3) +
                                         pow(xq,4)*(-53 + 37*y - 8*pow(y,2)) + 3*pow(xq,2)*(60 - 8*y - 39*pow(y,2) + 3*pow(y,3)) -
                                         8*pow(B,2)*(-5 + 7*pow(xq,4) + 9*y + 5*pow(y,2) - pow(xq,2)*(7 + 9*y)) +
                                         B*(32 + 11*pow(xq,6) + pow(xq,4)*(99 - 32*y) + 150*y + 3*pow(y,2) - 48*pow(y,3) + pow(xq,2)*(-246 - 54*y + 85*pow(y,2))))) +
                                   pow(xk,2)*(-1 + pow(xq,2))*(-(pow(xk,6)*(-1 + pow(xq,2) + pow(B,2)*(-1 + pow(xq,2) - y) + 3*y - 3*pow(xq,2)*y + B*(pow(xq,2) - y)*y +
                                           pow(y,2))) + (-1 + B)*pow(pow(xq,2) - y,2)*
                                       (-2 + pow(B,3) + 3*pow(B,2)*(-1 + y) + 4*y - 3*pow(y,2) + pow(y,3) + B*(4 - 6*y + 3*pow(y,2))) +
                                      pow(xk,4)*(-1 + pow(xq,4)*(3 - 4*y) + 7*y - 3*pow(y,2) - 3*pow(y,3) + pow(B,3)*(1 - 2*pow(xq,2) + 2*y) +
                                         pow(B,2)*(-1 + pow(xq,4) + pow(xq,2)*(2 - 6*y) - 2*y + 5*pow(y,2)) + pow(xq,2)*(-2 - 6*y + 8*pow(y,2)) +
                                         B*(1 - 3*y - 6*pow(y,2) + 3*pow(y,3) - 2*pow(xq,2)*(1 - 6*y + 2*pow(y,2)))) +
                                      pow(xk,2)*(pow(xq,6)*(-2 + y + B*y) + pow(xq,4)*(-3 + 2*pow(B,3) + 10*y - 5*pow(y,2) + pow(B,2)*(-7 + 4*y) + B*(10 - 12*y + pow(y,2))) +
                                         y*(-3 + pow(B,4) + 2*y + 3*pow(y,2) - 3*pow(y,3) + pow(B,3)*(-3 + 5*y) + pow(B,2)*(2 - 12*y + 7*pow(y,2)) +
                                            B*(3 + 7*y - 12*pow(y,2) + 3*pow(y,3))) - pow(xq,2)*
                                          (pow(B,4) + pow(B,3)*(-3 + 7*y) - pow(-1 + y,2)*(3 + 7*y) + pow(B,2)*(2 - 19*y + 11*pow(y,2)) + B*(3 + 17*y - 24*pow(y,2) + 5*pow(y,3)))
                                         )) + pow(rl,2)*(pow(xk,8)*(-1 + pow(xq,2))*(-4 - 3*B + 2*pow(xq,2) - 7*y) +
                                      8*(-1 + B)*(pow(xq,2) - y)*(pow(B,3)*(pow(xq,2) - y) + B*(-1 + y)*(5 + pow(xq,4) + pow(xq,2)*(-9 + y) + 5*y - pow(y,2)) +
                                         pow(B,2)*(2 + 3*y - 2*pow(y,2) + pow(xq,2)*(-5 + 2*y)) - (-1 + y)*(3 + pow(xq,4) + y - 3*pow(y,2) + pow(xq,2)*(-5 + 3*y))) +
                                      pow(xk,6)*(-15 + 8*pow(B,3) - pow(xq,6) - 9*y + 19*pow(y,2) + pow(xq,4)*(-26 + 7*y) + pow(B,2)*(-5 - 19*pow(xq,2) + 8*y) +
                                         pow(xq,2)*(50 + 26*y - 19*pow(y,2)) + B*(5 + 29*pow(xq,4) + y - pow(xq,2)*(26 + 33*y))) +
                                      pow(xk,4)*(34 + 8*pow(B,4) - pow(xq,8) - 4*pow(B,3)*(1 + 11*pow(xq,2) - 8*y) + 46*y + 3*pow(y,2) + 17*pow(y,3) + pow(xq,6)*(19 + 4*y) +
                                         pow(xq,4)*(46 - 63*y + 10*pow(y,2)) + 2*pow(B,2)*(-7 + 31*pow(xq,4) + pow(xq,2)*(24 - 38*y) - 26*y + 12*pow(y,2)) -
                                         pow(xq,2)*(106 + 51*y - 59*pow(y,2) + 17*pow(y,3)) -
                                         B*(28 + 19*pow(xq,6) + pow(xq,4)*(112 - 53*y) + 23*y + 42*pow(y,2) + pow(xq,2)*(-111 - 130*y + 54*pow(y,2)))) +
                                      pow(xk,2)*(-16 - 16*pow(B,4)*(pow(xq,2) - y) - pow(xq,8)*(-1 + y) - 54*y + 17*pow(y,2) + 40*pow(y,3) + 5*pow(y,4) +
                                         pow(xq,6)*(-4 + 3*y + pow(y,2)) + 4*pow(B,3)*(-4 + 9*pow(xq,4) + pow(xq,2)*(13 - 19*y) - 9*y + 10*pow(y,2)) +
                                         pow(xq,4)*(-67 + 90*y - 36*pow(y,2) + 5*pow(y,3)) + pow(xq,2)*(78 + 50*y - 134*pow(y,2) + 27*pow(y,3) - 5*pow(y,4)) +
                                         pow(B,2)*(16 - 19*pow(xq,6) - 14*y - 111*pow(y,2) + 24*pow(y,3) + 3*pow(xq,4)*(-41 + 20*y) + pow(xq,2)*(6 + 234*y - 73*pow(y,2))) +
                                         B*(16 + pow(xq,8) + pow(xq,6)*(23 - 6*y) + 84*y + 60*pow(y,2) - 72*pow(y,3) + pow(xq,4)*(148 - 142*y + 29*pow(y,2)) -
                                            pow(xq,2)*(116 + 208*y - 207*pow(y,2) + 24*pow(y,3))))) -
                                   rl*(8*pow(B,4)*(pow(xk,4)*pow(xq,2) - pow(xk,2)*(1 + pow(xq,2))*(pow(xq,2) - y) + pow(pow(xq,2) - y,2)) +
                                      8*pow(pow(xq,2) - y,2)*pow(-1 + y,2) + pow(xk,8)*(-1 + pow(xq,2))*(1 + pow(xq,2) - 10*y + 3*pow(xq,2)*y - pow(y,2)) +
                                      pow(xk,6)*(1 + pow(xq,6)*(2 - 4*y) - 30*y + 22*pow(y,2) + 3*pow(y,3) + pow(xq,4)*(-23 + 5*y + 8*pow(y,2)) +
                                         pow(xq,2)*(20 + 37*y - 30*pow(y,2) - 3*pow(y,3))) -
                                      pow(xk,2)*(-1 + y)*(pow(xq,8) + pow(xq,6)*(-3 - 3*y + pow(y,2)) - y*(16 + 8*y + 3*pow(y,2) + pow(y,3)) -
                                         pow(xq,4)*(30 - 21*y + 2*pow(y,2) + 2*pow(y,3)) + pow(xq,2)*(16 + 38*y - 15*pow(y,2) + 5*pow(y,3) + pow(y,4))) +
                                      pow(xk,4)*(pow(xq,8)*(-3 + y) + pow(xq,6)*(10 + 13*y - 5*pow(y,2)) + y*(37 - 32*y + 14*pow(y,2) + 3*pow(y,3)) +
                                         pow(xq,4)*(32 - 59*y + 4*pow(y,2) + 7*pow(y,3)) - pow(xq,2)*(31 + 24*y - 57*pow(y,2) + 21*pow(y,3) + 3*pow(y,4))) +
                                      pow(B,3)*(pow(xk,6)*(-7 + 15*pow(xq,2)) + 16*pow(pow(xq,2) - y,2)*(-2 + y) + pow(xk,4)*(14 - 22*pow(xq,4) + 8*pow(xq,2)*(-5 + 4*y)) +
                                         pow(xk,2)*(7*pow(xq,6) + pow(xq,4)*(47 - 24*y) + y*(-18 + 23*y) + pow(xq,2)*(18 - 70*y + 17*pow(y,2)))) +
                                      pow(B,2)*(7*pow(xk,8)*(-1 + pow(xq,2)) + 8*pow(pow(xq,2) - y,2)*(6 - 6*y + pow(y,2)) +
                                         pow(xk,6)*(22 - 19*pow(xq,4) - 21*y + pow(xq,2)*(-19 + 29*y)) -
                                         pow(xk,2)*(pow(xq,8) + pow(xq,6)*(18 - 4*y) + (43 - 12*y)*pow(y,2) - 2*pow(xq,2)*y*(68 - 43*y + 6*pow(y,2)) +
                                            pow(xq,4)*(93 - 92*y + 15*pow(y,2))) + pow(xk,4)*
                                          (-24 + 9*pow(xq,6) + pow(xq,4)*(76 - 33*y) + 31*y - 10*pow(y,2) + pow(xq,2)*(19 - 94*y + 34*pow(y,2)))) -
                                      B*(3*pow(xk,8)*(-1 + pow(xq,2))*(2 + pow(xq,2) - 2*y) + 16*pow(pow(xq,2) - y,2)*(2 - 3*y + pow(y,2)) +
                                         pow(xk,6)*(-4*pow(xq,6) + pow(xq,4)*(-32 + 5*y) + pow(xq,2)*(14 + 53*y - 14*pow(y,2)) + 14*(1 - 3*y + pow(y,2))) +
                                         pow(xk,4)*(pow(xq,8) + pow(xq,6)*(13 + 4*y) + pow(xq,4)*(86 - 74*y + pow(y,2)) + 2*(-5 + 33*y - 18*pow(y,2) + 5*pow(y,3)) -
                                            pow(xq,2)*(42 + 92*y - 83*pow(y,2) + 10*pow(y,3))) -
                                         pow(xk,2)*(pow(xq,8)*y + pow(xq,6)*(14 - 5*y - 4*pow(y,2)) - 2*y*(-13 - 6*y + 8*pow(y,2) + pow(y,3)) +
                                            pow(xq,4)*(84 - 118*y + 40*pow(y,2) + pow(y,3)) + pow(xq,2)*(-26 - 96*y + 120*pow(y,2) - 33*pow(y,3) + 2*pow(y,4))))))))/
    (pow(-1 + B + rl,2)*pow(A + 2*rl,2)*pow(xk,4)*pow(-1 + pow(xq,2),2)*pow(pow(xk,2) - pow(xq,2) + y,2));



 
    
    
  double inth1= (-8*fk*pow(MK,5)*rl*(4*pow(rl,3)*(-2 + pow(xk,2) + 2*pow(xq,2))*
                                       (-1 + pow(xk,2) + y) +
                                      2*pow(A,2)*(-1 + B + rl)*(-1 + pow(xk,2) + y)*
                                       (-1 + 2*B + pow(xk,2) + y) -
                                      pow(xk,2)*(-1 + pow(xq,2))*
                                       (-(pow(xq,4)*(-1 + B + y)) +
                                         y*(3 - 2*B - pow(B,2) + pow(xk,4) +
                                            2*pow(xk,2)*(-2 + y) - 4*y + pow(y,2)) +
                                         pow(xq,2)*(-3 + pow(B,2) + pow(xk,2) + 3*y +
                                            B*(2 - pow(xk,2) + y))) -
                                      pow(rl,2)*(8*pow(B,2)*(pow(xk,2) - pow(xq,2) + y) +
                                         pow(xk,4)*(17 - 13*pow(xq,2) + 4*y) -
                                         4*(-1 + y)*(2 + pow(xq,2)*(-3 + y) + y - pow(y,2)) +
                                         2*pow(xk,2)*(-11 + pow(xq,4) + pow(xq,2)*(10 - 8*y) +
                                            4*y + 4*pow(y,2)) +
                                         4*B*(pow(xk,4) - (-2 + 5*pow(xq,2) - 3*y)*(-1 + y) +
                                            pow(xk,2)*(3 - 7*pow(xq,2) + 4*y))) +
                                      rl*(pow(xk,6)*(-1 + pow(xq,2)) +
                                         4*(-1 + B)*(2*pow(B,2) + 3*B*(-1 + y) + pow(-1 + y,2))*
                                          (pow(xq,2) - y) -
                                         pow(xk,4)*(-13 + 8*pow(B,2) + 12*pow(xq,2) +
                                            pow(xq,4) - 4*y + B*(3 - 11*pow(xq,2) + 4*y)) -
                                         pow(xk,2)*(10 + 8*pow(B,3) - 2*pow(xq,6) -
                                            4*pow(B,2)*(2 + 5*pow(xq,2) - 5*y) - y -
                                            9*pow(y,2) + pow(xq,4)*(3 + y) +
                                            pow(xq,2)*(-15 + 12*y + pow(y,2)) +
                                            B*(-10 + pow(xq,4) + pow(xq,2)*(33 - 15*y) - 17*y +
                                               8*pow(y,2)))) +
                                      A*(pow(xk,6)*(-1 + pow(xq,2)) +
                                         2*(-1 + B)*(2*pow(B,2) + 3*B*(-1 + y) + pow(-1 + y,2))*
                                          (pow(xq,2) - y) -
                                         pow(xk,4)*(-8 + B + 4*pow(B,2) + 7*pow(xq,2) -
                                            5*B*pow(xq,2) + pow(xq,4) + 2*B*y - 2*pow(xq,2)*y) +
                                         2*pow(rl,2)*(-1 + pow(xk,2) + y)*
                                          (4*B + 3*pow(xk,2) + 2*(-2 + pow(xq,2) + y)) +
                                         pow(xk,2)*(-4*pow(B,3) +
                                            2*pow(B,2)*(3 + 4*pow(xq,2) - 5*y) +
                                            (-1 + y)*(7 - pow(xq,4) + pow(xq,2)*(-8 + y) + 3*y) +
                                            B*(5 + 2*pow(xq,4) + 9*y - 4*pow(y,2) +
                                               pow(xq,2)*(-19 + 7*y))) +
                                         rl*(2*B*(pow(xk,4) + 2*pow(xk,2)*(-5 + 3*pow(xq,2)) +
                                               (-8 + 5*pow(xq,2) - y)*(-1 + y)) +
                                            4*pow(B,2)*(-2 + pow(xk,2) + pow(xq,2) + y) -
                                            (-1 + pow(xk,2) + y)*
                                             (pow(xk,2)*(13 - 7*pow(xq,2) + 2*y) +
                                               2*(-4 + 3*pow(xq,2) + y - pow(xq,2)*y + pow(y,2)))
                                            ))))/
                                  ((-1 + B + rl)*(A + 2*rl)*(-1 + pow(xq,2))*
                                   (pow(xk,2) - pow(xq,2) + y));



  double inth2 =(4*fk*pow(MK,5)*rl*(1 - pow(xq,2))*
                 (4*pow(rl,4)*pow(xk,2)*(-1 + pow(xq,2)) +
                   2*pow(A,3)*(-1 + B + rl)*(-1 + 2*rl - pow(xq,2))*
                    (-1 + pow(xk,2) + y) +
                   pow(rl,3)*(4*(1 - 2*B + pow(xq,2) - 2*y)*(pow(xq,2) - y)*
                       (-1 + B + y) + pow(xk,4)*(1 + 8*B - pow(xq,2) + 8*y) +
                      2*pow(xk,2)*(1 + 4*pow(B,2) + 2*pow(xq,2) +
                         pow(xq,4) - 2*B*(1 + 5*pow(xq,2) - 6*y) - 2*y -
                         10*pow(xq,2)*y + 8*pow(y,2))) +
                   pow(xk,2)*(-1 + pow(xq,2))*
                    (pow(y,2)*pow(-1 + pow(xk,2) + y,2) +
                      pow(B,3)*(-pow(xq,2) + y) -
                      pow(xq,2)*y*(-1 + pow(xk,2) + y)*(-3 + 2*y) +
                      pow(xq,4)*(2 - 3*y + pow(y,2)) +
                      pow(B,2)*(pow(xq,4) + y*(-2 + 2*pow(xk,2) + 3*y) -
                         pow(xq,2)*(-2 + pow(xk,2) + 4*y)) +
                      B*(pow(xq,4)*(-3 + 2*y) +
                         pow(xq,2)*(-1 + pow(xk,2) + 7*y - 3*pow(xk,2)*y -
                            5*pow(y,2)) +
                         y*(1 - 2*pow(xk,2) + pow(xk,4) - 4*y +
                            4*pow(xk,2)*y + 3*pow(y,2)))) +
                   pow(rl,2)*(-2*pow(xk,6)*(-1 + pow(xq,2)) -
                      2*(-3 + 2*B - pow(xq,2))*(pow(xq,2) - y)*(-1 + B + y)*
                       (-1 + 2*B - pow(xq,2) + 2*y) +
                      pow(xk,4)*(-7 + 8*pow(B,2) + pow(xq,2) +
                         6*pow(xq,4) - 2*y - 14*pow(xq,2)*y +
                         8*B*(-2*pow(xq,2) + y)) +
                      pow(xk,2)*(4 + 8*pow(B,3) - 4*pow(xq,6) -
                         8*pow(B,2)*(1 + 4*pow(xq,2) - 3*y) + y -
                         17*pow(y,2) + pow(xq,4)*(-11 + 19*y) +
                         pow(xq,2)*(-5 + 28*y - 15*pow(y,2)) +
                         B*(-6 + 21*pow(xq,4) + pow(xq,2)*(41 - 43*y) - 29*y +
                            16*pow(y,2)))) -
                   pow(A,2)*(1 - pow(xk,2) - pow(xk,4) + pow(xk,6) +
                      4*pow(xq,2) - 5*pow(xk,2)*pow(xq,2) +
                      4*pow(xk,4)*pow(xq,2) - pow(xk,6)*pow(xq,2) +
                      3*pow(xq,4) - 6*pow(xk,2)*pow(xq,4) +
                      pow(xk,4)*pow(xq,4) +
                      rl*(-4*pow(xk,4)*(1 + pow(xq,2)) +
                         pow(xk,2)*(1 + 5*pow(xq,4) + pow(xq,2)*(10 - 8*y) -
                            16*y) + (-1 + 3*pow(xq,4) + pow(xq,2)*(6 - 4*y) -
                            12*y)*(-1 + y)) - 5*y + 3*pow(xk,2)*y +
                      2*pow(xk,4)*y - 8*pow(xq,2)*y +
                      8*pow(xk,2)*pow(xq,2)*y - 2*pow(xk,4)*pow(xq,2)*y -
                      3*pow(xq,4)*y + pow(xk,2)*pow(xq,4)*y + 4*pow(y,2) +
                      pow(xk,2)*pow(y,2) + 4*pow(xq,2)*pow(y,2) -
                      pow(xk,2)*pow(xq,2)*pow(y,2) -
                      8*pow(rl,3)*(-1 + pow(xk,2) + y) +
                      2*pow(B,2)*(-(pow(xk,2)*(1 + 3*pow(xq,2))) +
                         (1 + pow(xq,2))*(1 + pow(xq,2) - 2*y) +
                         rl*(-2 + 4*pow(xk,2) - 2*pow(xq,2) + 4*y)) +
                      2*pow(rl,2)*(2*pow(xk,4) -
                         (-5 + pow(xq,2) - 4*y)*(-1 + y) +
                         pow(xk,2)*(4 - 2*pow(xq,2) + 6*y)) +
                      B*(-3 - 8*pow(xq,2) - 4*pow(xk,4)*pow(xq,2) -
                         5*pow(xq,4) - 4*pow(rl,2)*(-1 + pow(xq,2)) + 9*y +
                         12*pow(xq,2)*y + 3*pow(xq,4)*y - 4*pow(y,2) -
                         4*pow(xq,2)*pow(y,2) +
                         pow(xk,2)*(3 + 12*pow(xq,2) + 5*pow(xq,4) - 4*y -
                            8*pow(xq,2)*y) +
                         2*rl*(2 + 2*pow(xk,4) + 5*pow(xq,2) + pow(xq,4) -
                            9*y - 3*pow(xq,2)*y + 4*pow(y,2) +
                            pow(xk,2)*(-5 - 5*pow(xq,2) + 6*y)))) +
                   rl*(2*(-1 + B)*(1 + pow(xq,2))*(pow(xq,2) - y)*(-1 + B + y)*
                       (-1 + 2*B - pow(xq,2) + 2*y) -
                      pow(xk,6)*(-1 + pow(xq,2))*(-2 + 5*B + 3*y) +
                      pow(xk,4)*(4 + pow(B,2)*(5 - 13*pow(xq,2)) - 6*y +
                         4*pow(y,2) + pow(xq,4)*(-8 + 3*y) +
                         pow(xq,2)*(4 + 11*y - 4*pow(y,2)) +
                         B*(-11 + 9*pow(xq,4) + 8*y - 2*pow(xq,2)*(-5 + 8*y)))\
                       + pow(xk,2)*(-2 - 8*pow(B,3)*pow(xq,2) +
                         4*pow(xq,6) +
                         4*pow(B,2)*(-1 + 4*pow(xq,4) -
                            5*pow(xq,2)*(-1 + y) - y) + y + 2*pow(y,2) +
                         pow(y,3) + pow(xq,4)*(7 - 18*y + pow(y,2)) -
                         pow(xq,2)*(1 + 7*y - 13*pow(y,2) + pow(y,3)) +
                         B*(6 - 4*pow(xq,6) + y - 3*pow(y,2) +
                            pow(xq,4)*(-25 + 17*y) +
                            pow(xq,2)*(-9 + 30*y - 13*pow(y,2))))) -
                   A*(pow(xq,2) + 2*pow(xk,2)*pow(xq,2) -
                      3*pow(xk,4)*pow(xq,2) + 2*pow(xq,4) -
                      3*pow(xk,2)*pow(xq,4) + 3*pow(xk,4)*pow(xq,4) +
                      pow(xq,6) - 3*pow(xk,2)*pow(xq,6) - y +
                      3*pow(xk,4)*y - 2*pow(xk,6)*y - 5*pow(xq,2)*y -
                      4*pow(xk,4)*pow(xq,2)*y + 2*pow(xk,6)*pow(xq,2)*y -
                      5*pow(xq,4)*y + 11*pow(xk,2)*pow(xq,4)*y -
                      3*pow(xk,4)*pow(xq,4)*y - pow(xq,6)*y +
                      pow(xk,2)*pow(xq,6)*y + 3*pow(y,2) +
                      pow(xk,2)*pow(y,2) - 4*pow(xk,4)*pow(y,2) +
                      6*pow(xq,2)*pow(y,2) -
                      6*pow(xk,2)*pow(xq,2)*pow(y,2) +
                      4*pow(xk,4)*pow(xq,2)*pow(y,2) +
                      3*pow(xq,4)*pow(y,2) -
                      3*pow(xk,2)*pow(xq,4)*pow(y,2) - 2*pow(y,3) -
                      2*pow(xk,2)*pow(y,3) - 2*pow(xq,2)*pow(y,3) +
                      2*pow(xk,2)*pow(xq,2)*pow(y,3) -
                      2*pow(B,3)*(pow(xq,4) - y -
                         pow(xq,2)*(-1 + 2*pow(xk,2) + y)) +
                      2*pow(rl,3)*(4*pow(xk,4) -
                         2*(1 + 3*pow(xq,2) - 4*y)*(-1 + y) +
                         B*(-4 + 8*pow(xk,2) - 4*pow(xq,2) + 8*y) +
                         pow(xk,2)*(-1 - 11*pow(xq,2) + 12*y)) +
                      pow(rl,2)*(4*pow(B,2)*
                          (-2 + 3*pow(xk,2) - pow(xq,2) + 3*y) -
                         pow(xk,4)*(3 + 13*pow(xq,2) + 4*y) +
                         pow(xk,2)*(2 + 15*pow(xq,4) +
                            pow(xq,2)*(27 - 11*y) - 25*y - 8*pow(y,2)) -
                         2*(-1 + y)*(-3 - 2*pow(xq,4) + pow(xq,2)*(-9 + y) +
                            11*y + 2*pow(y,2)) +
                         2*B*(8 + 2*pow(xk,4) + 11*pow(xq,2) + pow(xq,4) -
                            2*pow(xk,2)*(5 + 6*pow(xq,2) - 3*y) - 19*y -
                            5*pow(xq,2)*y + 4*pow(y,2))) +
                      pow(B,2)*(pow(xk,4)*(-1 + 5*pow(xq,2)) +
                         (1 + pow(xq,2))*
                          (pow(xq,4) - 5*pow(xq,2)*(-1 + y) + y*(-5 + 4*y)) +
                         pow(xk,2)*(1 - 8*pow(xq,4) + y +
                            pow(xq,2)*(-9 + 11*y))) +
                      B*(pow(xk,6)*(-1 + pow(xq,2)) + pow(xq,6)*(-2 + y) -
                         2*pow(xk,4)*(-1 + pow(xq,4) + pow(xq,2)*(2 - 4*y) +
                            2*y) + pow(xq,4)*(-6 + 10*y - 3*pow(y,2)) +
                         y*(4 - 7*y + 2*pow(y,2)) +
                         pow(xq,2)*(-4 + 13*y - 10*pow(y,2) + 2*pow(y,3)) +
                         pow(xk,2)*(-1 + pow(xq,6) - pow(y,2) -
                            2*pow(xq,4)*(-7 + 5*y) +
                            pow(xq,2)*(2 - 14*y + 9*pow(y,2)))) +
                      rl*(-3*pow(xk,6)*(-1 + pow(xq,2)) +
                         pow(B,2)*(-4*pow(xk,4) +
                            2*pow(xk,2)*(1 + pow(xq,2) - 6*y) +
                            4*(-1 + pow(xq,2) - 2*y)*(-1 + y)) -
                         4*pow(B,3)*(pow(xk,2) - pow(xq,2) + y) +
                         pow(xk,4)*(-5 + 10*pow(xq,2) + 3*pow(xq,4) + 8*y) +
                         (-1 + y)*(-2 + pow(xq,6) + 5*y + 6*pow(y,2) -
                            pow(xq,4)*(2 + 3*y) +
                            pow(xq,2)*(-5 - 2*y + 2*pow(y,2))) +
                         pow(xk,2)*(-(pow(xq,4)*(9 + 5*y)) + y*(1 + 11*y) +
                            pow(xq,2)*(-7 + 4*y + 5*pow(y,2))) +
                         B*(-6 + pow(xq,6) +
                            pow(xk,4)*(6 - 6*pow(xq,2) - 4*y) + 7*y +
                            10*pow(y,2) - 4*pow(y,3) - pow(xq,4)*(2 + y) +
                            pow(xq,2)*(-5 - 2*y + 2*pow(y,2)) +
                            pow(xk,2)*(3 + 5*pow(xq,4) + 10*y - 8*pow(y,2) +
                               2*pow(xq,2)*(2 + y)))))))/
               ((-1 + B + rl)*(A + 2*rl)*pow(-1 + pow(xq,2),2)*
                (pow(xk,2) - pow(xq,2) + y));

 
  
 double intfA =  (4*fk*pow(MK,5)*rl*(4*pow(A,3)*(-1 + B + rl)*
                                       pow(-1 + pow(xk,2) + y,2) +
                                      4*pow(rl,3)*(pow(xk,6) +
                                         2*pow(-1 + pow(xq,2),2)*(-1 + y) +
                                         pow(xk,4)*(-4 + 3*pow(xq,2) + y) +
                                         2*pow(xk,2)*(-1 + pow(xq,2))*(-3 + pow(xq,2) + 2*y)) +
                                      pow(rl,2)*(8*pow(B,2)*
                                          (pow(pow(xq,2) - y,2) + pow(xk,2)*(-1 + y)) +
                                         4*pow(xk,6)*(-4 + 3*pow(xq,2) + y) +
                                         pow(xk,4)*(35 + 12*pow(xq,4) - 23*y + 16*pow(y,2) -
                                            pow(xq,2)*(39 + y)) +
                                         8*(-1 + y)*(-1 + pow(xq,4)*(-2 + y) - pow(y,2) +
                                            pow(y,3) + pow(xq,2)*(2 + 2*y - 2*pow(y,2))) +
                                         2*B*(2*pow(xk,6) +
                                            pow(xk,4)*(-15 + 7*pow(xq,2) + 8*y) +
                                            4*(-1 + y)*(1 - 2*pow(xq,2) + 3*pow(xq,4) -
                                               4*pow(xq,2)*y + 2*pow(y,2)) +
                                            2*pow(xk,2)*(8 - 8*pow(xq,2) + 7*pow(xq,4) - 8*y -
                                               6*pow(xq,2)*y + 7*pow(y,2))) -
                                         2*pow(xk,2)*(15 + pow(xq,6) + pow(xq,4)*(11 - 10*y) -
                                            12*y + 10*pow(y,2) - 10*pow(y,3) +
                                            pow(xq,2)*(-21 + 16*pow(y,2)))) +
                                      pow(A,2)*(-(pow(xk,6)*(-3 + 2*B + pow(xq,2))) +
                                         8*(-1 + B)*(pow(xq,2) - y)*(-1 + y)*(-1 + B + y) +
                                         8*pow(rl,2)*pow(-1 + pow(xk,2) + y,2) -
                                         2*rl*pow(-1 + pow(xk,2) + y,2)*
                                          (4 + pow(xk,2) - 4*pow(xq,2) + 4*y) +
                                         pow(xk,4)*(-5 - 4*pow(B,2) + pow(xq,4) +
                                            4*B*(2 + 2*pow(xq,2) - 3*y) + 14*y -
                                            2*pow(xq,2)*(4 + y)) +
                                         4*B*rl*(pow(xk,4) + 2*(-1 + pow(xq,2))*(-1 + y) +
                                            pow(xk,2)*(-2 + pow(xq,2) + y)) +
                                         pow(xk,2)*(4*pow(B,2)*(2 + pow(xq,2) - 3*y) +
                                            (-1 + y)*(-2 + pow(xq,4) + 19*y -
                                               pow(xq,2)*(17 + y)) +
                                            2*B*(-5 + 16*y - 9*pow(y,2) + 2*pow(xq,2)*(-5 + 4*y))))
                                        - pow(xk,2)*(-1 + pow(xq,2))*
                                       (-2*pow(B,3)*pow(xk,2) - pow(xq,6)*(-1 + y) +
                                         pow(xq,4)*(-(pow(xk,2)*(-2 + y)) + 3*(-1 + y)*y) +
                                         pow(B,2)*(-2*pow(xk,4) + pow(xq,4) +
                                            pow(xq,2)*(3 - 5*y) +
                                            pow(xk,2)*(2 + 3*pow(xq,2) - y) + y*(-3 + 4*y)) +
                                         pow(xq,2)*(3 - 8*y + 10*pow(y,2) - 5*pow(y,3) +
                                            pow(xk,4)*(1 + y) + pow(xk,2)*(-4 + 6*y - 4*pow(y,2))
                                            ) + y*(-3 + pow(xk,6) + 5*pow(xk,4)*(-1 + y) + 8*y -
                                            8*pow(y,2) + 3*pow(y,3) +
                                            pow(xk,2)*(7 - 13*y + 7*pow(y,2))) -
                                         B*(pow(xq,6) + pow(xq,4)*(3 + 2*pow(xk,2) - 6*y) +
                                            pow(xq,2)*(6 - 5*pow(xk,2) + pow(xk,4) - 15*y +
                                               6*pow(xk,2)*y + 12*pow(y,2)) -
                                            y*(6 + 2*pow(xk,4) - 12*y + 7*pow(y,2) +
                                               pow(xk,2)*(-8 + 9*y)))) +
                                      rl*(pow(xk,8)*(-1 + pow(xq,2)) +
                                         8*(-1 + B)*pow(pow(xq,2) - y,2)*pow(-1 + B + y,2) +
                                         pow(xk,6)*(12 + 2*pow(xq,2)*(-6 + y) - 6*y +
                                            B*(-13 + 13*pow(xq,2) + 4*y)) +
                                         pow(xk,4)*(-17 - 15*pow(xq,4) + pow(xq,6) + 12*y -
                                            12*pow(y,2) + 4*pow(B,2)*(-4 + pow(xq,2) + 3*y) +
                                            pow(xq,2)*(27 + 8*y - 4*pow(y,2)) +
                                            B*(31 + 9*pow(xq,4) - 30*y + 16*pow(y,2) -
                                               2*pow(xq,2)*(12 + y))) +
                                         pow(xk,2)*(2*pow(xq,8) + 8*pow(B,3)*(-1 + y) -
                                            2*pow(xq,6)*(1 + 2*y) +
                                            pow(xq,4)*(21 - 22*y + 7*pow(y,2)) +
                                            4*pow(B,2)*(4 + 5*pow(xq,4) - 6*y - 10*pow(xq,2)*y +
                                               7*pow(y,2)) -
                                            pow(xq,2)*(15 + 15*y - 33*pow(y,2) + 5*pow(y,3)) -
                                            3*(-2 + y - 4*pow(y,2) + 5*pow(y,3)) +
                                            B*(-14 - 2*pow(xq,6) + 19*y - 43*pow(y,2) +
                                               20*pow(y,3) + pow(xq,4)*(-39 + 25*y) +
                                               pow(xq,2)*(15 + 56*y - 37*pow(y,2))))) +
                                      A*(pow(xk,8)*(-1 + pow(xq,2)) +
                                         4*(-1 + B)*pow(pow(xq,2) - y,2)*pow(-1 + B + y,2) +
                                         pow(xk,6)*(7 - 8*y + B*(-3 + 3*pow(xq,2) + 2*y) +
                                            pow(xq,2)*(-7 + 6*y)) -
                                         pow(xk,4)*(11 + pow(xq,6) +
                                            pow(B,2)*(2 + 4*pow(xq,2) - 6*y) - 20*y +
                                            17*pow(y,2) + pow(xq,4)*(4 + 6*y) +
                                            pow(xq,2)*(-14 + 4*y - 9*pow(y,2)) -
                                            2*B*(6 + 5*pow(xq,4) + pow(xq,2)*(-7 + y) - 9*y +
                                               4*pow(y,2))) -
                                         2*pow(rl,2)*(pow(xk,6) +
                                            pow(xk,4)*(-1 + 4*B - 10*pow(xq,2) + 11*y) -
                                            2*pow(xk,2)*(3 - 13*pow(xq,2) + pow(xq,4) +
                                               2*B*(2 + pow(xq,2) - 3*y) + 7*y + 11*pow(xq,2)*y -
                                               9*pow(y,2)) -
                                            2*(-1 + y)*(1 + pow(xq,4) + 4*y - 4*B*y - 4*pow(y,2) +
                                               pow(xq,2)*(-6 + 4*B + 4*y))) +
                                         pow(xk,2)*(4*pow(B,3)*(-1 + y) +
                                            2*pow(B,2)*(4 + pow(xq,2) + 4*pow(xq,4) - 7*y -
                                               9*pow(xq,2)*y + 7*pow(y,2)) +
                                            B*(-9 + pow(xq,6) + 21*y - 31*pow(y,2) +
                                               10*pow(y,3) + pow(xq,4)*(-14 + 3*y) +
                                               pow(xq,2)*(2 + 26*y - 9*pow(y,2))) +
                                            (-1 + y)*(-5 + 2*pow(xq,6) + 7*y - 14*pow(y,2) -
                                               6*pow(xq,4)*(1 + y) +
                                               pow(xq,2)*(3 + 15*y + 4*pow(y,2)))) +
                                         rl*(pow(xk,6)*(-1 + 3*pow(xq,2) + 2*y) +
                                            pow(xk,4)*(9 + 11*pow(xq,4) +
                                               2*pow(xq,2)*(-20 + y) + 10*y + 8*pow(y,2)) -
                                            4*pow(B,2)*(2*pow(xk,4) - pow(xq,4) -
                                               2*pow(xq,2)*(-2 + y) + y*(-4 + 3*y) +
                                               pow(xk,2)*(-3 - 2*pow(xq,2) + 5*y)) +
                                            pow(xk,2)*(-17 + 4*pow(xq,4)*(-2 + y) - 16*y +
                                               19*pow(y,2) + 10*pow(y,3) +
                                               pow(xq,2)*(55 - 38*y - 9*pow(y,2))) +
                                            4*(-1 + y)*(-1 + pow(xq,4)*(-2 + y) - 4*y +
                                               3*pow(y,2) + pow(y,3) -
                                               2*pow(xq,2)*(-3 + y + pow(y,2))) -
                                            2*B*(pow(xk,6) +
                                               pow(xk,4)*(-4 - 8*pow(xq,2) + 8*y) +
                                               2*(-1 + y)*(-1 + 10*pow(xq,2) - 3*pow(xq,4) - 8*y +
                                                  2*pow(y,2)) +
                                               pow(xk,2)*(-6*pow(xq,4) + pow(xq,2)*(29 - 13*y) +
                                                  y*(-21 + 11*y)))))))/
                                  ((-1 + B + rl)*(A + 2*rl)*pow(xk,2)*(-1 + pow(xq,2))*
                                   (pow(xk,2) - pow(xq,2) + y));


 double intfV=  (4*fk*pow(MK,5)*rl*(16*pow(rl,4)*pow(xk,2) +
                                      pow(A,2)*pow(xk,2)*
                                       (-4*pow(B,2) - (-2 + 2*rl - pow(xk,2) + pow(xq,2) - y)*
                                          (-1 + pow(xk,2) + y) - 2*B*(-3 + 2*rl + pow(xk,2) + y))\
                                       + 8*pow(rl,3)*(2*pow(xk,4) + (-1 + pow(xq,2))*(-1 + y) +
                                         pow(xk,2)*(-4 + 2*B - pow(xq,2) + y)) +
                                      pow(rl,2)*(8*pow(B,2)*pow(xk,2) + 4*pow(xk,6) +
                                         2*B*(9*pow(xk,4) -
                                            2*pow(xk,2)*(6 + 3*pow(xq,2) - 5*y) +
                                            4*(-1 + pow(xq,2))*(-1 + y)) +
                                         pow(xk,4)*(-23 - 4*pow(xq,2) + 3*y) +
                                         8*(-1 + pow(xq,2) + y - pow(xq,2)*y) +
                                         2*pow(xk,2)*(11 + 2*pow(xq,2) + pow(xq,4) - 8*y +
                                            2*pow(y,2))) +
                                      rl*pow(xk,2)*(-6 + 8*pow(B,3) + pow(xk,6) +
                                         5*pow(xq,2) - 4*pow(xq,4) - 2*pow(xq,6) + 7*y -
                                         2*pow(xq,2)*y + 8*pow(xq,4)*y -
                                         5*pow(xq,2)*pow(y,2) - pow(y,3) +
                                         4*pow(B,2)*(-4 + pow(xk,2) - pow(xq,2) + 3*y) -
                                         2*pow(xk,4)*(2 + pow(xq,2) + 3*y) +
                                         pow(xk,2)*(9 + 3*pow(xq,4) + 2*pow(xq,2)*(-1 + y) +
                                            4*y - 8*pow(y,2)) +
                                         B*(14 + pow(xk,4) + 2*pow(xq,4) +
                                            pow(xk,2)*(-11 + pow(xq,2) - 2*y) +
                                            pow(xq,2)*(-1 + y) - 19*y + 3*pow(y,2))) +
                                      pow(xk,2)*(2*pow(B,3)*pow(xk,2) - pow(xq,6)*(-1 + y) +
                                         pow(xq,4)*(4 + (-9 + pow(xk,2))*y + 5*pow(y,2)) +
                                         pow(xq,2)*(-1 + y)*
                                          (3 + pow(xk,4) + 3*y - 5*pow(y,2) -
                                            4*pow(xk,2)*(1 + y)) +
                                         pow(B,2)*(2*pow(xk,4) +
                                            3*(-1 + pow(xq,2))*(pow(xq,2) - y) +
                                            pow(xk,2)*(-2 - 3*pow(xq,2) + y)) +
                                         y*(3 - pow(xk,6) - pow(xk,4)*(-5 + y) - 4*y +
                                            pow(y,3) + pow(xk,2)*(-7 + 5*y + pow(y,2))) +
                                         B*(-pow(xq,6) + pow(xq,4)*(-5 + 6*y) +
                                            pow(xq,2)*(6 - 5*pow(xk,2) + pow(xk,4) + y -
                                               6*pow(y,2)) +
                                            y*(-6 + 8*pow(xk,2) - 2*pow(xk,4) + 4*y -
                                               pow(xk,2)*y + pow(y,2)))) +
                                      A*(8*pow(rl,3)*pow(xk,2) +
                                         2*pow(rl,2)*(5*pow(xk,4) -
                                            2*pow(xk,2)*(3 + pow(xq,2)) +
                                            2*(-1 + pow(xq,2))*(-1 + y)) +
                                         pow(xk,2)*(4*pow(B,3) + pow(xk,6) -
                                            pow(xk,4)*(7 + 2*pow(xq,2)) +
                                            2*pow(B,2)*(-4 + pow(xk,2) - pow(xq,2) + 3*y) +
                                            pow(xk,2)*(11 + 7*pow(xq,2) + pow(xq,4) - 8*y +
                                               2*pow(xq,2)*y - 3*pow(y,2)) +
                                            B*(9 + 3*pow(xk,4) - pow(xq,2) - pow(xq,4) -
                                               4*pow(xk,2)*(3 + pow(xq,2) - y) - 13*y +
                                               5*pow(xq,2)*y + pow(y,2)) -
                                            (-1 + y)*(-5 + 2*pow(xq,2) + 2*pow(xq,4) + 3*y -
                                               4*pow(xq,2)*y + 2*pow(y,2))) +
                                         rl*(-4*pow(B,2)*pow(xk,2) + 9*pow(xk,6) +
                                            2*B*(6*pow(xk,4) - 3*pow(xk,2)*(pow(xq,2) - y) +
                                               2*(-1 + pow(xq,2))*(-1 + y)) +
                                            pow(xk,4)*(-21 - 11*pow(xq,2) + 12*y) +
                                            4*(-1 + pow(xq,2) + y - pow(xq,2)*y) +
                                            pow(xk,2)*(9 - 8*y + 3*pow(y,2) +
                                               2*pow(xq,2)*(1 + y))))))/
                                  ((-1 + B + rl)*(A + 2*rl)*pow(xk,2)*
                                   (pow(xk,2) - pow(xq,2) + y));



 double inthe1 =  (8*fk*pow(MK,5)*rl*(-(pow(A,2)*(-1 + B + rl)*(-1 + pow(xk,2) + y)*(-3 + 2*B + 4*rl + pow(xk,2) - 2*pow(xq,2) + y)) +
                                       A*(-2*pow(xk,6)*(-1 + pow(xq,2)) -
                                          8*pow(rl,3)*(-1 + pow(xk,2) + y) -
                                          (-1 + B)*(pow(xq,2) - y)*(-1 + B + y)*
                                           (-3 + 2*B - 2*pow(xq,2) + y) +
                                          pow(xk,4)*(-8 + 2*pow(B,2) + 2*pow(xq,4) -
                                             4*pow(xq,2)*(-2 + y) + 3*y + B*(4 - 8*pow(xq,2) + y))\
                                           - pow(rl,2)*(pow(xk,4) +
                                             pow(xk,2)*(-19 + 4*pow(xq,2) - y) +
                                             2*(-8 + pow(xq,2) - y)*(-1 + y) +
                                             4*B*(-3 + 2*pow(xk,2) + pow(xq,2) + 2*y)) +
                                          pow(xk,2)*(6 + 2*pow(B,3) + 2*pow(xq,4)*(-3 + y) -
                                             4*y + pow(B,2)*(-11*pow(xq,2) + 5*y) +
                                             pow(xq,2)*(-7 + 11*y - 2*pow(y,2)) +
                                             B*(-8 + 4*pow(xq,4) + pow(xq,2)*(20 - 11*y) - 3*y +
                                                2*pow(y,2))) +
                                          rl*(pow(xk,4)*(7 - 8*pow(xq,2) + y) +
                                             2*pow(B,2)*(2 + pow(xk,2) - 3*pow(xq,2) + y) +
                                             (-1 + y)*(-8 + 5*pow(xq,2) + 2*pow(xq,4) - 5*y -
                                                3*pow(xq,2)*y + pow(y,2)) +
                                             pow(xk,2)*(-17 + 13*pow(xq,2) + 4*pow(xq,4) + y -
                                                11*pow(xq,2)*y + 2*pow(y,2)) +
                                             B*(-12 + pow(xk,4) + 11*pow(xq,2) + 2*pow(xq,4) +
                                                y - 7*pow(xq,2)*y + 5*pow(y,2) +
                                                pow(xk,2)*(13 - 15*pow(xq,2) + 6*y)))) +
                                       2*(pow(xk,2)*(-1 + pow(xq,2))*
                                           (-(pow(xq,4)*(-1 + B + y)) -
                                             (-1 + B)*pow(xq,2)*(-1 + B + pow(xk,2) + y) +
                                             y*pow(-1 + B + pow(xk,2) + y,2)) +
                                          pow(rl,3)*(pow(xk,4) +
                                             pow(xk,2)*(7 + 4*B - 12*pow(xq,2) + 5*y) -
                                             2*(1 + y - 2*B*y - 2*pow(y,2) +
                                                pow(xq,2)*(-3 + 2*B + 3*y))) +
                                          pow(rl,2)*(pow(xk,4)*(6 - 9*pow(xq,2) + y) +
                                             6*pow(B,2)*(pow(xk,2) - pow(xq,2) + y) +
                                             (-1 + y)*(-2 + 9*pow(xq,2) + 2*pow(xq,4) - 7*y -
                                                3*pow(xq,2)*y + pow(y,2)) +
                                             pow(xk,2)*(-13 + 15*pow(xq,2) + 10*pow(xq,4) - y -
                                                13*pow(xq,2)*y + 2*pow(y,2)) +
                                             B*(-2 + 3*pow(xk,4) + 15*pow(xq,2) + 2*pow(xq,4) -
                                                11*y - 11*pow(xq,2)*y + 7*pow(y,2) +
                                                pow(xk,2)*(5 - 23*pow(xq,2) + 10*y))) +
                                          rl*(-(pow(xk,6)*(-1 + pow(xq,2))) -
                                             (-1 + B)*(pow(xq,2) - y)*(-1 + B + y)*
                                              (-3 + 2*B - 2*pow(xq,2) + y) +
                                             pow(xk,4)*(-7 + 2*pow(B,2) + pow(xq,4) -
                                                2*pow(xq,2)*(-4 + y) + y + B*(3 - 7*pow(xq,2) + y))\
                                              + pow(xk,2)*(6 + 2*pow(B,3) - 2*pow(xq,6) - 3*y -
                                                pow(y,2) + pow(xq,4)*(-7 + 3*y) +
                                                pow(B,2)*(-11*pow(xq,2) + 5*y) -
                                                pow(xq,2)*(4 - 9*y + pow(y,2)) +
                                                B*(-8 + 17*pow(xq,2) + 7*pow(xq,4) - 4*y -
                                                   10*pow(xq,2)*y + 2*pow(y,2)))))))/
                                   ((-1 + B + rl)*pow(xk,2)*(-1 + pow(xq,2))*
                                    (pow(xk,2) - pow(xq,2) + y));


 double inthe2= (4*fk*pow(MK,5)*(1 - B - rl)*rl*
                 (pow(A,3)*(-1 + B + rl)*pow(-1 + pow(xk,2) + y,2) +
                   (1 + B)*pow(xk,2)*(-1 + pow(xq,2))*(-1 + B + 2*y)*
                    (pow(xq,2)*(-1 + B + y) - y*(-1 + B + pow(xk,2) + y)) -
                   4*pow(rl,4)*(2 - 3*y + pow(xk,2)*y + pow(y,2) +
                      pow(xq,2)*(-1 - pow(xk,2) + y) +
                      B*(pow(xk,2) - pow(xq,2) + y)) -
                   pow(A,2)*(-1 + B + rl)*
                    ((-1 + pow(xk,2) + y)*
                       (-2*pow(rl,2) +
                         rl*(5 - 5*pow(xk,2) + pow(xq,2) - 4*y) -
                         (pow(xq,2) - 2*y)*(-1 + y) +
                         pow(xk,2)*(1 - pow(xq,2) + 2*y)) +
                      B*(pow(xk,4) - (-1 + pow(xq,2) - 2*y)*(-1 + y) +
                         pow(xk,2)*(-2 + pow(xq,2) + 3*y))) -
                   pow(rl,3)*(4*pow(B,2)*(pow(xk,2) - pow(xq,2) + y) +
                      pow(xk,4)*(3 - 3*pow(xq,2) + 6*y) -
                      2*(-1 + y)*(-6 + pow(xq,2) + pow(xq,4) + 7*y -
                         pow(xq,2)*y - 2*pow(y,2)) +
                      2*pow(xk,2)*(-1 + 2*pow(xq,2) + 3*pow(xq,4) - 6*y -
                         3*pow(xq,2)*y + 5*pow(y,2)) +
                      2*B*(4 + 3*pow(xk,4) + 3*pow(xq,2) + pow(xq,4) - 11*y -
                         pow(xq,2)*y + 4*pow(y,2) +
                         pow(xk,2)*(-7 - 4*pow(xq,2) + 7*y))) -
                   pow(rl,2)*(2*pow(B,2)*
                       (3*pow(xk,4) + pow(xq,4) -
                         2*pow(xk,2)*(3 + pow(xq,2) - 2*y) -
                         2*pow(xq,2)*(-2 + y) + (-4 + y)*y) -
                      2*(-1 + y)*(2 + pow(xq,2) - pow(xq,4) - 5*y +
                         2*pow(xq,2)*y + pow(y,2) - pow(xq,2)*pow(y,2) +
                         pow(y,3)) + pow(xk,4)*
                       (-2 - 7*y - 2*pow(y,2) + pow(xq,2)*(2 + y)) +
                      pow(xk,2)*(2 + 7*y - 5*pow(y,2) - 4*pow(y,3) +
                         pow(xq,4)*(-9 + 5*y) + pow(xq,2)*(3 + pow(y,2))) +
                      B*(pow(xk,4)*(-2 - 4*pow(xq,2) + 4*y) +
                         pow(xk,2)*(6 + 9*pow(xq,4) + pow(xq,2)*(5 - 3*y) -
                            17*y + 4*pow(y,2)) -
                         2*(2 - 11*y + pow(xq,4)*y + 8*pow(y,2) +
                            pow(xq,2)*(5 - 2*y - 3*pow(y,2))))) +
                   rl*(2*(pow(xq,2) - y)*pow(-1 + y,2)*y +
                      2*pow(B,3)*(1 + y)*(pow(xk,2) - pow(xq,2) + y) +
                      pow(xk,4)*(1 - 4*pow(y,2) +
                         pow(xq,2)*(-1 + 2*pow(y,2))) +
                      pow(B,2)*(pow(xk,4)*(3 - 3*pow(xq,2) + 2*y) -
                         4*(pow(xq,2) - y)*(-1 - y + pow(y,2)) +
                         2*pow(xk,2)*(-3*pow(xq,2) + pow(xq,4) + y -
                            4*pow(xq,2)*y + 3*pow(y,2))) -
                      2*pow(xk,2)*(pow(xq,4)*(2 - 4*y + pow(y,2)) +
                         y*(1 - 4*y + 3*pow(y,2)) -
                         pow(xq,2)*(2 - 4*y + pow(y,3))) -
                      2*B*((pow(xq,2) - y)*
                          (1 + 2*y - 4*pow(y,2) + pow(y,3)) +
                         pow(xk,4)*(-2 - pow(y,2) + pow(xq,2)*(2 + y)) +
                         pow(xk,2)*(1 + pow(xq,4)*(-5 + y) - 3*y +
                            4*pow(y,2) - 2*pow(y,3) +
                            pow(xq,2)*(3 - y + 2*pow(y,2))))) +
                   A*(4*pow(rl,4)*(-1 + pow(xk,2) + y) +
                      pow(rl,3)*(6*pow(xk,4) -
                         2*(3 + 2*pow(xq,2) - y)*(-1 + y) +
                         2*B*(-2 + pow(xk,2) + pow(xq,2) + y) +
                         pow(xk,2)*(-13 - 3*pow(xq,2) + 8*y)) +
                      pow(xk,4)*(1 - 3*y - pow(y,2) +
                         pow(B,2)*(1 - pow(xq,2) + y) +
                         pow(xq,2)*(-1 + 3*y) +
                         B*(2 + pow(xq,2)*(-2 + y) - 2*y + pow(y,2))) -
                      (-1 + B)*(pow(xq,2) - y)*
                       (pow(-1 + y,2)*y + pow(B,2)*(1 + y) +
                         B*(-1 - y + 2*pow(y,2))) -
                      pow(rl,2)*(2*pow(B,2)*(pow(xk,2) - pow(xq,2) + y) +
                         pow(xk,4)*(8 - 2*pow(xq,2) + 7*y) -
                         (-1 + y)*(pow(xq,2) + pow(xq,4) + 7*y +
                            pow(xq,2)*y - 6*pow(y,2)) +
                         pow(xk,2)*(-8 - 7*y - 2*pow(xq,2)*y + 13*pow(y,2)) +
                         B*(-4 - pow(xk,4) + 3*pow(xq,2) + pow(xq,4) - 3*y -
                            pow(xq,2)*y + 4*pow(y,2) +
                            pow(xk,2)*(pow(xq,2) + 3*y))) +
                      pow(xk,2)*(pow(B,3)*(pow(xq,2) + y) -
                         pow(B,2)*(-1 + 2*pow(xq,2) + pow(xq,4) + 3*y -
                            3*pow(y,2)) -
                         (-1 + y)*(-1 + pow(xq,4) + 3*y - 4*pow(xq,2)*y +
                            2*pow(y,2)) +
                         B*(pow(xq,4)*(4 - 3*y) + pow(xq,2)*(-3 + 4*y) +
                            2*y*(1 - 3*y + pow(y,2)))) +
                      rl*(pow(xk,4)*(1 + 10*y + pow(y,2) -
                            pow(xq,2)*(1 + 3*y)) -
                         pow(B,2)*(-2 + 5*pow(xk,4) + pow(xq,4) +
                            pow(xq,2)*(6 - 4*y) - 6*y + 5*pow(y,2) +
                            pow(xk,2)*(-7 - 3*pow(xq,2) + 10*y)) +
                         pow(xk,2)*(2 + pow(xq,4)*(-1 + y) - 19*y +
                            14*pow(y,2) + 2*pow(y,3) +
                            pow(xq,2)*(3 + 2*y - 4*pow(y,2))) +
                         (-1 + y)*(2 - pow(xq,4) - 9*y + 5*pow(y,2) +
                            pow(y,3) - pow(xq,2)*(-3 + pow(y,2))) +
                         B*(pow(xk,4)*(5*pow(xq,2) - 6*y) + pow(xq,4)*y +
                            y*(-17 + 20*y - 4*pow(y,2)) -
                            pow(xq,2)*(-9 + 8*y + pow(y,2)) -
                            pow(xk,2)*(7 + 5*pow(xq,4) +
                               pow(xq,2)*(2 - 6*y) - 20*y + 10*pow(y,2)))))))/
               (pow(-1 + B + rl,2)*pow(xk,2)*(-1 + pow(xq,2))*
                (pow(xk,2) - pow(xq,2) + y));

 double intfeA = (-4*fk*pow(MK,5)*rl*(-(pow(A,3)*(-1 + B + rl)*
                                          (1 + 2*B + 2*rl - 2*pow(xk,2) - pow(xq,2) - 2*y)*
                                          (-1 + pow(xk,2) + y)) -
                                       pow(A,2)*(3 - 12*pow(xk,2) + 13*pow(xk,4) -
                                          4*pow(xk,6) + 3*pow(xq,2) - 5*pow(xk,2)*pow(xq,2) -
                                          3*pow(xk,4)*pow(xq,2) + 2*pow(xk,6)*pow(xq,2) +
                                          5*pow(xk,2)*pow(xq,4) - 2*pow(xk,4)*pow(xq,4) +
                                          2*pow(B,3)*(-1 + pow(xq,2)) - 9*y + 24*pow(xk,2)*y -
                                          13*pow(xk,4)*y - 6*pow(xq,2)*y +
                                          4*pow(xk,4)*pow(xq,2)*y - 2*pow(xk,2)*pow(xq,4)*y +
                                          11*pow(y,2) - 14*pow(xk,2)*pow(y,2) +
                                          3*pow(xq,2)*pow(y,2) +
                                          2*pow(xk,2)*pow(xq,2)*pow(y,2) - 5*pow(y,3) +
                                          14*pow(rl,3)*(-1 + pow(xk,2) + y) -
                                          pow(rl,2)*(3*pow(xk,4) +
                                             (-1 + y)*(23 + 4*pow(xq,2) + 3*y) +
                                             pow(xk,2)*(23 + pow(xq,2) + 6*y)) +
                                          pow(B,2)*(7 + 4*pow(xk,4) - pow(xq,4) - 8*y -
                                             2*pow(xq,2)*y + 4*pow(y,2) +
                                             2*rl*(-9 + 7*pow(xk,2) + 2*pow(xq,2) + 7*y) +
                                             pow(xk,2)*(-17 + 7*pow(xq,2) + 8*y)) +
                                          rl*(2*pow(xk,6) + pow(xk,4)*(-8 + 3*pow(xq,2) + 9*y) +
                                             (-1 + y)*(12 + pow(xq,2)*(7 - 3*y) - 3*y +
                                                5*pow(y,2)) +
                                             pow(xk,2)*(21 + 4*pow(xq,2) - 3*pow(xq,4) - 16*y +
                                                12*pow(y,2))) +
                                          B*(-8 + 2*pow(xk,6) - 5*pow(xq,2) + pow(xq,4) + 17*y +
                                             8*pow(xq,2)*y - 15*pow(y,2) -
                                             3*pow(xq,2)*pow(y,2) + 5*pow(y,3) +
                                             3*pow(xk,4)*(-5 + pow(xq,2) + 3*y) +
                                             2*pow(rl,2)*(-15 + 14*pow(xk,2) + pow(xq,2) +
                                                14*y) + pow(xk,2)*
                                              (29 - 4*pow(xq,2) - 3*pow(xq,4) - 30*y +
                                                12*pow(y,2)) +
                                             rl*(30 + pow(xk,4) + 4*pow(xq,2) - pow(xq,4) -
                                                28*y - 6*pow(xq,2)*y + pow(y,2) +
                                                2*pow(xk,2)*(-20 + 3*pow(xq,2) + y)))) +
                                       A*(-20*pow(rl,4)*(-1 + pow(xk,2) + y) +
                                          2*pow(xk,6)*(pow(B,2) + (-4 + 3*pow(xq,2))*y +
                                             B*(1 - 2*pow(xq,2) + y)) +
                                          pow(rl,3)*(pow(xk,4) -
                                             4*(-13 + 4*pow(xq,2) - 2*y)*(-1 + y) +
                                             pow(xk,2)*(61 - 26*pow(xq,2) + 9*y) -
                                             2*B*(-22 + 15*pow(xk,2) + 7*pow(xq,2) + 15*y)) +
                                          pow(xk,4)*(2 + 4*pow(B,3) + pow(xq,4)*(6 - 8*y) -
                                             2*pow(B,2)*(1 + 7*pow(xq,2) - 5*y) + 19*y -
                                             19*pow(y,2) + pow(xq,2)*(-10 - y + 12*pow(y,2)) +
                                             B*(-4 + 2*pow(xq,4) - 21*y + 7*pow(y,2) +
                                                pow(xq,2)*(16 + y))) +
                                          pow(rl,2)*(-4*pow(xk,6) +
                                             pow(xk,4)*(30 - 23*pow(xq,2) - 16*y) +
                                             2*pow(B,2)*(14 + pow(xk,2) - 15*pow(xq,2) + y) +
                                             pow(xk,2)*(-72 + 41*pow(xq,2) + 17*pow(xq,4) +
                                                28*y - 23*pow(xq,2)*y - 19*pow(y,2)) -
                                             (-1 + y)*(40 - 23*pow(xq,2) - 6*pow(xq,4) + 3*y +
                                                7*pow(xq,2)*y + 7*pow(y,2)) +
                                             B*(-68 - pow(xk,4) + 49*pow(xq,2) + 8*pow(xq,4) +
                                                31*y - 35*pow(xq,2)*y + 15*pow(y,2) +
                                                pow(xk,2)*(90 - 78*pow(xq,2) + 14*y))) +
                                          pow(xk,2)*(-2 + 2*pow(B,4) -
                                             2*pow(B,3)*(2 + 6*pow(xq,2) - 5*y) +
                                             2*pow(xq,6)*(-1 + y) - 12*y + 27*pow(y,2) -
                                             14*pow(y,3) + pow(xq,4)*(-5 + 16*y - 8*pow(y,2)) +
                                             pow(xq,2)*(15 - 27*y + 4*pow(y,2) + 6*pow(y,3)) -
                                             2*B*(-2 + pow(xq,6) - 19*pow(xq,2)*(-1 + y) - 15*y +
                                                19*pow(y,2) - 4*pow(y,3) + pow(xq,4)*(-3 + 5*y))\
                                              + pow(B,2)*(3*pow(xq,4) + pow(xq,2)*(35 - 19*y) +
                                                y*(-28 + 15*y))) -
                                          (-1 + B)*(pow(xq,6)*(1 + B - y) -
                                             pow(xq,4)*(2 + 3*pow(B,2) - 3*y + pow(y,2) +
                                                B*(-5 + 6*y)) -
                                             y*(-3 + 2*pow(B,3) + 8*y - 8*pow(y,2) + 3*pow(y,3) +
                                                pow(B,2)*(-7 + 6*y) + B*(8 - 14*y + 7*pow(y,2))) +
                                             pow(xq,2)*(-3 + 2*pow(B,3) + 10*y - 12*pow(y,2) +
                                                5*pow(y,3) + pow(B,2)*(-7 + 9*y) +
                                                B*(8 - 19*y + 12*pow(y,2)))) +
                                          rl*(2*pow(xk,6)*(5 - 3*pow(xq,2) + y) +
                                             2*pow(B,3)*(2 + 7*pow(xk,2) - 9*pow(xq,2) + 7*y) +
                                             pow(xk,4)*(-37 + 2*pow(xq,4) +
                                                pow(xq,2)*(29 - 7*y) + 13*y + 7*pow(y,2)) +
                                             pow(B,2)*(-16 + 2*pow(xk,4) + 11*pow(xq,4) +
                                                pow(xq,2)*(42 - 28*y) - 22*y + 13*pow(y,2) +
                                                pow(xk,2)*(25 - 64*pow(xq,2) + 15*y)) +
                                             pow(xk,2)*(33 - 29*y + 2*pow(y,2) + 8*pow(y,3) -
                                                2*pow(xq,4)*(7 + 3*y) +
                                                pow(xq,2)*(-26 + 42*y - 6*pow(y,2))) +
                                             (-1 + y)*(8 + pow(xq,6) + pow(xq,4)*(-8 + y) - 2*y +
                                                2*pow(y,2) + 3*pow(y,3) +
                                                pow(xq,2)*(-10 + 14*y - 5*pow(y,2))) -
                                             B*(-20 + 2*pow(xk,6) + pow(xq,6) +
                                                pow(xq,4)*(19 - 12*y) + 2*y + 9*pow(y,2) +
                                                pow(xk,4)*(-31 + 33*pow(xq,2) + 6*y) -
                                                2*pow(xk,2)*
                                                 (-36 + 8*pow(xq,4) + pow(xq,2)*(47 - 19*y) + 5*y -
                                                   2*pow(y,2)) +
                                                pow(xq,2)*(34 - 52*y + 19*pow(y,2))))) +
                                       2*(-2*pow(xk,2)*(-1 + pow(xq,2))*
                                           pow(pow(xq,2)*(-1 + B + y) -
                                             y*(-1 + B + pow(xk,2) + y),2) +
                                          pow(rl,4)*(3*pow(xk,4) +
                                             pow(xk,2)*(21 + 10*B - 34*pow(xq,2) + 13*y) +
                                             10*(-1 + B*y + pow(y,2) - pow(xq,2)*(-2 + B + 2*y))) +
                                          pow(rl,3)*(22*pow(B,2)*(pow(xk,2) - pow(xq,2) + y) +
                                             pow(xk,4)*(16 - 27*pow(xq,2) + 2*y) +
                                             (-1 + y)*(-12 + 33*pow(xq,2) + 6*pow(xq,4) - 17*y -
                                                13*pow(xq,2)*y + 3*pow(y,2)) +
                                             pow(xk,2)*(-40 + 49*pow(xq,2) + 29*pow(xq,4) - 2*y -
                                                41*pow(xq,2)*y + 5*pow(y,2)) +
                                             B*(-12 + 9*pow(xk,4) + 53*pow(xq,2) + 6*pow(xq,4) -
                                                29*y - 43*pow(xq,2)*y + 25*pow(y,2) +
                                                pow(xk,2)*(28 - 84*pow(xq,2) + 34*y))) +
                                          pow(rl,2)*(14*pow(B,3)*(pow(xk,2) - pow(xq,2) + y) +
                                             pow(xk,6)*(4 - 4*pow(xq,2) + 2*y) +
                                             pow(xk,4)*(-23 + 6*pow(xq,4) + 3*y + 7*pow(y,2) -
                                                3*pow(xq,2)*(-9 + 5*y)) +
                                             pow(B,2)*(-2 + 10*pow(xk,4) + 9*pow(xq,4) +
                                                pow(xq,2)*(42 - 32*y) - 38*y + 21*pow(y,2) +
                                                pow(xk,2)*(5 - 64*pow(xq,2) + 31*y)) +
                                             (-1 + y)*(2 + pow(xq,6) + pow(xq,4)*(-8 + y) + 10*y -
                                                8*pow(y,2) + 3*pow(y,3) +
                                                pow(xq,2)*(-16 + 20*y - 5*pow(y,2))) -
                                             pow(xk,2)*(-23 + 8*pow(xq,6) +
                                                pow(xq,4)*(24 - 14*y) + 9*y + 10*pow(y,2) -
                                                8*pow(y,3) + 2*pow(xq,2)*(11 - 23*y + 9*pow(y,2)))
                                               + B*(4 + 2*pow(xk,6) - pow(xq,6) + 32*y -
                                                39*pow(y,2) + 10*pow(y,3) +
                                                pow(xq,4)*(-17 + 12*y) +
                                                pow(xk,4)*(13 - 39*pow(xq,2) + 12*y) +
                                                pow(xq,2)*(-44 + 68*y - 25*pow(y,2)) +
                                                pow(xk,2)*(-42 + 38*pow(xq,4) +
                                                   pow(xq,2)*(86 - 66*y) - 22*y + 20*pow(y,2)))) +
                                          rl*(2*pow(xk,6)*(pow(B,2) + (-2 + pow(xq,2))*y +
                                                B*(1 - 2*pow(xq,2) + y)) +
                                             pow(xk,4)*(4 + 4*pow(B,3) -
                                                2*pow(B,2)*(8*pow(xq,2) - 5*y) + 7*y -
                                                7*pow(y,2) - pow(xq,2)*(8 + y) +
                                                pow(xq,4)*(2 + 4*y) +
                                                B*(-8 + 6*pow(xq,4) + pow(xq,2)*(16 - 11*y) - 9*y +
                                                   7*pow(y,2))) +
                                             pow(xk,2)*(-4 + 2*pow(B,4) -
                                                2*pow(B,3)*(1 + 7*pow(xq,2) - 5*y) -
                                                6*pow(xq,6)*(-1 + y) - 2*y + 11*pow(y,2) -
                                                6*pow(y,3) + pow(xq,4)*(-3 + 8*pow(y,2)) +
                                                pow(B,2)*(-6 + 31*pow(xq,2) + 13*pow(xq,4) -
                                                   18*y - 29*pow(xq,2)*y + 15*pow(y,2)) +
                                                pow(xq,2)*(7 - 13*y + 4*pow(y,2) - 2*pow(y,3)) -
                                                2*B*(-5 + 5*pow(xq,6) + pow(xq,4)*(3 - 7*y) - 5*y +
                                                   11*pow(y,2) - 4*pow(y,3) +
                                                   pow(xq,2)*(12 - 17*y + 8*pow(y,2)))) -
                                             (-1 + B)*(pow(xq,6)*(1 + B - y) -
                                                pow(xq,4)*(2 + 3*pow(B,2) - 3*y + pow(y,2) +
                                                   B*(-5 + 6*y)) -
                                                y*(-3 + 2*pow(B,3) + 8*y - 8*pow(y,2) +
                                                   3*pow(y,3) + pow(B,2)*(-7 + 6*y) +
                                                   B*(8 - 14*y + 7*pow(y,2))) +
                                                pow(xq,2)*(-3 + 2*pow(B,3) + 10*y - 12*pow(y,2) +
                                                   5*pow(y,3) + pow(B,2)*(-7 + 9*y) +
                                                   B*(8 - 19*y + 12*pow(y,2))))))))/
                                   ((-1 + B + rl)*(A + 2*rl)*pow(xk,2)*(-1 + pow(xq,2))*
                                    (pow(xk,2) - pow(xq,2) + y));

 double intfeV = (4*fk*pow(MK,5)*rl*(-(pow(A,3)*
                                         (3 + 2*B - 2*rl - 2*pow(xk,2) + pow(xq,2) - 2*y)*
                                         (-1 + pow(xk,2) + y)) +
                                      pow(A,2)*(-3 + 4*pow(xk,2) - pow(xk,4) + 5*pow(xq,2) -
                                         3*pow(xk,2)*pow(xq,2) + pow(xk,4)*pow(xq,2) -
                                         3*pow(xk,2)*pow(xq,4) +
                                         rl*(11*pow(xk,4) - 3*pow(xk,2)*(9 + pow(xq,2) - 6*y) -
                                            (13 + 6*pow(xq,2) - 7*y)*(-1 + y)) + y +
                                         4*pow(xk,2)*y - 3*pow(xk,4)*y - 10*pow(xq,2)*y +
                                         6*pow(xk,2)*pow(xq,2)*y + 5*pow(y,2) -
                                         6*pow(xk,2)*pow(y,2) + 5*pow(xq,2)*pow(y,2) -
                                         3*pow(y,3) + 10*pow(rl,2)*(-1 + pow(xk,2) + y) +
                                         pow(B,2)*(-2 + 4*pow(xk,2) - 2*pow(xq,2) + 4*y) +
                                         B*(5 + pow(xk,2) + 2*pow(xk,4) -
                                            7*pow(xk,2)*pow(xq,2) - pow(xq,4) +
                                            rl*(4 - 6*pow(xk,2) + 2*pow(xq,2) - 6*y) - 6*y +
                                            4*pow(xk,2)*y + 2*pow(y,2))) -
                                      2*rl*(2*B*pow(xk,6) + 3*pow(xq,2) - 8*B*pow(xq,2) +
                                         7*pow(B,2)*pow(xq,2) - 2*pow(B,3)*pow(xq,2) -
                                         2*pow(xq,4) + B*pow(xq,4) + pow(B,2)*pow(xq,4) -
                                         pow(xq,6) - B*pow(xq,6) - 3*y + 8*B*y - 7*pow(B,2)*y +
                                         2*pow(B,3)*y - 2*pow(xq,2)*y + 7*B*pow(xq,2)*y -
                                         5*pow(B,2)*pow(xq,2)*y + 5*pow(xq,4)*y +
                                         pow(xq,6)*y + 4*pow(y,2) - 8*B*pow(y,2) +
                                         4*pow(B,2)*pow(y,2) - 4*pow(xq,2)*pow(y,2) -
                                         3*pow(xq,4)*pow(y,2) + B*pow(y,3) +
                                         3*pow(xq,2)*pow(y,3) - pow(y,4) +
                                         pow(xk,4)*(-4 + 4*pow(B,2) + 4*pow(xq,2) -
                                            10*B*pow(xq,2) + y + 6*B*y - pow(xq,2)*y - pow(y,2))\
                                          + pow(xk,2)*(4 + 2*pow(B,3) - 2*y + pow(y,2) -
                                            2*pow(y,3) - pow(xq,4)*(1 + 4*y) +
                                            pow(B,2)*(-12*pow(xq,2) + 8*y) +
                                            pow(xq,2)*(-3 + 5*y + 2*pow(y,2)) +
                                            B*(-6 + 7*pow(xq,4) + pow(xq,2)*(15 - 13*y) - 6*y +
                                               5*pow(y,2))) +
                                         pow(rl,2)*(pow(xk,4) +
                                            pow(xk,2)*(9 + 6*B - 16*pow(xq,2) + 7*y) +
                                            6*(pow(-1 + y,2) + B*(-pow(xq,2) + y))) +
                                         rl*(2*pow(xk,6) + pow(xk,4)*(7 - 17*pow(xq,2) + 10*y) +
                                            (-1 + y)*(-2 - 4*pow(xq,4) + pow(xq,2)*(17 - 5*y) -
                                               11*y + 5*pow(y,2)) +
                                            pow(xk,2)*(-15 + 13*pow(xq,4) -
                                               21*pow(xq,2)*(-1 + y) - 11*y + 13*pow(y,2)) +
                                            B*(-2 + 8*pow(xk,4) + 4*pow(xq,4) -
                                               11*pow(xq,2)*(-1 + y) - 7*y + 5*pow(y,2) +
                                               pow(xk,2)*(-11 - 10*pow(xq,2) + 13*y)))) +
                                      A*(-2*B*pow(xk,6) - 3*pow(xq,2) + 8*B*pow(xq,2) -
                                         7*pow(B,2)*pow(xq,2) + 2*pow(B,3)*pow(xq,2) +
                                         2*pow(xq,4) - B*pow(xq,4) - pow(B,2)*pow(xq,4) +
                                         pow(xq,6) + B*pow(xq,6) + 3*y - 8*B*y + 7*pow(B,2)*y -
                                         2*pow(B,3)*y + 2*pow(xq,2)*y - 7*B*pow(xq,2)*y +
                                         5*pow(B,2)*pow(xq,2)*y - 5*pow(xq,4)*y -
                                         pow(xq,6)*y - 4*pow(y,2) + 8*B*pow(y,2) -
                                         4*pow(B,2)*pow(y,2) + 4*pow(xq,2)*pow(y,2) +
                                         3*pow(xq,4)*pow(y,2) - B*pow(y,3) -
                                         3*pow(xq,2)*pow(y,3) + pow(y,4) +
                                         12*pow(rl,3)*(-1 + pow(xk,2) + y) +
                                         pow(xk,4)*(2 - 4*pow(B,2) + B*(2 + 8*pow(xq,2) - 6*y) +
                                            pow(xq,2)*(-2 + y) - y + pow(y,2)) +
                                         pow(rl,2)*(8 + 13*pow(xk,4) +
                                            2*pow(xq,2)*(4 + 5*B - 4*y) - 8*y - 10*B*y +
                                            pow(xk,2)*(-41 - 10*B + 12*pow(xq,2) + 13*y)) -
                                         pow(xk,2)*(2 + 2*pow(B,3) + pow(xq,4)*(1 - 4*y) -
                                            2*pow(B,2)*(1 + 5*pow(xq,2) - 4*y) + pow(y,2) -
                                            2*pow(y,3) + pow(xq,2)*(-3 + 3*y + 2*pow(y,2)) +
                                            B*(-2 + 5*pow(xq,4) + pow(xq,2)*(13 - 11*y) - 8*y +
                                               5*pow(y,2))) +
                                         rl*(-2*pow(xk,6) +
                                            pow(xk,4)*(-7 + 17*pow(xq,2) - 16*y) +
                                            pow(B,2)*(-4 + 8*pow(xk,2) - 4*pow(xq,2) + 8*y) +
                                            (-1 + y)*(8 + 4*pow(xq,4) + 15*y - 11*pow(y,2) +
                                               3*pow(xq,2)*(-9 + 5*y)) -
                                            B*(-12 + 4*pow(xk,4) + 6*pow(xq,4) -
                                               11*pow(xq,2)*(-1 + y) + 5*y + pow(y,2) +
                                               pow(xk,2)*(-17 + 8*pow(xq,2) + 5*y)) +
                                            pow(xk,2)*(19 - 17*pow(xq,4) + 21*y - 25*pow(y,2) +
                                               pow(xq,2)*(-25 + 31*y))))))/
                                  ((A + 2*rl)*pow(xk,2)*(-1 + pow(xq,2))*
                                   (pow(xk,2) - pow(xq,2) + y));
    
    
 double sdh1=  8*pow(MK,6)*pow(xk,4)*(-2*pow(rl,2) + 2*rl*pow(xq,2) +
                                          (pow(xq,2) - y)*(-1 + B + pow(xk,2) + y) +
                                          A*(-1 + 2*B + pow(xk,2) + y));


 double sdh2= (4*pow(MK,6)*rl*pow(xk,4)*pow(1 - pow(xq,2),2)*
               (rl - pow(xq,2))*(pow(A,2) + pow(B,2) + pow(xq,2) +
                 A*(1 - 2*B - pow(xk,2) + pow(xq,2) - 2*y) - y +
                 pow(xk,2)*y - pow(xq,2)*y + pow(y,2) +
                 B*(-1 + pow(xk,2) - pow(xq,2) + 2*y)))/
    pow(-1 + pow(xq,2),2);

 
 double sdfA =2*pow(MK,6)*
    (pow(xk,2)*pow(xq,2) - 2*B*pow(xk,2)*pow(xq,2) +
      2*pow(B,2)*pow(xk,2)*pow(xq,2) - 2*pow(xk,4)*pow(xq,2) +
      2*B*pow(xk,4)*pow(xq,2) + pow(xk,6)*pow(xq,2) +
      2*pow(xq,4) - 4*B*pow(xq,4) + 2*pow(B,2)*pow(xq,4) -
      2*pow(xk,2)*pow(xq,4) + 2*B*pow(xk,2)*pow(xq,4) +
      2*pow(xk,4)*pow(xq,4) + pow(xk,2)*pow(xq,6) -
      2*pow(rl,2)*pow(-1 + pow(xk,2) + pow(xq,2),2) -
      4*pow(xq,2)*y + 8*B*pow(xq,2)*y - 4*pow(B,2)*pow(xq,2)*y +
      4*pow(xk,2)*pow(xq,2)*y - 4*B*pow(xk,2)*pow(xq,2)*y -
      4*pow(xq,4)*y + 4*B*pow(xq,4)*y + 2*pow(y,2) -
      4*B*pow(y,2) + 2*pow(B,2)*pow(y,2) -
      4*pow(xk,2)*pow(y,2) + 4*B*pow(xk,2)*pow(y,2) +
      2*pow(xk,4)*pow(y,2) + 8*pow(xq,2)*pow(y,2) -
      8*B*pow(xq,2)*pow(y,2) -
      4*pow(xk,2)*pow(xq,2)*pow(y,2) + 2*pow(xq,4)*pow(y,2) -
      4*pow(y,3) + 4*B*pow(y,3) + 4*pow(xk,2)*pow(y,3) -
      4*pow(xq,2)*pow(y,3) + 2*pow(y,4) +
      2*pow(A,2)*(pow(xk,4) + pow(-1 + y,2) +
         pow(xk,2)*(-2 - rl + pow(xq,2) + 2*y)) +
      rl*(pow(xk,6) - 2*pow(xk,4)*(1 + B - pow(xq,2)) +
         2*pow(xq,2)*pow(-1 + pow(xq,2),2) +
         pow(xk,2)*(1 - 2*pow(B,2) - 6*pow(xq,2) + 3*pow(xq,4) +
            2*B*(1 + pow(xq,2) - 2*y) + 4*pow(xq,2)*y - 2*pow(y,2))
         ) + 2*A*(pow(xk,2)*pow(xq,4) -
         2*y*(-1 + pow(xk,2) + y)*(-1 + B + pow(xk,2) + y) +
         rl*pow(xk,2)*(-1 + 2*B + pow(xk,2) - pow(xq,2) + 2*y) +
         pow(xq,2)*(pow(xk,4) + 2*B*(-1 + y) + 2*pow(-1 + y,2) +
                      pow(xk,2)*(-3 + 2*y))));
 

 double sdfV = 2*pow(MK,6)*(8*pow(rl,3)*pow(xk,2) + pow(xk,2)*pow(xq,2) -
                              2*B*pow(xk,2)*pow(xq,2) +
                              2*pow(B,2)*pow(xk,2)*pow(xq,2) - 2*pow(xk,4)*pow(xq,2) +
                              2*B*pow(xk,4)*pow(xq,2) + pow(xk,6)*pow(xq,2) -
                              2*pow(xq,4) + 4*B*pow(xq,4) - 2*pow(B,2)*pow(xq,4) +
                              2*pow(xk,2)*pow(xq,4) - 2*B*pow(xk,2)*pow(xq,4) -
                              2*pow(xk,4)*pow(xq,4) + pow(xk,2)*pow(xq,6) +
                              2*pow(rl,2)*(3*pow(xk,4) -
                                 (-1 + pow(xq,2))*(3 + pow(xq,2) - 4*y) -
                                 2*pow(xk,2)*(3 + pow(xq,2) - 2*y)) + 4*pow(xq,2)*y -
                              8*B*pow(xq,2)*y + 4*pow(B,2)*pow(xq,2)*y -
                              8*pow(xk,2)*pow(xq,2)*y + 8*B*pow(xk,2)*pow(xq,2)*y +
                              4*pow(xk,4)*pow(xq,2)*y + 4*pow(xq,4)*y -
                              4*B*pow(xq,4)*y - 4*pow(xk,2)*pow(xq,4)*y - 2*pow(y,2) +
                              4*B*pow(y,2) - 2*pow(B,2)*pow(y,2) +
                              4*pow(xk,2)*pow(y,2) - 4*B*pow(xk,2)*pow(y,2) -
                              2*pow(xk,4)*pow(y,2) - 8*pow(xq,2)*pow(y,2) +
                              8*B*pow(xq,2)*pow(y,2) +
                              8*pow(xk,2)*pow(xq,2)*pow(y,2) - 2*pow(xq,4)*pow(y,2) +
                              4*pow(y,3) - 4*B*pow(y,3) - 4*pow(xk,2)*pow(y,3) +
                              4*pow(xq,2)*pow(y,3) - 2*pow(y,4) -
                              2*pow(A,2)*(pow(xk,4) + pow(-1 + y,2) +
                                 pow(xk,2)*(-2 + rl - pow(xq,2) + 2*y)) +
                              rl*(pow(xk,6) + 2*pow(xq,2)*pow(1 + pow(xq,2) - 2*y,2) +
                                 2*pow(xk,4)*(-1 + 3*B + pow(xq,2) + 2*y) +
                                 pow(xk,2)*(1 + 6*pow(B,2) - 2*pow(xq,2) - 5*pow(xq,4) -
                                    2*B*(3 + pow(xq,2) - 4*y) - 4*y + 8*pow(xq,2)*y +
                                    2*pow(y,2))) - 2*A*
                               (-(pow(xk,2)*pow(xq,4)) +
                                 rl*pow(xk,2)*(-1 + 2*B + pow(xk,2) + pow(xq,2)) -
                                 2*y*(-1 + pow(xk,2) + y)*(-1 + B + pow(xk,2) + y) +
                                 pow(xq,2)*(pow(xk,4) + 2*B*(-1 + y) + 2*pow(-1 + y,2) +
                                              pow(xk,2)*(-3 + 4*y))));
 

 double sdhe1 = 8*pow(MK,6)*pow(A + 2*rl,2)*
    (2*pow(rl,2) - rl*pow(xk,2) - pow(xq,2) - 2*rl*pow(xq,2) +
      2*pow(xk,2)*pow(xq,2) + B*(2*rl + pow(xq,2) - y) + y -
      pow(xk,2)*y + pow(xq,2)*y - pow(y,2) +
     A*(-1 + B + pow(xk,2) + y));

 double sdhe2 = (-4*B*pow(MK,6)*pow(1 - B - rl,2)*rl*pow(A + 2*rl,2)*
                 (y*(-1 + B + y) + rl*(1 + B + y) - A*(-1 + rl + y)))/ pow(-1 + B + rl,2);

 double sdfeA = 2*pow(MK,6)*(pow(A,3)*B +
                               2*pow(A,2)*(pow(B,2) +
                                  B*(-1 + 3*rl + pow(xk,2) + pow(xq,2)) +
                                  pow(-1 + pow(xk,2) + y,2)) +
                               2*(2*pow(B,3)*rl + 6*pow(rl,4) -
                                  4*pow(rl,3)*(1 + 3*pow(xq,2) - y) -
                                  2*rl*(-2*pow(xq,2)*(pow(xq,2) - y)*(-1 + y) +
                                     pow(xk,4)*y + pow(xk,2)*(-1 + y)*y +
                                     pow(xk,2)*pow(xq,2)*(1 + y)) +
                                  pow(pow(xq,2) - pow(xq,2)*y + y*(-1 + pow(xk,2) + y),
                                   2) + pow(B,2)*(12*pow(rl,2) + pow(pow(xq,2) - y,2) +
                                     2*rl*(-2 + pow(xk,2) - pow(xq,2) + 2*y)) +
                                  pow(rl,2)*(2 + 3*pow(xk,4) + 6*pow(xq,4) -
                                     8*pow(xq,2)*(-1 + y) - 8*y + 6*pow(y,2) +
                                     pow(xk,2)*(-4 + 8*y)) +
                                  2*B*(9*pow(rl,3) +
                                     pow(rl,2)*(-5 + pow(xk,2) - 3*pow(xq,2) + 4*y) +
                                     (pow(xq,2) - y)*
                                      (pow(xq,2)*(-1 + y) - y*(-1 + pow(xk,2) + y)) +
                                     rl*(1 + pow(xk,4) + pow(xq,2) + 3*pow(xq,4) - 2*y -
                                        4*pow(xq,2)*y + 2*pow(y,2) +
                                        pow(xk,2)*(-1 - pow(xq,2) + y)))) +
                               A*(pow(B,3) + 2*pow(B,2)*
                                   (-1 + 5*rl + pow(xk,2) + pow(xq,2)) +
                                  B*(1 + 17*pow(rl,2) - 2*pow(xk,2) + 2*pow(xk,4) -
                                     6*pow(xq,2) + 2*pow(xq,4) +
                                     2*rl*(-2 + pow(xk,2) + 5*pow(xq,2) - 2*y) + 4*y -
                                     4*pow(xk,2)*y + 4*pow(xq,2)*y - 4*pow(y,2)) +
                                  2*(pow(rl,3) - 2*pow(rl,2)*
                                      (-1 + pow(xk,2) + pow(xq,2) + y) -
                                     2*(-1 + pow(xk,2) + y)*
                                      (pow(xq,2) - pow(xq,2)*y + y*(-1 + pow(xk,2) + y)) +
                                     rl*(3*pow(xk,4) + pow(-1 + pow(xq,2) + y,2) +
                                         2*pow(xk,2)*(-2 + pow(xq,2) + 2*y)))));
    
    
 double sdfeV = 2*pow(MK,6)*(pow(A,3)*B -
                               2*pow(A,2)*(pow(B,2) -
                                  B*(1 + rl - pow(xk,2) + pow(xq,2) - 2*y) +
                                  pow(-1 + pow(xk,2) + y,2)) +
                               A*(pow(B,3) - 2*pow(B,2)*
                                   (1 + 3*rl - pow(xk,2) + pow(xq,2) - 2*y) +
                                  B*(1 + pow(rl,2) - 2*pow(xk,2) + 2*pow(xk,4) +
                                     6*pow(xq,2) + 2*pow(xq,4) - 8*y + 8*pow(xk,2)*y -
                                     8*pow(xq,2)*y + 8*pow(y,2) -
                                     2*rl*(2 + 5*pow(xk,2) - pow(xq,2) + 2*y)) +
                                  2*(pow(rl,3) + rl*(-pow(xk,4) - 2*pow(xk,2)*pow(xq,2) +
                                        pow(1 + pow(xq,2) - y,2)) +
                                     2*pow(rl,2)*(-1 + pow(xk,2) - pow(xq,2) + y) +
                                     2*(-1 + pow(xk,2) + y)*
                                      (pow(xq,2) - pow(xq,2)*y + y*(-1 + pow(xk,2) + y))))\
                                + 2*(2*pow(B,3)*rl + 2*pow(rl,4) +
                                  4*pow(rl,3)*(-1 + 2*pow(xk,2) - pow(xq,2) + y) +
                                  2*rl*(-2*pow(xq,2)*(pow(xq,2) - y)*(-1 + y) +
                                     pow(xk,4)*y + pow(xk,2)*(-1 + y)*y +
                                     pow(xk,2)*pow(xq,2)*(1 + y)) -
                                  pow(pow(xq,2) - pow(xq,2)*y + y*(-1 + pow(xk,2) + y),
                                   2) + pow(B,2)*(-pow(pow(xq,2) - y,2) +
                                     rl*(-4 + 6*pow(xk,2) - 2*pow(xq,2) + 4*y)) +
                                  pow(rl,2)*(5*pow(xk,4) -
                                     4*pow(xk,2)*(3 + 2*pow(xq,2) - 2*y) +
                                     2*(3 + pow(xq,4) - 4*y + pow(y,2))) +
                                  2*B*(pow(rl,3) + pow(rl,2)*
                                      (-5 + pow(xk,2) - pow(xq,2) + 2*y) -
                                     (pow(xq,2) - y)*
                                      (pow(xq,2)*(-1 + y) - y*(-1 + pow(xk,2) + y)) +
                                     rl*(1 + 3*pow(xk,4) + pow(xq,2) + pow(xq,4) - 2*y -
                                        2*pow(xq,2)*y + 2*pow(y,2) +
                                         pow(xk,2)*(-3 - 3*pow(xq,2) + 5*y)))));


 double sdh1h2 = (-8*pow(MK,6)*rl*pow(xk,4)*(1 - pow(xq,2))*
                  (-2*pow(B,2) + rl*pow(xk,2) - pow(xq,2) +
                    B*(2 - 2*pow(xk,2) + pow(xq,2) - 3*y) + y - pow(xk,2)*y +
                    pow(xq,2)*y - pow(y,2) + A*(-1 + 2*B + pow(xk,2) + y)))/ (-1 + pow(xq,2));

 double sdh1fA = -4*pow(MK,6)*pow(xk,2)*(pow(xq,2) - 3*B*pow(xq,2) +
                                             2*pow(B,2)*pow(xq,2) - 3*pow(xk,2)*pow(xq,2) +
                                             3*B*pow(xk,2)*pow(xq,2) + 2*pow(xk,4)*pow(xq,2) -
                                             pow(xq,4) + B*pow(xq,4) + 2*pow(xk,2)*pow(xq,4) -
                                             4*pow(rl,2)*(-1 + pow(xk,2) + pow(xq,2)) - y + 3*B*y -
                                             2*pow(B,2)*y + 2*pow(xk,2)*y - 3*B*pow(xk,2)*y -
                                             pow(xk,4)*y + B*pow(xq,2)*y + pow(xq,4)*y + pow(y,2) -
                                             2*B*pow(y,2) - pow(xk,2)*pow(y,2) - pow(xq,2)*pow(y,2) +
                                             2*pow(A,2)*(-1 + pow(xk,2) + y) +
                                             A*(2*B*(-1 + pow(xk,2) + pow(xq,2)) +
                                                (-1 + pow(xk,2) + 3*pow(xq,2) - 2*y)*(-1 + pow(xk,2) + y))
                                               + rl*(pow(xk,4) + 4*pow(xq,2)*(-1 + pow(xq,2)) +
                                                     pow(xk,2)*(-1 + 3*pow(xq,2) + 2*y)));

 double sdh1fV =4*pow(MK,6)*pow(xk,2)*(-pow(xq,2) + 3*B*pow(xq,2) -
                                           2*pow(B,2)*pow(xq,2) + 3*pow(xk,2)*pow(xq,2) -
                                           3*B*pow(xk,2)*pow(xq,2) - 2*pow(xk,4)*pow(xq,2) -
                                           pow(xq,4) + B*pow(xq,4) + 2*pow(xk,2)*pow(xq,4) -
                                           4*pow(rl,2)*(-1 + pow(xk,2) + pow(xq,2)) +
                                           rl*(pow(xk,2) - pow(xk,4) - 5*pow(xk,2)*pow(xq,2) +
                                              4*pow(xq,2)*(1 + pow(xq,2) - 2*y)) + y - 3*B*y +
                                           2*pow(B,2)*y - 2*pow(xk,2)*y + 3*B*pow(xk,2)*y +
                                           pow(xk,4)*y + 4*pow(xq,2)*y - 5*B*pow(xq,2)*y -
                                           6*pow(xk,2)*pow(xq,2)*y + pow(xq,4)*y - 3*pow(y,2) +
                                           4*B*pow(y,2) + 3*pow(xk,2)*pow(y,2) -
                                           3*pow(xq,2)*pow(y,2) + 2*pow(y,3) +
                                           2*pow(A,2)*(-1 + pow(xk,2) + y) -
                                           A*(2*B*(-1 + pow(xk,2) - pow(xq,2) + 2*y) +
                                              (-1 + pow(xk,2) + y)*(-1 + pow(xk,2) - 3*pow(xq,2) + 4*y)));
    
    
 double sdh1he1 =16*pow(MK,6)*(A + 2*rl)*pow(xk,2)*
    (-pow(rl,2) + A*(-1 + B + pow(xk,2) + y) +
      (pow(xq,2) - y)*(-1 + B + pow(xk,2) + y) +
     rl*(-1 + 2*B + pow(xk,2) + pow(xq,2) + y));
    

 double sdh1he2 =(-4*pow(MK,6)*(1 - B - rl)*rl*(A + 2*rl)*pow(xk,2)*
                  (2*pow(B,2) - 2*rl + 2*pow(rl,2) + rl*pow(xk,2) +
                    pow(xq,2) - 2*rl*pow(xq,2) + y + 2*rl*y - pow(xk,2)*y -
                    pow(xq,2)*y - pow(y,2) + A*(-1 - 2*B + pow(xk,2) + y) +
                   B*(-2 + 2*rl - pow(xq,2) + 3*y)))/(-1 + B + rl);
    
 double sdh1feA = -4*pow(MK,6)*pow(xk,2)*(-2*pow(rl,2) - 8*pow(rl,3) -
                                              rl*pow(xk,2) + 6*pow(rl,2)*pow(xk,2) + rl*pow(xk,4) +
                                              pow(xq,2) - 6*rl*pow(xq,2) + 6*pow(rl,2)*pow(xq,2) -
                                              pow(xk,2)*pow(xq,2) + 3*rl*pow(xk,2)*pow(xq,2) -
                                              pow(xq,4) + 2*rl*pow(xq,4) +
                                              pow(B,2)*(4*rl + pow(xq,2) - y) - y + 4*rl*y +
                                              4*pow(rl,2)*y + 2*pow(xk,2)*y - 4*rl*pow(xk,2)*y -
                                              pow(xk,4)*y + 4*rl*pow(xq,2)*y + pow(xq,4)*y + pow(y,2) -
                                              4*rl*pow(y,2) - pow(xk,2)*pow(y,2) -
                                              pow(xq,2)*pow(y,2) + pow(A,2)*(-1 + 2*B + pow(xk,2) + y) +
                                              B*(8*pow(rl,2) + (-2 + pow(xk,2))*pow(xq,2) + pow(xq,4) +
                                                 rl*(-4 + 3*pow(xk,2) + 8*pow(xq,2) - 4*y) -
                                                 y*(-2 + 2*pow(xk,2) + y)) +
                                              A*(1 + 2*pow(B,2) + 10*B*rl - 2*pow(rl,2) - 2*pow(xk,2) +
                                                 pow(xk,4) - 2*pow(xq,2) + pow(xk,2)*pow(xq,2) +
                                                 3*B*(-1 + pow(xk,2) + pow(xq,2)) + 2*pow(xq,2)*y -
                                                 pow(y,2) + rl*(-6 + 7*pow(xk,2) + 2*pow(xq,2) + 6*y)));

 double sdh1feV = 4*pow(MK,6)*pow(xk,2)*(6*pow(rl,2) - 4*pow(rl,3) +
                                             rl*pow(xk,2) - 6*pow(rl,2)*pow(xk,2) - rl*pow(xk,4) -
                                             pow(xq,2) + 2*rl*pow(xq,2) + 2*pow(rl,2)*pow(xq,2) +
                                             pow(xk,2)*pow(xq,2) - rl*pow(xk,2)*pow(xq,2) -
                                             pow(xq,4) + 2*rl*pow(xq,4) + y - 4*pow(rl,2)*y -
                                             2*pow(xk,2)*y - 2*rl*pow(xk,2)*y + pow(xk,4)*y +
                                             4*pow(xq,2)*y - 4*rl*pow(xq,2)*y -
                                             2*pow(xk,2)*pow(xq,2)*y + pow(xq,4)*y - 3*pow(y,2) +
                                             3*pow(xk,2)*pow(y,2) - 3*pow(xq,2)*pow(y,2) +
                                             2*pow(y,3) + pow(A,2)*(-1 + 2*B + pow(xk,2) + y) +
                                             pow(B,2)*(-4*rl - pow(xq,2) + y) +
                                             B*(-4*pow(rl,2) + pow(xq,4) +
                                                rl*(4 - 7*pow(xk,2) + 4*pow(xq,2) - 8*y) +
                                                y*(-2 + 2*pow(xk,2) + 3*y) -
                                                pow(xq,2)*(-2 + pow(xk,2) + 4*y)) -
                                             A*(1 + 2*pow(B,2) + 2*pow(rl,2) - 2*pow(xk,2) +
                                                pow(xk,4) + 2*pow(xq,2) - pow(xk,2)*pow(xq,2) - 4*y +
                                                4*pow(xk,2)*y - 2*pow(xq,2)*y + 3*pow(y,2) +
                                                rl*(-2 + pow(xk,2) - 2*pow(xq,2) + 2*y) +
                                                B*(-3 - 2*rl + 3*pow(xk,2) - 3*pow(xq,2) + 6*y)));

 double sdh2fA= (4*pow(MK,6)*rl*pow(xk,2)*(1 - pow(xq,2))*
                 (-pow(xq,2) + 3*B*pow(xq,2) - 2*pow(B,2)*pow(xq,2) +
                   pow(xk,2)*pow(xq,2) - 3*B*pow(xk,2)*pow(xq,2) -
                   pow(xq,4) + B*pow(xq,4) +
                   rl*pow(xk,2)*(-1 + pow(xk,2) + pow(xq,2)) + y - 3*B*y +
                   2*pow(B,2)*y - 2*pow(xk,2)*y + 3*B*pow(xk,2)*y +
                   pow(xk,4)*y + 4*pow(xq,2)*y - 5*B*pow(xq,2)*y -
                   2*pow(xk,2)*pow(xq,2)*y + pow(xq,4)*y - 3*pow(y,2) +
                   4*B*pow(y,2) + 3*pow(xk,2)*pow(y,2) -
                   3*pow(xq,2)*pow(y,2) + 2*pow(y,3) +
                   2*pow(A,2)*(-1 + pow(xk,2) + y) -
                   A*(2*B*(-1 + pow(xk,2) - pow(xq,2) + 2*y) +
                      (-1 + pow(xk,2) + y)*
                       (-1 + pow(xk,2) - 3*pow(xq,2) + 4*y))))/ (-1 + pow(xq,2));

 double sdh2he1 = (4*pow(MK,6)*rl*(A + 2*rl)*pow(xk,2)*(1 - pow(xq,2))*
                   (2*pow(B,2) - 2*rl + rl*pow(xk,2) + 3*pow(xq,2) -
                     2*rl*pow(xq,2) - 2*pow(xk,2)*pow(xq,2) + 2*pow(xq,4) -
                     y + 4*rl*y + pow(xk,2)*y - 5*pow(xq,2)*y + pow(y,2) -
                     A*(-1 + 2*B + 4*rl + pow(xk,2) - 4*pow(xq,2) + y) +
                     B*(-2 + 4*rl + 2*pow(xk,2) - 5*pow(xq,2) + 3*y)))/ (-1 + pow(xq,2));

 double sdh2he2 = (2*pow(MK,6)*(1 - B - rl)*rl*(A + 2*rl)*pow(xk,2)*(1 - pow(xq,2))*
                   (2*pow(rl,2) - 2*rl*pow(xq,2) + 2*rl*y - 2*pow(rl,2)*y -
                     3*rl*pow(xk,2)*y + pow(xq,2)*y + 2*rl*pow(xq,2)*y -
                     pow(y,2) - 2*rl*pow(y,2) + pow(xk,2)*pow(y,2) -
                     pow(xq,2)*pow(y,2) + pow(y,3) +
                     pow(A,2)*(-1 + pow(xk,2) + y) +
                     pow(B,2)*(pow(xq,2) + y) -
                     B*(2*pow(rl,2) + pow(xq,2) - y*(-1 + pow(xk,2) + 2*y) +
                        rl*(-2 + 3*pow(xk,2) - 2*pow(xq,2) + 2*y)) -
                     A*(-2*pow(rl,2) + pow(xq,2) +
                        rl*(2 - 3*pow(xk,2) + 2*pow(xq,2) - 2*y) - 2*y +
                        2*pow(xk,2)*y - pow(xq,2)*y + 2*pow(y,2) +
                        B*(-1 + pow(xk,2) + pow(xq,2) + 2*y))))/ ((-1 + B + rl)*(-1 + pow(xq,2)));

 double sdh2feA =(2*pow(MK,6)*rl*pow(xk,2)*(1 - pow(xq,2))*
                  (-2*pow(B,3) - 2*rl + 4*pow(rl,2) + 2*rl*pow(xk,2) -
                    3*pow(rl,2)*pow(xk,2) + 2*pow(xq,2) - 9*rl*pow(xq,2) +
                    6*pow(rl,2)*pow(xq,2) - 2*pow(xk,2)*pow(xq,2) +
                    6*rl*pow(xk,2)*pow(xq,2) + 2*pow(xq,4) -
                    6*rl*pow(xq,4) + 2*pow(A,2)*
                     (1 + B + rl - pow(xk,2) - pow(xq,2) - y) - 2*y + 5*rl*y -
                    10*pow(rl,2)*y + 4*pow(xk,2)*y - 2*rl*pow(xk,2)*y -
                    2*pow(xk,4)*y - 7*pow(xq,2)*y + 15*rl*pow(xq,2)*y +
                    4*pow(xk,2)*pow(xq,2)*y - 2*pow(xq,4)*y + 5*pow(y,2) -
                    3*rl*pow(y,2) - 5*pow(xk,2)*pow(y,2) +
                    5*pow(xq,2)*pow(y,2) - 3*pow(y,3) -
                    2*pow(B,2)*(-2 + 6*rl + 2*pow(xk,2) - 3*pow(xq,2) + 3*y) +
                    A*(2 + 10*pow(rl,2) - 4*pow(xk,2) + 2*pow(xk,4) +
                       4*pow(xq,2) - 2*pow(xk,2)*pow(xq,2) - 7*y +
                       7*pow(xk,2)*y - 2*pow(xq,2)*y + 5*pow(y,2) +
                       rl*(-3 + pow(xk,2) - 10*pow(xq,2) + y) +
                       2*B*(-2 + 5*rl + 2*pow(xk,2) - pow(xq,2) + 2*y)) -
                    B*(2 + 10*pow(rl,2) + 2*pow(xk,4) + 8*pow(xq,2) +
                       4*pow(xq,4) - 8*y - 11*pow(xq,2)*y + 7*pow(y,2) +
                       pow(xk,2)*(-4 - 6*pow(xq,2) + 8*y) +
                       rl*(-10 + 6*pow(xk,2) - 17*pow(xq,2) + 15*y))))/ (-1 + pow(xq,2));
    
    
 double sdh2feV = (2*pow(MK,6)*rl*pow(xk,2)*(1 - pow(xq,2))*
                   (-2*pow(B,3) - 2*rl + 4*pow(rl,2) + 4*rl*pow(xk,2) -
                     pow(rl,2)*pow(xk,2) - 2*rl*pow(xk,4) - 7*rl*pow(xq,2) +
                     2*pow(rl,2)*pow(xq,2) + 4*rl*pow(xk,2)*pow(xq,2) -
                     2*rl*pow(xq,4) - 2*pow(A,2)*
                      (1 + B - rl - pow(xk,2) + pow(xq,2) - y) + 7*rl*y -
                     6*pow(rl,2)*y - 8*rl*pow(xk,2)*y + pow(xq,2)*y +
                     9*rl*pow(xq,2)*y - pow(y,2) - 5*rl*pow(y,2) +
                     pow(xk,2)*pow(y,2) - pow(xq,2)*pow(y,2) + pow(y,3) -
                     4*pow(B,2)*(-1 + pow(xk,2) - pow(xq,2) + y) +
                     A*(4*pow(B,2) + 6*pow(rl,2) - 2*pow(xq,2) -
                        2*B*(1 + rl - pow(xk,2) + pow(xq,2) - y) + 3*y -
                        3*pow(xk,2)*y + 4*pow(xq,2)*y - 3*pow(y,2) +
                        rl*(-5 + 7*pow(xk,2) - 6*pow(xq,2) + 3*y)) -
                     B*(2 + 6*pow(rl,2) + 2*pow(xk,4) + 4*pow(xq,2) +
                        2*pow(xq,4) - 4*pow(xk,2)*(1 + pow(xq,2) - y) - 4*y -
                        3*pow(xq,2)*y + pow(y,2) +
                        rl*(-6 + 8*pow(xk,2) - 7*pow(xq,2) + 5*y))))/ (-1 + pow(xq,2));

 double sdfAfV = 4*pow(MK,6)*(2*pow(rl,2)*pow(-1 + pow(xk,2) + pow(xq,2),2) +
                                pow(xk,2)*pow(xq,2)*
                                 (-2*pow(A,2) + 2*pow(B,2) - 2*A*(pow(xq,2) - y) +
                                   2*B*(-1 + pow(xk,2) + y) +
                                   (-1 + pow(xk,2) + pow(xq,2))*
                                    (-1 + pow(xk,2) - pow(xq,2) + 2*y)) +
                                rl*(pow(xk,6) - 2*pow(xq,2)*(-1 + pow(xq,2))*
                                    (1 + pow(xq,2) - 2*y) +
                                   2*pow(xk,4)*(-1 - A + B + pow(xq,2) + y) +
                                   pow(xk,2)*(1 + 2*pow(A,2) + 2*pow(B,2) - 4*pow(xq,2) +
                                      pow(xq,4) + A*(2 - 4*B + 2*pow(xq,2) - 4*y) -
                                      2*B*(1 + pow(xq,2) - 2*y) - 2*y + 2*pow(xq,2)*y +
                                                2*pow(y,2))));
    
    
 double sdfAhe1 = -4*pow(MK,6)*(A + 2*rl)*(pow(xq,2) - 2*B*pow(xq,2) +
                                             pow(B,2)*pow(xq,2) - 3*pow(xk,2)*pow(xq,2) +
                                             3*B*pow(xk,2)*pow(xq,2) + 2*pow(xk,4)*pow(xq,2) -
                                             pow(xq,4) + B*pow(xq,4) + 2*pow(xk,2)*pow(xq,4) -
                                             2*pow(rl,2)*(-1 + pow(xk,2) + pow(xq,2)) - y + 2*B*y -
                                             pow(B,2)*y + 2*pow(xk,2)*y - 2*B*pow(xk,2)*y -
                                             pow(xk,4)*y + pow(xq,4)*y + pow(y,2) - B*pow(y,2) -
                                             pow(xk,2)*pow(y,2) - pow(xq,2)*pow(y,2) +
                                             pow(A,2)*(-1 + pow(xk,2) + y) +
                                             A*(1 - 2*pow(xk,2) + pow(xk,4) - 2*pow(xq,2) +
                                                3*pow(xk,2)*pow(xq,2) +
                                                B*(-1 + pow(xk,2) + pow(xq,2)) + 2*pow(xq,2)*y -
                                                pow(y,2) + rl*(-4 + 3*pow(xk,2) + 4*y)) +
                                             rl*(pow(xk,4) + pow(xk,2)*(-1 + B + 5*pow(xq,2) - 2*y) +
                                                2*(pow(xq,4) - 2*y*(-1 + B + y) +
                                                   pow(xq,2)*(-3 + 2*B + 2*y))));
    
    
 double sdfAhe2 =(-2*pow(MK,6)*(1 - B - rl)*rl*(A + 2*rl)*
                  (-2*pow(B,2)*pow(xk,2) + 3*B*pow(xq,2) -
                    pow(B,2)*pow(xq,2) - pow(xk,2)*pow(xq,2) +
                    B*pow(xk,2)*pow(xq,2) - pow(xq,4) + B*pow(xq,4) -
                    2*pow(rl,2)*(-1 + pow(xk,2) + pow(xq,2)) + 2*y - 3*B*y +
                    pow(B,2)*y - 3*pow(xk,2)*y + pow(xk,4)*y +
                    2*pow(xq,2)*y - 5*B*pow(xq,2)*y +
                    2*pow(xk,2)*pow(xq,2)*y + pow(xq,4)*y - 5*pow(y,2) +
                    4*B*pow(y,2) + 4*pow(xk,2)*pow(y,2) -
                    2*pow(xq,2)*pow(y,2) + 3*pow(y,3) +
                    pow(A,2)*(-1 + pow(xk,2) + y) +
                    A*(-2 + 3*pow(xk,2) - pow(xk,4) - 2*pow(xq,2) -
                       pow(xk,2)*pow(xq,2) + rl*(2 + pow(xk,2) - 2*y) +
                       B*(1 + pow(xk,2) + pow(xq,2) - 2*y) + 6*y -
                       5*pow(xk,2)*y + 2*pow(xq,2)*y - 4*pow(y,2)) -
                    rl*(pow(xk,4) + pow(xk,2)*(-2 + 3*B - pow(xq,2) + 3*y) -
                       2*(pow(xq,4) - pow(xq,2)*y + (-1 + y)*y +
                          B*(-pow(xq,2) + y)))))/(-1 + B + rl);
    
 double sdfAfeA = 2*pow(MK,6)*(-pow(xq,2) + 3*B*pow(xq,2) -
                                 3*pow(B,2)*pow(xq,2) + pow(B,3)*pow(xq,2) +
                                 2*pow(xk,2)*pow(xq,2) - 4*B*pow(xk,2)*pow(xq,2) +
                                 2*pow(B,2)*pow(xk,2)*pow(xq,2) - pow(xk,4)*pow(xq,2) +
                                 B*pow(xk,4)*pow(xq,2) - 2*pow(xk,2)*pow(xq,4) +
                                 2*B*pow(xk,2)*pow(xq,4) - pow(xq,6) + B*pow(xq,6) -
                                 8*pow(rl,3)*(-1 + pow(xk,2) + pow(xq,2)) + y - 3*B*y +
                                 3*pow(B,2)*y - pow(B,3)*y - 3*pow(xk,2)*y +
                                 6*B*pow(xk,2)*y - 3*pow(B,2)*pow(xk,2)*y + 3*pow(xk,4)*y -
                                 3*B*pow(xk,4)*y - pow(xk,6)*y + 3*pow(xq,2)*y -
                                 6*B*pow(xq,2)*y + 3*pow(B,2)*pow(xq,2)*y -
                                 2*pow(xk,2)*pow(xq,2)*y + 2*B*pow(xk,2)*pow(xq,2)*y -
                                 pow(xk,4)*pow(xq,2)*y + 3*pow(xq,4)*y - 3*B*pow(xq,4)*y +
                                 pow(xk,2)*pow(xq,4)*y + pow(xq,6)*y - 3*pow(y,2) +
                                 6*B*pow(y,2) - 3*pow(B,2)*pow(y,2) +
                                 6*pow(xk,2)*pow(y,2) - 6*B*pow(xk,2)*pow(y,2) -
                                 3*pow(xk,4)*pow(y,2) - 6*pow(xq,2)*pow(y,2) +
                                 6*B*pow(xq,2)*pow(y,2) +
                                 2*pow(xk,2)*pow(xq,2)*pow(y,2) - 3*pow(xq,4)*pow(y,2) +
                                 4*pow(y,3) - 4*B*pow(y,3) - 4*pow(xk,2)*pow(y,3) +
                                 4*pow(xq,2)*pow(y,3) - 2*pow(y,4) +
                                 pow(A,3)*(-1 + pow(xk,2) + y) +
                                 pow(A,2)*(-3*pow(xq,2) + 2*pow(xk,2)*pow(xq,2) + 3*y -
                                    3*pow(xk,2)*y + 3*pow(xq,2)*y - 3*pow(y,2) +
                                    B*(-2 + 2*pow(xk,2) + pow(xq,2) + y) +
                                    rl*(-6 + 7*pow(xk,2) + 6*y)) +
                                 pow(rl,2)*(3*pow(xk,4) +
                                    pow(xk,2)*(-2 + 15*pow(xq,2) - 9*y) -
                                    2*B*(-1 + pow(xk,2) - 4*pow(xq,2) + 5*y) +
                                    2*(-1 - 7*pow(xq,2) + 3*pow(xq,4) + 5*y + 5*pow(xq,2)*y -
                                       5*pow(y,2))) + rl*
                                  (pow(xk,6) + 2*pow(xq,6) +
                                    pow(xk,4)*(-2 + 2*B + 4*pow(xq,2) - y) +
                                    pow(xq,4)*(-5 + 5*B + y) +
                                    pow(xk,2)*(1 + pow(B,2) - 7*pow(xq,2) + 5*pow(xq,4) +
                                       B*(-2 + 15*pow(xq,2) - 7*y) + 5*y) +
                                    2*pow(xq,2)*(3 - 7*B + 4*pow(B,2) - y + 2*B*y -
                                       pow(y,2)) + y*(-4 + 12*B - 8*pow(B,2) + 3*y - 7*B*y +
                                       pow(y,2))) + A*(-1 + 3*pow(xk,2) - 3*pow(xk,4) +
                                    pow(xk,6) - 2*pow(xk,2)*pow(xq,2) +
                                    2*pow(xk,4)*pow(xq,2) - 3*pow(xq,4) +
                                    pow(xk,2)*pow(xq,4) +
                                    pow(B,2)*(-1 + pow(xk,2) + 2*pow(xq,2) - y) + 3*y -
                                    6*pow(xk,2)*y + 3*pow(xk,4)*y + 6*pow(xq,2)*y -
                                    2*pow(xk,2)*pow(xq,2)*y + 3*pow(xq,4)*y - 6*pow(y,2) +
                                    6*pow(xk,2)*pow(y,2) - 6*pow(xq,2)*pow(y,2) +
                                    4*pow(y,3) + 2*pow(rl,2)*
                                     (-4 + 4*pow(xk,2) - pow(xq,2) + 5*y) +
                                    2*B*(1 + pow(xk,4) - pow(xq,2) + pow(xq,4) - y -
                                       pow(xq,2)*y + pow(y,2) +
                                       rl*(-4 + 3*pow(xk,2) + 3*pow(xq,2) + y) +
                                       pow(xk,2)*(-2 + 4*pow(xq,2) + y)) +
                                    rl*(4 + 3*pow(xk,4) + 2*pow(xq,4) +
                                       pow(xk,2)*(-7 + 13*pow(xq,2) - 6*y) + 3*y -
                                        7*pow(y,2) + pow(xq,2)*(-11 + 9*y))));

 double sdfAfeV = -2*pow(MK,6)*(pow(xq,2) - 3*B*pow(xq,2) +
                                  3*pow(B,2)*pow(xq,2) - pow(B,3)*pow(xq,2) -
                                  2*pow(xk,2)*pow(xq,2) + 4*B*pow(xk,2)*pow(xq,2) -
                                  2*pow(B,2)*pow(xk,2)*pow(xq,2) + pow(xk,4)*pow(xq,2) -
                                  B*pow(xk,4)*pow(xq,2) - pow(xq,6) + B*pow(xq,6) -
                                  4*pow(rl,3)*(-1 + pow(xk,2) + pow(xq,2)) - y + 3*B*y -
                                  3*pow(B,2)*y + pow(B,3)*y + 3*pow(xk,2)*y -
                                  6*B*pow(xk,2)*y + 3*pow(B,2)*pow(xk,2)*y - 3*pow(xk,4)*y +
                                  3*B*pow(xk,4)*y + pow(xk,6)*y - 3*pow(xq,2)*y +
                                  6*B*pow(xq,2)*y - 3*pow(B,2)*pow(xq,2)*y +
                                  4*pow(xk,2)*pow(xq,2)*y - 4*B*pow(xk,2)*pow(xq,2)*y -
                                  pow(xk,4)*pow(xq,2)*y + 3*pow(xq,4)*y - 3*B*pow(xq,4)*y -
                                  pow(xk,2)*pow(xq,4)*y + pow(xq,6)*y + 3*pow(y,2) -
                                  6*B*pow(y,2) + 3*pow(B,2)*pow(y,2) -
                                  6*pow(xk,2)*pow(y,2) + 6*B*pow(xk,2)*pow(y,2) +
                                  3*pow(xk,4)*pow(y,2) - 3*pow(xq,4)*pow(y,2) -
                                  2*pow(y,3) + 2*B*pow(y,3) + 2*pow(xk,2)*pow(y,3) +
                                  2*pow(xq,2)*pow(y,3) + pow(A,3)*(-1 + pow(xk,2) + y) -
                                  pow(rl,2)*(2 + 7*pow(xk,4) - 6*pow(xq,2) - 2*pow(xq,4) +
                                     B*(2 + 4*pow(xk,2) + 4*pow(xq,2) - 6*y) + 6*y +
                                     6*pow(xq,2)*y - 6*pow(y,2) +
                                     pow(xk,2)*(-10 + 3*pow(xq,2) + 3*y)) -
                                  rl*(pow(xk,6) - 2*pow(xq,6) +
                                     pow(xk,4)*(-2 + 6*B + 2*pow(xq,2) - y) +
                                     pow(xq,4)*(3 - 3*B + y) +
                                     pow(xk,2)*(1 + 5*pow(B,2) - 7*pow(xq,2) - pow(xq,4) +
                                        B*(-6 + 3*pow(xq,2) - y) + 5*y + 12*pow(xq,2)*y -
                                        10*pow(y,2)) + 2*pow(xq,2)*
                                      (3 - 4*B + pow(B,2) - 9*y + 7*B*y + 5*pow(y,2)) -
                                     y*(4 - 6*B + 2*pow(B,2) - 11*y + 9*B*y + 7*pow(y,2))) +
                                  pow(A,2)*(rl*pow(xk,2) +
                                     pow(xq,2)*(-3 + B + 2*pow(xk,2) + 3*y) -
                                     y*(B + 3*(-1 + pow(xk,2) + y))) -
                                  A*(-1 + 3*pow(xk,2) - 3*pow(xk,4) + pow(xk,6) +
                                     3*pow(xq,4) - pow(xk,2)*pow(xq,4) -
                                     2*B*rl*(1 + pow(xk,2) - y) + 3*y - 6*pow(xk,2)*y +
                                     3*pow(xk,4)*y - 6*pow(xq,2)*y +
                                     4*pow(xk,2)*pow(xq,2)*y - 3*pow(xq,4)*y +
                                     6*pow(xq,2)*pow(y,2) - 2*pow(y,3) +
                                     pow(B,2)*(-1 + pow(xk,2) + y) +
                                     2*B*(-1 + pow(xk,2) + pow(xq,2))*
                                      (-1 + pow(xk,2) - pow(xq,2) + 2*y) +
                                     2*pow(rl,2)*(-4 + pow(xk,2) + pow(xq,2) + 3*y) +
                                     rl*(4 + 3*pow(xk,4) + 7*pow(xq,2) - 2*pow(xq,4) - 11*y -
                                        5*pow(xq,2)*y + 7*pow(y,2) +
                                         pow(xk,2)*(-7 - 3*pow(xq,2) + 12*y))));

 double sdfVhe1 = -4*pow(MK,6)*(A + 2*rl)*(pow(xq,2) - 2*B*pow(xq,2) +
                                             pow(B,2)*pow(xq,2) - 3*pow(xk,2)*pow(xq,2) +
                                             3*B*pow(xk,2)*pow(xq,2) + 2*pow(xk,4)*pow(xq,2) +
                                             pow(xq,4) - B*pow(xq,4) - 2*pow(xk,2)*pow(xq,4) +
                                             2*pow(rl,2)*(-1 + pow(xk,2) + pow(xq,2)) +
                                             rl*(pow(xk,4) + pow(xk,2)*(-1 + B + 3*pow(xq,2)) -
                                                2*pow(xq,2)*(1 + pow(xq,2) - 2*y)) - y + 2*B*y -
                                             pow(B,2)*y + 2*pow(xk,2)*y - 2*B*pow(xk,2)*y -
                                             pow(xk,4)*y - 4*pow(xq,2)*y + 4*B*pow(xq,2)*y +
                                             6*pow(xk,2)*pow(xq,2)*y - pow(xq,4)*y + 3*pow(y,2) -
                                             3*B*pow(y,2) - 3*pow(xk,2)*pow(y,2) +
                                             3*pow(xq,2)*pow(y,2) - 2*pow(y,3) -
                                             pow(A,2)*(-1 + pow(xk,2) + y) +
                                             A*(pow(xk,4) - (1 + 2*pow(xq,2) - 3*y)*(-1 + y) +
                                                B*(-1 + pow(xk,2) - pow(xq,2) + 2*y) +
                                                pow(xk,2)*(-2 + rl - 3*pow(xq,2) + 4*y)));

 double sdfVhe2 =(2*pow(MK,6)*(1 - B - rl)*rl*(A + 2*rl)*
                  (2*pow(B,2)*pow(xk,2) - 2*pow(xq,2) + B*pow(xq,2) -
                    pow(B,2)*pow(xq,2) + pow(xk,2)*pow(xq,2) -
                    B*pow(xk,2)*pow(xq,2) - pow(xq,4) + B*pow(xq,4) -
                    2*pow(rl,2)*(-1 + pow(xq,2)) +
                    A*(pow(xk,4) + B*(1 - 3*pow(xk,2) + pow(xq,2) - 2*y) +
                       pow(xk,2)*(-1 + 3*rl - pow(xq,2) - y) +
                       2*(pow(xq,2) - y)*(-1 + y)) - B*y + pow(B,2)*y +
                    pow(xk,2)*y + 2*B*pow(xk,2)*y - pow(xk,4)*y +
                    6*pow(xq,2)*y - 3*B*pow(xq,2)*y + pow(xq,4)*y -
                    pow(y,2) + 2*B*pow(y,2) - 4*pow(xq,2)*pow(y,2) +
                    pow(y,3) + pow(A,2)*(-1 + pow(xk,2) + y) +
                    rl*(-3*B*pow(xk,2) + pow(xk,4) +
                       2*pow(xq,2)*(1 + pow(xq,2) - 2*y) -
                        3*pow(xk,2)*(pow(xq,2) + y))))/(-1 + B + rl);
    
    
 double sdfVfeA = -2*pow(MK,6)*(pow(xq,2) - 3*B*pow(xq,2) +
                                  3*pow(B,2)*pow(xq,2) - pow(B,3)*pow(xq,2) -
                                  2*pow(xk,2)*pow(xq,2) + 4*B*pow(xk,2)*pow(xq,2) -
                                  2*pow(B,2)*pow(xk,2)*pow(xq,2) + pow(xk,4)*pow(xq,2) -
                                  B*pow(xk,4)*pow(xq,2) - pow(xq,6) + B*pow(xq,6) +
                                  pow(rl,3)*(8 - 6*pow(xk,2) - 8*pow(xq,2)) - y + 3*B*y -
                                  3*pow(B,2)*y + pow(B,3)*y + 3*pow(xk,2)*y -
                                  6*B*pow(xk,2)*y + 3*pow(B,2)*pow(xk,2)*y - 3*pow(xk,4)*y +
                                  3*B*pow(xk,4)*y + pow(xk,6)*y - 3*pow(xq,2)*y +
                                  6*B*pow(xq,2)*y - 3*pow(B,2)*pow(xq,2)*y +
                                  4*pow(xk,2)*pow(xq,2)*y - 4*B*pow(xk,2)*pow(xq,2)*y -
                                  pow(xk,4)*pow(xq,2)*y + 3*pow(xq,4)*y - 3*B*pow(xq,4)*y -
                                  pow(xk,2)*pow(xq,4)*y + pow(xq,6)*y + 3*pow(y,2) -
                                  6*B*pow(y,2) + 3*pow(B,2)*pow(y,2) -
                                  6*pow(xk,2)*pow(y,2) + 6*B*pow(xk,2)*pow(y,2) +
                                  3*pow(xk,4)*pow(y,2) - 3*pow(xq,4)*pow(y,2) -
                                  2*pow(y,3) + 2*B*pow(y,3) + 2*pow(xk,2)*pow(y,3) +
                                  2*pow(xq,2)*pow(y,3) + pow(A,3)*(-1 + pow(xk,2) + y) -
                                  pow(rl,2)*(2 + 3*pow(xk,4) - 12*pow(xq,2) - 6*pow(xq,4) +
                                     2*B*(-1 + 3*pow(xk,2) + pow(xq,2)) + 16*pow(xq,2)*y +
                                     pow(xk,2)*(-2 + 9*pow(xq,2) + y)) +
                                  pow(A,2)*(pow(xq,2)*(-3 + B + 2*pow(xk,2) + 3*y) +
                                     rl*(-4 + 5*pow(xk,2) + 4*y) - y*(B + 3*(-1 + pow(xk,2) + y))
                                     ) + rl*(-pow(xk,6) + 2*pow(xq,6) +
                                     pow(xq,4)*(-5 + 5*B + y) +
                                     pow(xk,4)*(2 + 2*B - 2*pow(xq,2) + 5*y) -
                                     2*pow(xq,2)*(3 - 5*B + 2*pow(B,2) - 10*y + 9*B*y +
                                        6*pow(y,2)) + y*
                                      (4 - 8*B + 4*pow(B,2) - 11*y + 11*B*y + 7*pow(y,2)) +
                                     pow(xk,2)*(-1 + 3*pow(B,2) + 11*pow(xq,2) + pow(xq,4) -
                                        9*y - 18*pow(xq,2)*y + 14*pow(y,2) +
                                        B*(-2 - 15*pow(xq,2) + 15*y))) -
                                  A*(-1 + 3*pow(xk,2) - 3*pow(xk,4) + pow(xk,6) +
                                     3*pow(xq,4) - pow(xk,2)*pow(xq,4) +
                                     2*pow(rl,2)*(-1 + 2*pow(xk,2) + pow(xq,2)) + 3*y -
                                     6*pow(xk,2)*y + 3*pow(xk,4)*y - 6*pow(xq,2)*y +
                                     4*pow(xk,2)*pow(xq,2)*y - 3*pow(xq,4)*y +
                                     6*pow(xq,2)*pow(y,2) - 2*pow(y,3) +
                                     pow(B,2)*(-1 + pow(xk,2) + y) +
                                     2*B*((-1 + pow(xk,2) + pow(xq,2))*
                                         (-1 + pow(xk,2) - pow(xq,2) + 2*y) +
                                        rl*(-2 + 6*pow(xk,2) - 2*pow(xq,2) + 4*y)) +
                                     rl*(4 + 7*pow(xk,4) + 5*pow(xq,2) - 2*pow(xq,4) - 15*y -
                                        3*pow(xq,2)*y + 11*pow(y,2) +
                                         pow(xk,2)*(-11 - 9*pow(xq,2) + 20*y))));

 double sdfVfeV=2*pow(MK,6)*(-pow(xq,2) + 3*B*pow(xq,2) -
                               3*pow(B,2)*pow(xq,2) + pow(B,3)*pow(xq,2) +
                               2*pow(xk,2)*pow(xq,2) - 4*B*pow(xk,2)*pow(xq,2) +
                               2*pow(B,2)*pow(xk,2)*pow(xq,2) - pow(xk,4)*pow(xq,2) +
                               B*pow(xk,4)*pow(xq,2) + 2*pow(xk,2)*pow(xq,4) -
                               2*B*pow(xk,2)*pow(xq,4) - pow(xq,6) + B*pow(xq,6) +
                               pow(rl,3)*(4 + 6*pow(xk,2) - 4*pow(xq,2)) + y - 3*B*y +
                               3*pow(B,2)*y - pow(B,3)*y - 3*pow(xk,2)*y +
                               6*B*pow(xk,2)*y - 3*pow(B,2)*pow(xk,2)*y + 3*pow(xk,4)*y -
                               3*B*pow(xk,4)*y - pow(xk,6)*y + 3*pow(xq,2)*y -
                               6*B*pow(xq,2)*y + 3*pow(B,2)*pow(xq,2)*y -
                               6*pow(xk,2)*pow(xq,2)*y + 6*B*pow(xk,2)*pow(xq,2)*y +
                               3*pow(xk,4)*pow(xq,2)*y + 3*pow(xq,4)*y -
                               3*B*pow(xq,4)*y - 3*pow(xk,2)*pow(xq,4)*y + pow(xq,6)*y -
                               3*pow(y,2) + 6*B*pow(y,2) - 3*pow(B,2)*pow(y,2) +
                               6*pow(xk,2)*pow(y,2) - 6*B*pow(xk,2)*pow(y,2) -
                               3*pow(xk,4)*pow(y,2) - 6*pow(xq,2)*pow(y,2) +
                               6*B*pow(xq,2)*pow(y,2) +
                               6*pow(xk,2)*pow(xq,2)*pow(y,2) - 3*pow(xq,4)*pow(y,2) +
                               4*pow(y,3) - 4*B*pow(y,3) - 4*pow(xk,2)*pow(y,3) +
                               4*pow(xq,2)*pow(y,3) - 2*pow(y,4) +
                               pow(A,3)*(-1 + pow(xk,2) + y) +
                               pow(rl,2)*(3*pow(xk,4) +
                                  B*(-2 - 4*pow(xk,2) + 2*pow(xq,2)) +
                                  2*(3 + pow(xq,4) - 4*y) +
                                  pow(xk,2)*(-10 - 3*pow(xq,2) + y)) +
                               pow(A,2)*(-3*pow(xq,2) + 2*pow(xk,2)*pow(xq,2) +
                                  B*(2 - 2*pow(xk,2) + pow(xq,2) - 3*y) + 3*y -
                                  3*pow(xk,2)*y + 3*pow(xq,2)*y - 3*pow(y,2) +
                                  rl*(-2 + 3*pow(xk,2) + 2*y)) +
                               rl*(pow(xk,6) + 2*pow(xq,6) + pow(xk,4)*(-2 + 2*B - y) +
                                  pow(xq,4)*(-3 + 3*B - y) -
                                  2*pow(xq,2)*(-3 + pow(B,2) + 2*y + 2*B*(1 + y)) +
                                  y*(-4 + 2*pow(B,2) + 3*y + pow(y,2) + B*(2 + 3*y)) +
                                  pow(xk,2)*(1 + pow(B,2) - 7*pow(xq,2) - 3*pow(xq,4) +
                                     5*y + 2*pow(xq,2)*y + B*(-2 - 9*pow(xq,2) + 3*y))) +
                               A*(-1 + 3*pow(xk,2) - 3*pow(xk,4) + pow(xk,6) +
                                  2*pow(xk,2)*pow(xq,2) - 2*pow(xk,4)*pow(xq,2) -
                                  3*pow(xq,4) + pow(xk,2)*pow(xq,4) +
                                  pow(rl,2)*(2 + 6*pow(xk,2) - 2*pow(xq,2)) + 3*y -
                                  6*pow(xk,2)*y + 3*pow(xk,4)*y + 6*pow(xq,2)*y -
                                  6*pow(xk,2)*pow(xq,2)*y + 3*pow(xq,4)*y - 6*pow(y,2) +
                                  6*pow(xk,2)*pow(y,2) - 6*pow(xq,2)*pow(y,2) +
                                  4*pow(y,3) + pow(B,2)*
                                   (-1 + pow(xk,2) - 2*pow(xq,2) + 3*y) +
                                  rl*(4 + 3*pow(xk,4) + 2*pow(xq,4) - y - 3*pow(y,2) -
                                     pow(xk,2)*(7 + pow(xq,2) + 2*y) + pow(xq,2)*(-5 + 3*y))
                                    + 2*B*(1 + pow(xk,4) + pow(xq,2) + pow(xq,4) +
                                     rl*(1 - 2*pow(xk,2) + pow(xq,2) - 2*y) - 3*y -
                                     3*pow(xq,2)*y + 3*pow(y,2) +
                                           pow(xk,2)*(-2 - 4*pow(xq,2) + 3*y))));
    
 double sdhe1he2 =(8*pow(MK,6)*(1 - B - rl)*rl*pow(A + 2*rl,2)*
                   (-2*pow(rl,2) - pow(xq,2) + B*pow(xq,2) - y + B*y +
                     pow(xk,2)*y + pow(xq,2)*y + pow(y,2) -
                     A*(-1 + pow(xk,2) + y) -
                    rl*(-2 + pow(xk,2) - 2*pow(xq,2) + 2*y)))/(-1 + B + rl);
    
 double sdhe1feA = -4*pow(MK,6)*(A + 2*rl)*(-2*pow(rl,2) + 10*pow(rl,3) -
                                              rl*pow(xk,2) - 3*pow(rl,2)*pow(xk,2) + 2*rl*pow(xk,4) +
                                              pow(xq,2) + rl*pow(xq,2) - 14*pow(rl,2)*pow(xq,2) -
                                              2*pow(xk,2)*pow(xq,2) + 6*rl*pow(xk,2)*pow(xq,2) -
                                              2*pow(xq,4) + 4*rl*pow(xq,4) +
                                              pow(B,2)*(6*rl + pow(xq,2) - y) - y - 3*rl*y +
                                              4*pow(rl,2)*y + 3*pow(xk,2)*y + 3*rl*pow(xk,2)*y -
                                              2*pow(xk,4)*y + pow(xq,2)*y - 3*rl*pow(xq,2)*y +
                                              2*pow(xq,4)*y + pow(y,2) + 3*rl*pow(y,2) -
                                              2*pow(xk,2)*pow(y,2) - 2*pow(xq,2)*pow(y,2) +
                                              pow(A,2)*(-1 + 2*B + pow(xk,2) + y) +
                                              A*(1 + 2*pow(B,2) + rl + 8*B*rl - 3*pow(xk,2) +
                                                 2*pow(xk,4) - 3*pow(xq,2) + 2*pow(xk,2)*pow(xq,2) +
                                                 3*B*(-1 + pow(xk,2) + pow(xq,2)) - rl*y + pow(xk,2)*y +
                                                 3*pow(xq,2)*y - pow(y,2)) +
                                              B*(20*pow(rl,2) + 2*pow(xq,4) +
                                                 pow(xq,2)*(-2 + 2*pow(xk,2) - y) -
                                                 y*(-2 + 3*pow(xk,2) + y) +
                                                 rl*(-6 + pow(xk,2) - pow(xq,2) + 3*y)));

 double sdhe1feV = 4*pow(MK,6)*(A + 2*rl)*(10*pow(rl,2) - 10*pow(rl,3) +
                                             rl*pow(xk,2) - 13*pow(rl,2)*pow(xk,2) - pow(xq,2) +
                                             3*rl*pow(xq,2) + 10*pow(rl,2)*pow(xq,2) + y - rl*y -
                                             8*pow(rl,2)*y - pow(xk,2)*y - rl*pow(xk,2)*y +
                                             3*pow(xq,2)*y - 5*rl*pow(xq,2)*y - 3*pow(y,2) +
                                             rl*pow(y,2) + 2*pow(xk,2)*pow(y,2) -
                                             2*pow(xq,2)*pow(y,2) + 2*pow(y,3) +
                                             pow(A,2)*(-1 + 2*B + pow(xk,2) + y) +
                                             pow(B,2)*(-6*rl - pow(xq,2) + y) +
                                             A*(-1 - 2*pow(B,2) - 4*pow(rl,2) + pow(xk,2) - pow(xq,2) +
                                                B*(3 + 4*rl - 3*pow(xk,2) + 3*pow(xq,2) - 6*y) +
                                                rl*(3 - 2*pow(xk,2) + 4*pow(xq,2) - 3*y) + 4*y -
                                                3*pow(xk,2)*y + pow(xq,2)*y - 3*pow(y,2)) +
                                             B*(-4*pow(rl,2) + rl*(6 - 13*pow(xk,2) + 9*pow(xq,2) -
                                                                     11*y) + pow(xq,2)*(2 - 3*y) + y*(-2 + pow(xk,2) + 3*y)));
 

 double sdhe2feA = (4*pow(MK,6)*(1 - B - rl)*rl*(A + 2*rl)*
                    (-6*pow(rl,2) + 6*pow(rl,3) + rl*pow(xk,2) +
                      3*pow(rl,2)*pow(xk,2) - pow(xq,2) + 3*rl*pow(xq,2) -
                      6*pow(rl,2)*pow(xq,2) + y - rl*y + 8*pow(rl,2)*y -
                      pow(xk,2)*y - rl*pow(xk,2)*y + 3*pow(xq,2)*y -
                      5*rl*pow(xq,2)*y - 3*pow(y,2) + rl*pow(y,2) +
                      2*pow(xk,2)*pow(y,2) - 2*pow(xq,2)*pow(y,2) +
                      2*pow(y,3) + pow(A,2)*(-1 + pow(xk,2) + y) +
                      pow(B,2)*(2*rl - pow(xq,2) + y) +
                      A*(-1 + pow(xk,2) - pow(xq,2) +
                         B*(1 - 2*rl + pow(xk,2) - pow(xq,2) - 2*y) + 4*y -
                         3*pow(xk,2)*y + pow(xq,2)*y - 3*pow(y,2) +
                         rl*(-1 + 2*pow(xk,2) + y)) +
                      B*(4*pow(rl,2) + pow(xq,2)*(2 - 3*y) +
                         y*(-2 + pow(xk,2) + 3*y) +
                         rl*(-2 + 3*pow(xk,2) - 7*pow(xq,2) + 5*y))))/(-1 + B + rl);

 

 double sdfeAfeV= -4*pow(MK,6)*(pow(A,3)*B +
                                  2*pow(A,2)*(rl*(-rl + pow(xq,2)) +
                                     B*(2*rl + pow(xq,2) - y)) -
                                  A*(pow(B,3) + 2*pow(B,2)*(-1 + rl + pow(xk,2) + y) +
                                     2*rl*(1 + 5*pow(rl,2) + pow(xk,2)*(-1 + y) - 2*y +
                                        pow(xq,2)*y + pow(y,2) +
                                        rl*(-4 + 4*pow(xk,2) - 5*pow(xq,2) + 3*y)) +
                                     B*(1 - 3*pow(rl,2) + 2*pow(xk,2)*(-1 + y) - 2*y +
                                        2*pow(xq,2)*y +
                                        2*rl*(-2 + 5*pow(xk,2) - 5*pow(xq,2) + 6*y))) -
                                  4*rl*(pow(B,3) + pow(B,2)*
                                      (-2 + 3*rl + 2*pow(xk,2) - pow(xq,2) + 2*y) +
                                     rl*(2 + 4*pow(rl,2) + pow(xq,2) + pow(xk,2)*(-2 + y) -
                                        3*y + pow(xq,2)*y + pow(y,2) +
                                        rl*(-6 + 5*pow(xk,2) - 4*pow(xq,2) + 4*y)) +
                                     B*(1 + 3*pow(rl,2) + pow(xq,2) + pow(xk,2)*(-2 + y) -
                                        2*y + pow(xq,2)*y +
                                        rl*(-5 + 7*pow(xk,2) - 5*pow(xq,2) + 6*y))));

 //pt

 pt *= pow(fk,2);




 //int
 inth1 *=real(h1);
 inth2 *=real(h2);
 intfA *= real(fA);
 intfV *= real(fV);
 inthe1 *= real(he1);
 inthe2 *= real(he2);
 intfeA *= real(feA);
 intfeV *= real(feV);


 //sd

 sdh1 *= norm(h1);
 sdh2 *= norm(h2);
 sdfA *= norm(fA);
 sdfV *= norm(fV);
 sdhe1 *= norm(he1);
 sdhe2 *= norm(he2);
 sdfeA *= norm(feA);
 sdfeV *= norm(feV);


 sdh1h2 *= real(h1*conj(h2));
 sdh1fA *= real(h1*conj(fA));
 sdh1fV *= real(h1*conj(fV));
 sdh1he1 *= real(h1*conj(he1));
 sdh1he2 *= real(h1*conj(he2));
 sdh1feA *= real(h1*conj(feA));
 sdh1feV *= real(h1*conj(feV));

 sdh2fA *= real(h2*conj(fA));
 sdh2he1 *= real(h2*conj(he1));
 sdh2he2 *= real(h2*conj(he2));
 sdh2feA *= real(h2*conj(feA));
 sdh2feV *= real(h2*conj(feV));

 sdfAfV *= real(fA*conj(fV));
 sdfAhe1 *= real(fA*conj(he1));
 sdfAhe2 *= real(fA*conj(he2));
 sdfAfeA *= real(fA*conj(feA));
 sdfAfeV *= real(fA*conj(feV));

 sdfVhe1 *= real(fV*conj(he1));
 sdfVhe2 *= real(fV*conj(he2));
 sdfVfeA *= real(fV*conj(feA));
 sdfVfeV *= real(fV*conj(feV));

 sdhe1he2 *= real(he1*conj(he2));
 sdhe1feA *= real(he1*conj(feA));
 sdhe1feV *= real(he1*conj(feV));

 sdhe2feA *= real(he2*conj(feA));

 sdfeAfeV *= real(feA*conj(feV));



 
 


 //Kahan sum the terms

 vector<double> Int_Kernel_To_Sum({inth1,inth2,intfA,intfV,inthe1,inthe2,intfeA,intfeV});
 vector<double> SD_Kernel_To_Sum({sdh1,sdh2,sdfA,sdfV,sdhe1,sdhe2,sdfeA,sdfeV,sdh1h2,sdh1fA,sdh1fV,sdh1he1,sdh1he2,sdh1feA,sdh1feV, sdh2fA,sdh2he1,sdh2he2,sdh2feA,sdh2feV, sdfAfV,sdfAhe1, sdfAhe2,sdfAfeA,sdfAfeV,sdfVhe1,sdfVhe2,sdfVfeA, sdfVfeV, sdhe1he2, sdhe1feA, sdhe1feV, sdhe2feA, sdfeAfeV});


 double SD_Kernel_To_Sum_naive= sdh1+sdh2+sdfA+sdfV+sdhe1+sdhe2+sdfeA+sdfeV+sdh1h2+sdh1fA+sdh1fV+sdh1he1+sdh1he2+sdh1feA+sdh1feV+ sdh2fA+sdh2he1+sdh2he2+sdh2feA+sdh2feV+sdfAfV+sdfAhe1+sdfAhe2+sdfAfeA+sdfAfeV+sdfVhe1+sdfVhe2+sdfVfeA+sdfVfeV+sdhe1he2+sdhe1feA+sdhe1feV+sdhe2feA+sdfeAfeV;


 double Int_Kernel= Kahan_sum(Int_Kernel_To_Sum);

 double SD_Kernel = Kahan_sum(SD_Kernel_To_Sum);

 double result;
 if(MODE=="PT") result= pt;
 else if(MODE=="INTERFERENCE") result= Int_Kernel;
 else if(MODE=="QUADRATIC") result= SD_Kernel;
 else if(MODE=="SD") result = SD_Kernel+Int_Kernel;
 else {
   vector<double> summ_vector({pt, Int_Kernel, SD_Kernel});
   result=Kahan_sum(summ_vector);
 }
 

 

 if(result <0 && ( MODE != "INTERFERENCE" && MODE != "SD")) {
   cout<<"Kernel extended is negative!"<<endl;
    cout<<"Printing info...."<<endl;
    cout<<"tot: "<<result<<endl;
    cout<<"point-like: "<<pt<<" "<<endl;
    cout<<"he1: "<<inthe1<<endl;
    cout<<"interference: "<<Int_Kernel<<endl;
    cout<<"Quadratic: "<<SD_Kernel<<"   Naive: "<<SD_Kernel_To_Sum_naive<<endl;
    cout<<"xk: "<<xk<<" xq: "<<xq<<endl;
    cout<<"#############"<<endl<<flush;
    cout<<sdh1<<" "<<sdh2<<" "<<sdfA<<" "<<sdfV<<" "<<sdhe1<<" "<<sdhe2<<" "<<sdfeA<<" "<<sdfeV<<" "<<sdh1h2<<" "<<sdh1fA<<" "<<sdh1fV<<" "<<sdh1he1<<" "<<sdh1he2<<" "<<sdh1feA<<" "<<sdh1feV<<" "<<sdh2fA<<" "<<sdh2he1<<" "<<sdh2he2<<" "<<sdh2feA<<" "<<sdh2feV<<" "<<sdfAfV<<" "<<sdfAhe1<<" "<<sdfAhe2<<" "<<sdfAfeA<<" "<<sdfAfeV<<" "<<sdfVhe1<<" "<<sdfVhe2<<" "<<sdfVfeA<<" "<<sdfVfeV<<" "<<sdhe1he2<<" "<<sdhe1feA<<" "<<sdhe1feV<<" "<<sdhe2feA<<" "<<sdfeAfeV;

    cout<<"###########"<<endl<<flush;
   
    crash("");
 }


 return result;
 
}











