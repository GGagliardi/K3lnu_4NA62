#ifndef __Num_integrate_l4_decay_rate__
#define __Num_integrate_l4_decay_rate__



#include "numerics.h"
#include "kernel.h"
#include "ChPT_form_factors.h"


using namespace std;


void display_results (char *title, double result, double error);


double Get_squared_mat_el(double xk, double xq, double y12, double y34, double phi, double rl, double rll, const function<complex<double>(double,double)> &H1, const function<complex<double>(double,double)> &H2, const function<complex<double>(double,double)> &FA, const function<complex<double>(double,double)> &FV, string mode);

double MonteCarlo_integration_extended_phase_space(const function<complex<double>(double,double)> &H1, const function<complex<double>(double,double)> &H2, const function<complex<double>(double,double)> &FA, const function<complex<double>(double,double)> &FV, double mk, double fk,string channel, double xk_inf,  bool same_lepton, int MySeed, string mode);

void Num_Integrate_Decay_Rate(string mode);







#endif
