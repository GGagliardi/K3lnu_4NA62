#include "../include/Num_integrate_l4_decay_rate.h"

using namespace std;
const bool Switch_Bij_INPUT=true;
const double MkPh=  (Switch_Bij_INPUT)?0.493646:0.493677; // (0.493677 is the PDG value)
const double r_mu =  pow(0.10565837/MkPh,2);  //squared ratio (m_mu / m_K)^2
const double r_el = pow( 0.000510998950/MkPh,2); //1.1e-6 //squared ratio (m_el / mK)^2
const double Gf = 1.1663787*1e-5; //Fermi constant in GeV^-2
const double Vus=  (Switch_Bij_INPUT)?0.22:0.2243; // PDG-average: 0.2243
const double dVus = 0.00085; // PDG-err (not used)
const double Gamma =  (Switch_Bij_INPUT)?5.342174672593326e-17:5.317*1e-17; // Total decay width of the K^+ in GeV  error is 0.009 1e-17  (5.317 is the PDG value)
const double alpha = 1/137.036; //alpha-em constant
const double fkPh =   (Switch_Bij_INPUT)?(0.1136*sqrt(2.0)):0.1565 ; // ( 0.1565 is the ETMC value from https://arxiv.org/abs/2104.06747 (ETMC) which uses a definition of isoQCD compatible with ours
double cut_ee = pow(0.140/MkPh,2); //cut on e+e-
double cut_mumu= 4.0*r_mu;
double NORM= pow(Gf*Vus*alpha,2)*(1.0/Gamma)*(1.0/(pow(M_PI,4)*pow(2.0,12)));
const int SEED_VAL= 366432;


void display_results(char *title, double result, double error)
{
  //printf ("%s ==================\n", title);
  //printf ("result = % .16f\n", result);
  printf ("relative error   = % .16f %% \n ", 100*error/result);
}


double Get_squared_mat_el(double xk, double xq, double y12, double y34, double phi, double rl, double rll, const function<complex<double>(double,double)> &H1, const function<complex<double>(double,double)> &H2, const function<complex<double>(double,double)> &FA, const function<complex<double>(double,double)> &FV, string mode) {

  bool same_lepton=false;
  if(fabs(rl-rll) < 1e-16) same_lepton=true;

  //set input
  double mk=MkPh;
  double fk= fkPh;
 

  auto lambda = [&](double xk, double xq) -> double { return sqrt( pow((1.0 -pow(xk,2)-pow(xq,2)),2) -4.0*pow(xk*xq,2));};

  auto l_12 = [&](double xk) -> double { return sqrt( 1-4.0*rll/pow(xk,2));};
  
  auto l_34 = [&](double xq) -> double { return 1-(rl/pow(xq,2));};
  
  auto delta= [&](double xk, double xq) -> double { return pow(xk,2) -pow(xq,2);};
  
  auto delta34 = [&](double xq) -> double {return - rl/pow(xq,2);};
  
  auto y= [&](double xk, double xq, double y12, double y34, double phi) -> double { return ((1.0 - delta(xk,xq))*(1.0 -delta34(xq)) - lambda(xk,xq)*y34)/2.0 -rl;};
  
  auto A= [&](double xk, double xq, double y12, double y34, double phi) -> double {
    
    double E1 = ((1.0 + delta(xk,xq)) + lambda(xk,xq)*y12)/4.0;
    double E4 =  ((1.0 - delta(xk,xq))*(1.0 - delta34(xq)) - lambda(xk,xq)*y34)/4.0;
    
    double p1x = -(xk/2.0)*sqrt( pow(l_12(xk),2) - pow(y12,2));
    double p1y = (lambda(xk,xq) + (1.0 + delta(xk,xq))*y12)/4.0;
    
    double p4x = (xq/2.0)*sqrt( pow(l_34(xq),2) - pow(y34,2))*cos(phi);
    double p4y = -(lambda(xk,xq)*(1.0 - delta34(xq)) - (1.0 - delta(xk,xq))*y34)/4.0;
    
    return 2.0*( E1*E4 - p1x*p4x - p1y*p4y);
  };
  
  auto B= [&](double xk, double xq, double y12, double y34, double phi) -> double {
    
    
    double E2 =  ((1.0 + delta(xk,xq)) - lambda(xk,xq)*y12)/4.0; 
    double E3 = ((1.0 - delta(xk,xq))*(1.0 + delta34(xq)) + lambda(xk,xq)*y34)/4.0;
    
    double p2x = (xk/2.0)*sqrt( pow(l_12(xk),2) - pow(y12,2));
    double p2y = (lambda(xk,xq) -(1.0 + delta(xk,xq))*y12)/4.0;
    
    double p3x = -(xq/2.0)*sqrt( pow(l_34(xq),2) - pow(y34,2))*cos(phi);
    double p3y = -(lambda(xk,xq)*(1.0 + delta34(xq)) + (1.0 - delta(xk,xq))*y34)/4.0;
    
    
    return 2.0*(E2*E3 - p2x*p3x - p2y*p3y);
    
  };

  

    
  double a= A(xk,xq, y12, y34, phi);
  double b= B(xk,xq, y12, y34, phi);
  double Y= y(xk,xq, y12, y34, phi);

  double xk_prime= sqrt(2*rll + a);
  double xq_prime= sqrt(rl + b);

  if(xk_prime < 2.0*sqrt(rll) && same_lepton) crash("invalid xk_prime generated");
  if(xk < 2.0*sqrt(rll)) crash("invalid xk generated");
  
  
  
  double jacobian= 4.0*xk*xq;
  
  double square_ampl;

  if(same_lepton) {
    square_ampl= Compute_square_amplitude_extended(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), H1(xk_prime,xq_prime), H2(xk_prime, xq_prime), FA(xk_prime, xq_prime), FV(xk_prime,xq_prime), mk, fk,rl,xk, xq, a,b,Y, mode )*(jacobian*lambda(xk,xq)/mk);
  }
  else {
    square_ampl = 2.0*(jacobian*lambda(xk,xq)/mk)*Compute_square_amplitude_different_lepton(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), mk, fk,rl,rll,xk, xq, a,b,Y, mode );
  }

  return square_ampl*NORM;

}


double MonteCarlo_integration_extended_phase_space(const function<complex<double>(double,double)> &H1, const function<complex<double>(double,double)> &H2, const function<complex<double>(double,double)> &FA, const function<complex<double>(double,double)> &FV,double mk, double fk, string channel, double xk_inf,  bool same_lepton, int MySeed, string mode) {

  VPfloat bounds;
  double bound_l[5];
  double bound_u[5];
  double rll,rl;
  if(channel=="e+e-") {
    if(same_lepton) {
    rll= r_el;
    rl=  r_el;
    }
    else {
      rl=r_mu;
      rll=r_el;
    }
    bounds = VPfloat{ { xk_inf,1-sqrt(rl)}, {sqrt(rl), 1.0}, {-1.0, 1.0}, {-1.0,1.0}, {0.0, 2.0*M_PI}};

    bound_l[0] = xk_inf;
    bound_l[1]= sqrt(rl);
    bound_l[2] = -1.0;
    bound_l[3] = -1.0;
    bound_l[4] = 0.0;
    bound_u[0] = 1.0-sqrt(rl);
    bound_u[1] = 1.0;
    bound_u[2] = 1.0;
    bound_u[3] = 1.0;
    bound_u[4] = 2.0*M_PI;
    

   
  }
  else if(channel=="mu+mu-") {
    if(same_lepton) {
    rll=r_mu;
    rl= r_mu;
    }
    else {
      rl= r_el;
      rll=r_mu;
    }
    bound_l[0] = xk_inf; //2.0*sqrt(rll);
    bound_l[1]= sqrt(rl);
    bound_l[2] = -1.0;
    bound_l[3] = -1.0;
    bound_l[4] = 0.0;
    bound_u[0] = 1.0-sqrt(rl);
    bound_u[1] = 1.0;
    bound_u[2] = 1.0;
    bound_u[3] = 1.0;
    bound_u[4] = 2.0*M_PI;
    bounds = VPfloat{ { xk_inf,1-sqrt(rl)}, {sqrt(rl), 1.0}, {-1.0, 1.0}, {-1.0,1.0}, {0.0, 2.0*M_PI} };

    
  }
  else crash("In Montecarlo_integration_extended_phase_space channel: "+channel+" not yet implemented");

  //define lambda function to define energies and momenta in terms of the 5 integration variables, xk, xq, y12, y34 and phi. All dimensionful quantities are normalized over the Kaon Mass

  auto lambda = [&](double xk, double xq) -> double { return sqrt( pow((1.0 -pow(xk,2)-pow(xq,2)),2) -4.0*pow(xk*xq,2));};

  auto l_12 = [&](double xk) -> double { return sqrt( 1-4.0*rll/pow(xk,2));};

  auto l_34 = [&](double xq) -> double { return 1-(rl/pow(xq,2));};

  auto delta= [&](double xk, double xq) -> double { return pow(xk,2) -pow(xq,2);};

  auto delta34 = [&](double xq) -> double {return - rl/pow(xq,2);};

  auto y= [&](double xk, double xq, double y12, double y34, double phi) -> double { return ((1.0 - delta(xk,xq))*(1.0 -delta34(xq)) - lambda(xk,xq)*y34)/2.0 -rl;};
  
  auto A= [&](double xk, double xq, double y12, double y34, double phi) -> double {
	    
	    double E1 = ((1.0 + delta(xk,xq)) + lambda(xk,xq)*y12)/4.0;
	    double E4 =  ((1.0 - delta(xk,xq))*(1.0 - delta34(xq)) - lambda(xk,xq)*y34)/4.0;
	   
	    double p1x = -(xk/2.0)*sqrt( pow(l_12(xk),2) - pow(y12,2));
	    double p1y = (lambda(xk,xq) + (1.0 + delta(xk,xq))*y12)/4.0;

	    double p4x = (xq/2.0)*sqrt( pow(l_34(xq),2) - pow(y34,2))*cos(phi);
	    double p4y = -(lambda(xk,xq)*(1.0 - delta34(xq)) - (1.0 - delta(xk,xq))*y34)/4.0;

	    return 2.0*( E1*E4 - p1x*p4x - p1y*p4y);
	  };
  
  auto B= [&](double xk, double xq, double y12, double y34, double phi) -> double {


	    double E2 =  ((1.0 + delta(xk,xq)) - lambda(xk,xq)*y12)/4.0; 
	    double E3 = ((1.0 - delta(xk,xq))*(1.0 + delta34(xq)) + lambda(xk,xq)*y34)/4.0;
	    
	    double p2x = (xk/2.0)*sqrt( pow(l_12(xk),2) - pow(y12,2));
	    double p2y = (lambda(xk,xq) -(1.0 + delta(xk,xq))*y12)/4.0;
	    
	    double p3x = -(xq/2.0)*sqrt( pow(l_34(xq),2) - pow(y34,2))*cos(phi);
	    double p3y = -(lambda(xk,xq)*(1.0 + delta34(xq)) + (1.0 - delta(xk,xq))*y34)/4.0;
	   

	    return 2.0*(E2*E3 - p2x*p3x - p2y*p3y);
	      
	  };

  
  //perform Monte Carlo integration


  auto square_amplitude = [&](vector<double> const &par) -> double {
			      
			     if((signed)par.size() != 5) crash("squared_amplitude in Num_Integrate_Decay_Rate_extendend_phase_space only accepts vector<double> of size 5. "+to_string( (signed)par.size())+" provided.");

			    
			     double xk= par[0];
			     double xq= par[1];
			     if(xq > 1.0-xk) return 0.0;
			     double y12= par[2]*l_12(xk);
			     double y34= par[3]*l_34(xq);
			     double phi=par[4];

			
			     
			     double a= A(xk,xq, y12, y34, phi);
			     double b= B(xk,xq, y12, y34, phi);
			     double Y= y(xk,xq, y12, y34, phi);

			     double xk_prime= sqrt(2*rll + a);
			     double xq_prime= sqrt(rl + b);

			     if(xk_prime < 2.0*sqrt(rll) && same_lepton) crash("invalid xk_prime generated");
			     if(xk < 2.0*sqrt(rll)) crash("invalid xk generated");

			     if(xk_prime < xk_inf && same_lepton) return 0.0;

			     assert(xk_prime + xq_prime <= 1);
			  			  
			 
			     double jacobian= 4.0*xk*xq;

			     double square_ampl;

			     if(same_lepton) {

			      			     
			      square_ampl= Compute_square_amplitude_extended(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), H1(xk_prime,xq_prime), H2(xk_prime, xq_prime), FA(xk_prime, xq_prime), FV(xk_prime,xq_prime), mk, fk,rl,xk, xq, a,b,Y, mode )*(jacobian*l_12(xk)*l_34(xq)*lambda(xk,xq)/mk);

			     }

			     else {
			       
			      
			       square_ampl = 2.0*(jacobian*l_12(xk)*l_34(xq)*lambda(xk,xq)/mk)*Compute_square_amplitude_different_lepton(H1(xk,xq), H2(xk,xq), FA(xk,xq), FV(xk,xq), mk, fk,rl,rll,xk, xq, a,b,Y, mode );

			     }
			     
			     return square_ampl;
			   };
   


    //VEGAS GLS INTEGRATION

    double res_vegas, err_vegas;
    size_t calls = 500000; 

    const gsl_rng_type *T;
     gsl_rng *r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

   
    //to enter a generic seed mySeed
    gsl_rng_set(r, MySeed);

    gsl_monte_function_pp<decltype(square_amplitude)> Fp(square_amplitude);

    gsl_monte_function *G = static_cast<gsl_monte_function*>(&Fp);


    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(5);

    int warmup_calls= 50000; 

    gsl_monte_vegas_integrate (G, bound_l, bound_u, 5, warmup_calls, r, s,
                               &res_vegas, &err_vegas);

  

    //printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (G,bound_l, bound_u, 5, calls/5, r, s,
                                   &res_vegas, &err_vegas);
        
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.1);

    display_results("vegas final", res_vegas, err_vegas);

    gsl_monte_vegas_free(s);

    gsl_rng_free(r);


    return res_vegas*NORM;

    
}


void  Num_Integrate_Decay_Rate(string mode) {


  string MODE=mode;
 

  function<complex<double>(double,double)> H1;
  function<complex<double>(double,double)> H2;
  function<complex<double>(double,double)> FA;
  function<complex<double>(double,double)> FV;


  //set input
  double mk=MkPh;
  double fk= fkPh;
  Compute_ChPT_form_factors(H1,H2,FA,FV);
    
  //#################################    FROM NOW ON VEGAS IS USED TO PERFORM NUMERICAL INTEGRATION ########################
   

  //e+e-e+

    
    

  double res_eee_MC = MonteCarlo_integration_extended_phase_space(H1, H2, FA, FV, mk, fk, "e+e-", sqrt(cut_ee),1,  SEED_VAL, MODE);
  cout<<"Branching ratio for e+e-e+nu_e: "<<res_eee_MC<<endl;
  cout<<"################"<<endl;
  
  //e+e-mu+
  
  double res_ee_MC_extended= MonteCarlo_integration_extended_phase_space(H1, H2, FA, FV, mk, fk, "e+e-", sqrt(cut_ee), 0, SEED_VAL, MODE);
  cout<<"Branching ratio for Carlo e+e-mu+nu_mu: "<<res_ee_MC_extended<<endl;
  cout<<"################"<<endl;
  
  //mu+mu-mu+
  
  double res_mumumu_MC = MonteCarlo_integration_extended_phase_space(H1, H2, FA, FV, mk, fk, "mu+mu-", sqrt(cut_mumu), 1, SEED_VAL, MODE);
  cout<<"Branching ratio for mu+mu-mu+nu_mu: "<<res_mumumu_MC<<endl;
  cout<<"################"<<endl;
  
  //mu+mu-e+
  
  double res_mumu_MC_extended= MonteCarlo_integration_extended_phase_space(H1, H2, FA, FV, mk, fk, "mu+mu-", sqrt(cut_mumu), 0, SEED_VAL, MODE);
  cout<<"Branching ratio for mu+mu-e+nu_e: "<<res_mumu_MC_extended<<endl;
  cout<<"################"<<endl;
  

  return;
}


