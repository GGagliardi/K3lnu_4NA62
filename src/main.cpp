#include "../include/numerics.h"
#include "../include/Num_integrate_l4_decay_rate.h"
#include "../include/ChPT_form_factors.h"

using namespace std;

/*
 * =========================================================================
 *  This file is ONLY an example driver.
 *
 *  Typical use case:
 *    - You will NOT normally use this main program.
 *    - Instead, include and call the functions
 *
 *        Num_Integrate_Decay_Rate(...)
 *        Get_squared_mat_el(...)
 *
 *      directly from your own code.
 *
 *  This file simply shows how to:
 *    (1) compute the total branching fractions
 *    (2) evaluate the squared matrix element at a chosen phase-space point
 * =========================================================================
 */

int main(int narg, char** argv)
{
  // This program takes no arguments. Simply run:
  //
  //     ./Klnull
  //
  // Extra arguments are treated as an error.
  if (narg != 1) {
    cout << "Usage: ./Klnull" << endl;
    exit(-1);
  }

  // ----------------------------------------------------------------------
  // 1) Total branching fractions: Num_Integrate_Decay_Rate(mode)
  // ----------------------------------------------------------------------
  //
  // The function
  //
  //     Num_Integrate_Decay_Rate(std::string mode)
  //
  // computes the branching fraction for the process
  //
  //     K^- → l'^+ l'^- l^- ν_l
  //
  // considering all possible leptonic channels.
  //
  // The argument `mode` selects which contribution to include:
  //
  //   "PT"
  //       Branching fraction in the point-like approximation
  //
  //   "TOTAL"
  //       Full branching fraction using ChPT form factors at O(p^4)
  //
  //   "QUADRATIC"
  //       Only terms quadratic in the structure-dependent form factors:
  //       ∝ H1^2, H2^2, FA^2, FV^2
  //
  //   "INTERFERENCE"
  //       Only terms interfering with the point-like contribution:
  //       ∝ fK * H1, fK * H2, fK * FA, fK * FV
  //       (with fK = kaon decay constant)
  //
  // Example:
  // This reproduces the result of hep-ph/9209261
  Num_Integrate_Decay_Rate("TOTAL");

  // ----------------------------------------------------------------------
  // Input parameters and kinematic cuts
  // ----------------------------------------------------------------------
  //
  // All numerical inputs (masses, CKM element Vus, etc.) are defined at
  // the top of:
  //
  //     Num_integrate_l4_decay_rate.cpp
  //
  // and can be freely modified.
  //
  // Setting
  //
  //     Switch_Bij_INPUT = true
  //
  // in Num_integrate_l4_decay_rate.cpp forces all inputs to match those
  // used in hep-ph/9209261.
  //
  // A lower cut in xk = sqrt(k^2)/MK is assumed. The values of these cuts
  // are controlled by:
  //
  //     cut_ee   : for e^+ e^- pairs
  //     cut_mumu : for μ^+ μ^- pairs
  //
  // also defined in Num_integrate_l4_decay_rate.cpp.

  // ----------------------------------------------------------------------
  // 2) Squared matrix element: Get_squared_mat_el(...)
  // ----------------------------------------------------------------------
  //
  // The function
  //
  //   Get_squared_mat_el(
  //       double xk,
  //       double xq,
  //       double y12,
  //       double y34,
  //       double phi,
  //       double rl,
  //       double rll,
  //       const function<complex<double>(double,double)> &H1,
  //       const function<complex<double>(double,double)> &H2,
  //       const function<complex<double>(double,double)> &FA,
  //       const function<complex<double>(double,double)> &FV,
  //       std::string mode
  //   )
  //
  // returns the quantity defined in Eq. (71) of
  //
  //     https://arxiv.org/pdf/2202.03833
  //
  // expressed in terms of the phase-space variables
  // (xk, xq, y12, y34, phi).
  //
  // Phase-space variables (MK = kaon mass):
  //
  //   xk = sqrt(k^2)/MK     invariant mass of l'^+ l'^- pair
  //   xq = sqrt(q^2)/MK     invariant mass of l^- ν_l pair
  //
  //   y12, y34, φ           defined in Eq. (68) of 2202.03833
  //
  // These variables can be computed directly from the four-momenta of
  // l'^+, l'^-, l^-, ν_l.
  //
  // IMPORTANT:
  // The measure transformation between
  //
  //     (dxk)(dxq)(dy12)(dy34)(dphi)
  //
  // and the four-body phase space dΦ4 is given in Eq. (67) of
  // 2202.03833.
  //
  // >>> If you want dΓ / dΦ4, you MUST multiply the result of
  //     Get_squared_mat_el by the appropriate Jacobian factor that
  //     converts (dxk dxq dy12 dy34 dphi) into dΦ4. <<<
  //
  // Lepton-mass ratios:
  //
  //   rl  = (m_l  / MK)^2
  //   rll = (m_l' / MK)^2

  // ----------------------------------------------------------------------
  // 3) Obtain ChPT form factors
  // ----------------------------------------------------------------------
  //
  // Compute_ChPT_form_factors fills the function objects H1, H2, FA, FV
  // with the O(p^4) ChPT predictions. You may replace them with any other
  // parametrisation if desired.
  function<complex<double>(double,double)> H1, H2, FA, FV;
  Compute_ChPT_form_factors(H1, H2, FA, FV);

  // ----------------------------------------------------------------------
  // 4) Example: K^- → e^+ e^- μ^- ν_μ
  // ----------------------------------------------------------------------
  //
  // Here:
  //   l'^+ l'^- = e^+ e^-
  //   l^-       = μ^-
  //
  // Mass ratios:
  double rll = pow(0.000510998950 / 0.493646, 2); // (m_e / MK)^2
  double rl  = pow(0.10565837     / 0.493646, 2); // (m_mu / MK)^2

  // Choose a kinematic point inside the allowed phase space
  double xk  = 0.4;
  double xq  = 0.25;

  // Allowed ranges (2202.03833):
  //
  //  -sqrt(1 - 4 rll / xk^2) <= y12 <= sqrt(1 - 4 rll / xk^2)
  //  -(1 - rl / xq^2)        <= y34 <=  (1 - rl / xq^2)
  //
  double y12 = 0.9 * sqrt(1 - 4 * rll / (xk * xk));
  double y34 = 0.9 * (1 - rl / (xq * xq));

  // φ ∈ [0, 2π)
  double phi = 0.5 * M_PI;

  // Evaluate the  squared matrix element in "TOTAL" mode
  double result =
      Get_squared_mat_el(xk, xq, y12, y34, phi, rl, rll,
                         H1, H2, FA, FV, "TOTAL");

  cout << "TEST for K^- -> e^+ e^- mu^- nu_mu:" << endl;
  cout << "xk:  " << xk  << endl;
  cout << "xq:  " << xq  << endl;
  cout << "y12: " << y12 << endl;
  cout << "y34: " << y34 << endl;
  cout << "phi: " << phi << endl;
  cout << "dGamma^5(xk,xq,y12,y34,phi)/(dxk*dxq*dy12*dy34*dphi) = "
       << result << endl;

  return 0;
}

