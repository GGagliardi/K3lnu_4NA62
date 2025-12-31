#ifndef __numerics__
#define __numerics__

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <map>
#include <string>
#include <functional>
#include <stdarg.h>
#include <numeric>
#include <sstream>
#include <cassert>
#include <utmpx.h>
#include <bits/stdc++.h> 
#include <ctime>
#include <chrono>
#include <complex>
#include <iomanip>
#include <math.h>
#include <stdexcept>
#include <utility>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

using namespace std;

typedef vector<int> Vint;
typedef vector<vector<int>> VVint;
typedef pair<int,int> Pint;
typedef pair<double,double> Pfloat;
typedef vector<pair<double,double>> VPfloat;
typedef vector<vector<pair<double,double>>> VVPfloat;
typedef vector<pair<int,int>> VPint;
typedef vector<double> Vfloat;
typedef vector<vector<double>> VVfloat;
typedef vector<long double> Vdouble;
typedef vector<vector<pair<int,int>>> VVPint;
typedef vector<vector<vector<vector<double>>>> VVVVfloat;
typedef vector<vector<vector<double>>> VVVfloat;
typedef vector<vector<vector<pair<double,double>>>> VVVPfloat;
typedef vector<vector<vector<vector<pair<double, double>>>>> VVVVPfloat;




//wrapper for lambdas with capture to gsl_monte_function
template< typename F >  class gsl_monte_function_pp : public gsl_monte_function {
public:
  gsl_monte_function_pp(const F& func) : _func(func) {
    f = &gsl_monte_function_pp::invoke;
    params=this;
    dim= 5;
  }
private:
  const F& _func;
  static double invoke(double x[], size_t dim, void *params) {
    vector<double> xv;
    for(int i=0; i<(signed)dim;i++) xv.push_back(x[i]);
    return static_cast<gsl_monte_function_pp*>(params)->_func(xv);
  }
};


//wrapper for lambdas with capture to gsl_function
template< typename F >  class gsl_function_pp : public gsl_function {
public:
  gsl_function_pp(const F& func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params=this;
  }
private:
  const F& _func;
  static double invoke(double x,void *params) {
    return static_cast<gsl_function_pp*>(params)->_func(x);
  }
};

//wrapper for lambdas with capture to gsl_function_fdf
template< typename F, typename F2, typename fdF >  class gsl_function_fdf_pp : public gsl_function_fdf {
public:
  gsl_function_fdf_pp(const F& func, const F2& df_func, const fdF& fdf_func) : _func(func), _df_func(df_func), _fdf_func(fdf_func) {
    f = &gsl_function_fdf_pp::invoke;
    df = &gsl_function_fdf_pp::invoke_df;
    fdf= &gsl_function_fdf_pp::invoke_fdf;
    params=this;
  }
private:
  const F& _func;
  const F2& _df_func;
  const fdF& _fdf_func;
  static double invoke(double x,void *params) {
    return static_cast<gsl_function_fdf_pp*>(params)->_func(x);
  }
   static double invoke_df(double x,void *params) {
    return static_cast<gsl_function_fdf_pp*>(params)->_df_func(x);
  }
  static void invoke_fdf(double x, void *params, double* f, double *df) {
    return static_cast<gsl_function_fdf_pp*>(params)->_fdf_func(x,f,df);
  }
};



void D(int k);
void crash(string Message);
double eps(int k);




template <typename T>
string to_string_with_precision(T a_value, const int n)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}

template <typename T>
void printV(const vector<T> &A, string B, bool mode) {
  int size = A.size();
  cout.precision(10);
  cout << B << endl;
  if(mode ==0 )for(int i = 0;i < size; i++) cout<< A[i] << " ";
  else for(int i=0;i<size;i++) cout<<i<<"  "<<A[i]<<endl;
  cout << endl;
  return;

}


template <typename T>
T Kahan_sum(const vector<T> &input) {

  T sum= 0.;
  T c = 0.;

  for (auto & val: input) {
    T y = val- c;
    T t = sum + y;
    c =(t-sum) -y;
    sum = t;
  }
   
  return sum;

};



#endif


