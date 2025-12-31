#include "../include/numerics.h"

using namespace std;


void D(int k) { 
  cout<<"D("<<k<<")"<<endl;
}

double eps(int k) {
  double epsilon= 1;
  for(int i=0;i<k;i++) epsilon/=10;

  return epsilon;
}

void crash(string Message) {
  cout << Message <<endl;
  exit(-1);
  return ;
}

