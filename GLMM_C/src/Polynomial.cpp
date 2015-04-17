#include "Polynomial.h"
double Polynomial::p1evl( double x, const std::vector<double> &coef, int N) 
   {
      double ans;
      ans = x + coef[0];
      for(int i=1; i<N; i++) { ans = ans*x+coef[i]; }
      return ans;
   }

double Polynomial::polevl( double x, const std::vector<double> &coef, int N) 
   {
      double ans;
      ans = coef[0];
      for(int i=1; i<=N; i++) ans = ans*x+coef[i];
      return ans;
   }
