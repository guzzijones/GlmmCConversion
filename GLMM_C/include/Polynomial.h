#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include "ConstantsForProb.h"
#include <vector>

class Polynomial:private ConstantsForProb{

   public:

   static double p1evl( double x, const std::vector<double> &coef, int N) ;
   static double polevl( double x, const std::vector<double> &coef, int N) ;


};
#endif
