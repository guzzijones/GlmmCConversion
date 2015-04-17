#ifndef GAMMADIST_H
#define GAMMADIST_H
#include "ConstantsForProb.h"
#include <cmath>
#include <stdexcept>
#include <vector>
#include "Polynomial.h"

class GammaDist:private ConstantsForProb{
   private:

   protected:

   public:

   static double powerSeries( double a, double b, double x ) ;
   static double incompleteBeta( double aa, double bb, double xx );
   static double logGamma(double x); 
   static double gamma(double x);
   static double stirlingFormula(double x) ;
   static double incompleteBetaFraction1( double a, double b, double x ) ;

   static double incompleteBetaFraction2( double a, double b, double x ) ;

};
#endif
