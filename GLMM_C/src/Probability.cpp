#include "Probability.h"
#include "GammaDist.h"
#include <exception>
#include <stdexcept>

double Probability::studentT(double k, double t) { 
   if( k <= 0 ) throw new std::runtime_error("k <=0");
   if( t == 0 ) return( 0.5 );
   
   double cdf = 0.5 * GammaDist::incompleteBeta( 0.5*k, 0.5, k / (k + t * t) );
   
   if (t >= 0) cdf = 1.0 - cdf; // fixes bug reported by stefan.bentink@molgen.mpg.de
    
   return cdf;
}


