#include <cmath>

#include "mathfunctions.hh"

namespace CoupledField {

Double gammaln(Double xx)
{
// Internal arithmetic will be done in double precision
// A nicety that you can omit if five-figure accuracy is good enough.
  Double x,y,tmp,ser;
  static Double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091,
		      -1.231739572450155, 0.1208650973866179e-2,-0.5395239384953e-5};
  y=xx;
  x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (Integer j=0; j<=5; j++)
    ser += cof[j]/++y;

  return -tmp+log(2.5066282746310005*ser/x);
}
 
} // end of namespace
