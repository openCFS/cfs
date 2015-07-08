#include <cmath>
#include <iostream>

#include "mathfunctions.hh"
#include "MatVec/Vector.hh"
#include "Utils/tools.hh"
#include "General/Environment.hh"
#include "MatVec/Matrix.hh"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>


namespace CoupledField {


//! =======================================================================
//! List of mathematical functions, primarily for the use in mathParser
//! to generate signals
//! =======================================================================


  //! Generate a sinus burst signal
  Double SinBurst( Double freq, Double numPeriods,
                Double nPerFadeIn, Double nPerFadeOut,
                Double t) {

    Double ret = 0.0;

    Double T = 1.0 / freq;

    if( t < nPerFadeIn * T ) {
      ret = sin( 2 * M_PI * freq * t) *
            pow( sin( 2 * M_PI * freq * t / 4 / nPerFadeIn), 2);
    }
    else if( t < (numPeriods - nPerFadeOut) * T )  {
      ret = sin( 2 * M_PI * freq * t);
    }
    else if( t < numPeriods * T) {
      ret = sin( 2 * M_PI * freq * t) *
            pow( sin( 2 * M_PI * freq *
                (t - numPeriods * T) / 4 / nPerFadeOut), 2);
    }

    return ret;
  }

  //! return a fade-in function
  Double FadeIn( Double fadeInTime, Double mode, Double t) {

    Double ret = 0.0;

    //! Helper Parameter for the exp functions
    Double T0 = 0.0;

    //! amount in percent that will be reached after fadeInTime
    Double perCent = 0.999;



    if( mode == 1) {
      //! using sin^2
      if ( t < fadeInTime ) {
        ret = sin(2 * M_PI * t / fadeInTime / 4) * sin(2 * M_PI * t / fadeInTime / 4);
      }
      else {
        ret = 1.0;
      }
    }
    else if( mode == 2) {
      //! using exp
      T0 = fadeInTime / ( -1.0 * log( 1 - perCent) );
      ret = 1 - exp( -(t/T0) );
    }
    else if( mode == 3) {
      //! using exp.^2
      T0 = fadeInTime / sqrt( -1.0 * log( 1 - perCent) );
      ret = 1 - exp( -(t/T0) * (t/T0) );
    }

    return ret;
  }

  //! Generate a spike signal
  Double Spike( Double duration, Double t ) {
    Double ret=0.0, riseT=duration/2;

    
    if (t <= duration)
      ret = (t <= riseT) ? t/riseT : 2.0-t/riseT;
    
    if (t < 0 )
       ret = 0.0;
    return ret;
  }

  Double Chirp(Double sweepTime, Double startFreq, Double stopFreq, 
               Double nPerFadeIn, Double nPerFadeOut, Double t) {

    if( t > sweepTime)
      return 0.0;

    Double f00 = (stopFreq + startFreq) / 2.0;

    Double t00 = nPerFadeIn / f00;
    Double t01 = sweepTime - (nPerFadeOut / f00);
    Double b = (stopFreq - startFreq) / 2.0 / sweepTime;
    Double h = 1.0;
    Double g = 1.0;
    if( t < t00 ) {
      h =  sin( 2 * M_PI * f00 * t / (4.*nPerFadeIn)) *
           sin( 2 * M_PI * f00 * t / (4.*nPerFadeIn));
    }
    if( t > t01 ) {
      g = 1 - sin(2 * M_PI * f00 * (t01-t) / (4.*nPerFadeOut)) *
              sin(2 * M_PI * f00 * (t01-t) / (4.*nPerFadeOut));
    }
    Double c = sin(2 * M_PI * ( b * t * t + startFreq*t));
    return h * c * g;

  }
  
  
  //! Generate a band-filtered spike signal
  Double SpikeBPF( Double cutOff, Double slewRate, Double t ) {
    return 0;
  }

  //! Generate a square pulse signal
  Double CosPulseComb( Double freq, Double pulseWidth, Double t ) {

    Double ret = 0.0;

    Double T = 1.0 / freq;

    if( (Mod(t,T) < pulseWidth/2.0) || (Mod(t,T) > T - pulseWidth/2.0)  )  {
      ret = 1 + cos( 2 * M_PI * t / pulseWidth);
    }

    return ret;
  }



  //! Generate a square pulse signal
  Double SquarePulse( Double freq, Double numPeriods, Double biPolarType,
                      Double pulseWidth, Double riseTime, Double t ) {

    PolarType pT = (PolarType) (Integer) biPolarType;

    Double ret = 0.0;

    Double T = 1.0 / freq;

    if( t > numPeriods * T )  {
      ret = sin( 2 * M_PI * freq * t );
    }

    switch ( pT ) {

      case UNI_POLAR:
        if( fmod(t,T) < riseTime ) {
          ret = 1.0/riseTime*fmod(t,T);
        }
        else if( fmod(t,T) < pulseWidth/100*T) {
          ret = 1.0;
        }
        else if( fmod(t,T) < pulseWidth/100*T+riseTime)
          ret = 1.0 - 1.0/riseTime*(fmod(t,T)-pulseWidth/100*T);
        else {
          ret = 0.0;
        }
        break;

      case BI_POLAR:
        if( fmod(t,T) < riseTime/2 ) {
          ret = 2.0/riseTime*fmod(t,T);
        }
        else if( fmod(t,T) < pulseWidth/100*T-riseTime/2) {
          ret = 1.0;
        }
        else if( fmod(t,T) < pulseWidth/100*T+riseTime/2) {
          ret = 1.0 - 2.0/riseTime*(fmod(t,T)-pulseWidth/100*T+riseTime/2);
        }
        else if( fmod(t,T) < T-riseTime/2) {
          ret = -1.0;
        }
        else {
          ret= -1.0 + 2.0/riseTime*(fmod(t,T)+riseTime/2-T);
        }
        break;

      default:
        EXCEPTION( "Polar Type '" <<  pT << "' not known!" );
    }

    if( t > numPeriods * T )  {
      ret = 0.0;
    }

    return ret;
  }

  //! Modulo function
  Double Mod( Double x, Double m ) {

    return x-floor(x/m)*m;
  }

  //! Generate a Gauss shaped signal
  //! normVal: 1   the returned max val equals to 1
  //!          0   the integral value of the returned values equals 1

  Double Gauss( Double mue, Double sigma, Double normVal, Double x ) {
    Double help = (x - mue) / sigma;
    if( normVal == 1) {
      return exp( -0.5 * help * help );
    }
    else {
      return 1.0 / ( sigma * sqrt( 2 * M_PI) ) *
             exp( -0.5 * help * help );
    }
  }

  //! Calculate cylindric bessel function of first kind
  Double BesselCylJ( Double x, Double v ) {
    return boost::math::cyl_bessel_j(v, x);
  }

  //! Calculate cylindric bessel function of second kind
  Double BesselCylY( Double x, Double v ) {
    return boost::math::cyl_neumann(v, x );
  }

  //! Calculate spherical bessel function of first kind
  Double BesselSphJ( Double x, Double v ) {
    return boost::math::sph_bessel( (UInt) v, x );
  }

  //! Calculate spherical bessel function of second kind
  Double BesselSphY( Double x, Double v ) {
    return boost::math::sph_neumann( (UInt) v, x );
  }

  //! Calculate cylindric Hankel function of first kind
  Complex HankelCyl1( Double x, Double v ) {
    return boost::math::cyl_hankel_1(v, x);
  }

  //! Calculate cylindric Hankel function of second kind
  Complex HankelCyl2( Double x, Double v ) {
    return boost::math::cyl_hankel_2(v, x);
  }

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
 
  /* this method determines an approximation of the eigenvalues lambda of a 
     symmetric positive definite matrix. (lambda * A = lambda x)
     Input: 
     mat: the matrix A whose eigenvalues are wanted. 
     Remark that, after the computation, the result of Givens rotations 
     is overwritten on the array "mat". 
     eps: error tolerance. 
     If the sum of the square of the non-diagonal elements is less than 
     the one of the initial matrix, the iteration is terminated.  
  */
  void eigenValues(Matrix<Double> & mat, Double eps, Vector<Double> & eigen)
  {
    Integer i, j;
    Integer k, kmax;
		// initialize or receive compiler warning
    Integer l_conv = 0;
    Double a2, eps2, dkmax;
    Double n, n2;
    Integer l_sort=1;

    Integer ndim = mat.GetNumRows();
  
    a2 = 0.0;
    for (i=0; i<ndim; ++i)
      for (j=0; j<i-1; ++j)
        a2 += mat[i][j] * mat[i][j];
    a2 *= 2.0;
   
    n = (Double)ndim;
    n2 = n * n;
    eps2 = eps * eps;
    dkmax = std::log(eps)/std::log((n2-n-2)/(n2-n));
    kmax = (Integer)std::ceil(dkmax);

   
    for (k=1; k<kmax; ++k)
      {
        givensRotation(ndim, mat);
        terminationCriterion(ndim, mat, a2, eps2, &l_conv);

        if (l_conv==1)
          break;
      }
    if (l_conv==0) 
      std::cout<<"\n Note: The Jacobi method did not converge!.\n";
  
    for (i=0; i<ndim; ++i)
      eigen[i] = mat[i][i];
    sortArray(ndim, l_sort, eigen);
  }


  void givensRotation(Integer ndim, Matrix<Double> & mat)
  {
    Matrix<Double> b(ndim,ndim);
    
    Double a2, max_a2;
		// initialize or receive compiler warning
    Integer i, j, p = 0, q = 0;
    Double z, t, c, s, u;
  
    max_a2 = 0.0;
    for (i=0; i<ndim; ++i)
      {
        for (j=0; j<=i-1; ++j)
          {
            a2 = mat[i][j];
            a2 = a2 * a2;
            if (a2 > max_a2)
              {
                p = i; 
                q = j;
                max_a2 = a2;
              }
          }
      }
   
    z = 0.5 * (mat[q][q] - mat[p][p]) / mat[p][q];
    t = std::fabs(z) + std::sqrt(1.0 + z*z);
    if (z < 0.0) t = - t;
    t = 1.0 / t;
    c = 1.0 / std::sqrt(1.0 + t*t);
    s = c * t;
    u = s / (1.0 + c);
   
    for (i=0; i<ndim; ++i)
      for (j=0; j<ndim; ++j)
        b[i][j] = mat[i][j];
  
    b[p][p] = mat[p][p] - t * mat[p][q];
    b[q][q] = mat[q][q] + t * mat[p][q];
    b[p][q] = b[q][p] = 0.0;
    for (j=0; j<ndim; ++j)
      if ((j!=p) && (j!=q))
        {
          b[p][j] = b[j][p] = mat[p][j] - s * (mat[q][j] + u*mat[p][j]);
          b[q][j] = b[j][q] = mat[q][j] + s * (mat[p][j] - u*mat[q][j]);
        }
  
    for (i=0; i<ndim; ++i)
      for (j=0; j<ndim; ++j){
        mat[i][j] = b[i][j];


      }
  
  }
  /*
    termination criterion
  
    The iteration is terminated if F(k) <= eps^2 * F(0), 
    where "eps" is the relative error tolerance, 
    F(k) is the sum of all the non-diagonal elements squared of the matrix 
    after the k-th iteration (Given's rotation).

    (input) 
    ndim:    the matrix size
    mat:     double [][NDIM] : the matrix after each Givens rotation. 
    a2:      double : the sum of all the non-diagonal elements squared 
    of the matrix A.
    eps2:    the square of the error tolerance "eps", i.e., eps^2.

    (output) 
    l_conv int * : l_conv = 1 if converged and 
    l_conv = 1 if not yet converged. 
  */
  void terminationCriterion(Integer ndim, Matrix<Double> & mat, Double a2, Double eps2, Integer *l_conv)
  {
    Double a_nd2;
    Integer i, j;
  
    *l_conv = 0;
  
    a_nd2 = 0.0;
    for (i=0; i<ndim; ++i)
      for (j=0; j<i-1; ++j)
        a_nd2 += mat[i][j] * mat[i][j];
    a_nd2 *= 2.0;
  
    if (a_nd2/a2 < eps2) *l_conv = 1;

  }
  /*
    sorting of a 1-dimensional array d(1), d(2), ..., d(n)
    Input: l_sort decides about the order, i.e. l_sort==1 -> sorting downwards
    lsort == 0 -> sorting upwards   */

  void sortArray(Integer ndim, Integer l_sort, Vector<Double> & d)
  {
    Double dv;
    Integer k, i;
 
    if (l_sort == 0)
      {
        for (k=0; k<ndim-1; ++k)
          for (i=k+1; i<ndim; ++i)
            if (d[i] > d[k])
              {
                dv = d[k];
                d[k] = d[i];
                d[i] = dv;
              }
      }
    if (l_sort == 1)
      {
        for (k=0; k<ndim-1; ++k)
          for (i=k+1; i<ndim; ++i)
            if (d[i] < d[k])
              {
                dv = d[k];
                d[k] = d[i];
                d[i] = dv;
              }
      }
  }

  /* Commented out due to refactoring.
   * Can be removed, if implementation in Vector class works.
  Double Normalize(Vector<Double>& vec)
  {
    Double norm = vec.NormL2();
    if(norm < 1e-20)
        vec *= 0.0;
    else
        vec *= 1.0 / norm;

    return norm;
  }
    
  void CrossProd(const Vector<Double>& a,
                 const Vector<Double>& b,
                 Vector<Double> &result)
  {
#ifdef CHECK_INDEX
    if (a.GetSize() != 3)
      EXCEPTION("Incompatible vector dimensions for CrossProd()");
    if (b.GetSize() != 3)
      EXCEPTION("Incompatible vector dimensions for CrossProd()");
#endif

    result.Resize(3);

    result[0] =  a[1] * b[2] - a[2] * b[1];
    result[1] = -a[0] * b[2] + a[2] * b[0];
    result[2] =  a[0] * b[1] - a[1] * b[0];
  }
    
  bool CoLinear(const Vector<Double> &a, const Vector<Double> &b)
  {
      Vector<Double> c;
      CrossProd(a, b, c);
      return (fabs(c.NormL2()) < 1e-12);
  }

  void GetBarycentricCoords(const Vector<Double>& p1,
                            const Vector<Double>& p2,
                            const Vector<Double>& p3,
                            const Vector<Double>& p,
                            Vector<Double>& b)
  {
    Vector<Double> u, w, v, diff, normal, n;
    Double area_factor;

    w = p2 - p1;
    //    v = p1 - p3;
    u = p3 - p2;

    CrossProd(w, u, normal);
    area_factor = 1.0 / Normalize(normal);

    diff = p - p1;
    CrossProd(w, diff, n);

    normal.Inner(n, b[2]);
    b[2] *= area_factor;

    diff = p - p2;
    CrossProd(u, diff, n);
    normal.Inner(n, b[0]);
    b[0] *= area_factor;

    //diff = p - p3;
    //v.Cross(diff, n);
    //normal.Inner(n, b[1]);
    //b[1] *= area_factor;

    b[1] = 1 - b[0] - b[2];
  }
*/

} // end of namespace
