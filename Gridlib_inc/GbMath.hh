/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBMATH_HH
#define  GBMATH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <math.h>
#include <iostream>
#include <typeinfo>

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
   This class encapsulates numerical stuff, that should be available
   in libm.so on most platforms. However, on some platforms there are
   subtle differences (signatures !) between the C++-library calls that 
   wrap the C-library functions. So, in case of doubt, implement it here.

   The fast math functions provided, use an approximating polynomial.
   Note that for performance reasons, there is no range checking on 
   the arguments.

   The approximation formulas are from:
   "Handbook of Mathematical Functions",
   Edited by M. Abramowitz and I.A. Stegun.
   Dover Publications, Inc.
   New York, NY
   Ninth printing, December 1972

   For the inverse tangent calls, all approximations are valid for |t| <= 1.
   To compute ATAN(t) for t > 1, use ATAN(t) = PI/2 - ATAN(1/t).  For t < -1,
   use ATAN(t) = -PI/2 - ATAN(1/t).

   Speedups were measured on a Pentium II 400 Mhz with a release build of the
   code. Comparisons are to the calls sin, cos, tan, and atan in the standard
   math library.
*/

template <class T>
class GbMath
{
public:
  static const T PI;
  static const T PIHALF;
  static const T TWOPI;
  static const T DEG2RAD;
  static const T RAD2DEG;

  INLINE static T Abs(T v);

  //! Return -1 if the input is negative, 0 if the input is zero, and +1 if the input is positive.
  INLINE static int Sign (int i);
  INLINE static T Sign(T v);

  INLINE static T Sqrt(T v);

  INLINE static T Sin(T v);
  INLINE static T Cos(T v);

  //! Just wrappers, but the input value is clamped to [-1,1] to avoid silent returns of NaN.
  INLINE static T ASin(T v);
  INLINE static T ACos(T v);
  INLINE static T ATan2(T v1, T v2);

  //! Generate a random number in [0,1).  The random number generator may be seeded by a first call to UnitRandom with a positive seed.
  INLINE static T UnitRandom(T seed = 0.0);

  //! Generate a random number in [-1,1).  The random number generator may be seeded by a first call to SymmetricRandom with a positive seed.
  INLINE static T SymmetricRandom(T seed = 0.0);

  //! Generate a random number in [min,max).  The random number generator may be seeded by a first call to IntervalRandom with a positive seed.
  INLINE static T IntervalRandom (T fMin, T fMax, T seed = 0.0);

  //! fast sine function: assert:  0 <= fT <= PI/2 maximum absolute error = 1.6415e-04 speedup = 1.91
  INLINE static T FastSin0 (T fT);
  //! fast sine function: assert:  0 <= fT <= PI/2 maximum absolute error = 2.3279e-09 speedup = 1.40
  INLINE static T FastSin1 (T fT);
  //! fast cosine function: assert:  0 <= fT <= PI/2 maximum absolute error = 1.1880e-03 speedup = 2.14
  INLINE static T FastCos0 (T fT);
  //! fast cosine function: assert:  0 <= fT <= PI/2 maximum absolute error = 2.3082e-09 speedup = 1.47
  INLINE static T FastCos1 (T fT);
  //! fast tangens function: assert:  0 <= fT <= PI/4 maximum absolute error = 8.0613e-04 speedup = 2.51
  INLINE static T FastTan0 (T fT);
  //! fast tangens function: assert:  0 <= fT <= PI/4 maximum absolute error = 1.8897e-08 speedup = 1.71
  INLINE static T FastTan1 (T fT);
  //! fast inverse sine function: assert:  0 <= fT <= 1 maximum absolute error = 6.7626e-05 speedup = 2.59
  INLINE static T FastInvSin0 (T fT);
  //! fast inverse cosine function: assert:  0 <= fT <= 1 maximum absolute error = 6.7626e-05 speedup = 2.59
  INLINE static T FastInvCos0 (T fT);
  //! fast inverse tangens function: assert:  |fT| <= 1 maximum absolute error = 4.8830e-03 speedup = 2.14
  INLINE static T FastInvTan0 (T fT);
  //! fast inverse tangens function: assert:  |fT| <= 1 maximum absolute error = 1.1492e-05 speedup = 2.16
  INLINE static T FastInvTan1 (T fT);
  //! fast inverse tangens function: assert:  |fT| <= 1 maximum absolute error = 1.3593e-08 speedup = 1.50
  INLINE static T FastInvTan2 (T fT);
};

#ifdef __GNUG__

#include "GbMath.in"
//#include "GbMath.T"

#else

#pragma instantiate GbMath<float>
#pragma instantiate GbMath<double>

#ifndef OUTLINE
#include "GbMath.in"
#endif

#endif  // g++

#endif  // GBMATH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/08/16 16:53:17  prkipfer
| improved type safety for template parameter
|
| Revision 1.2  2001/06/15 08:32:25  prkipfer
| added compiler hints
|
| Revision 1.1  2001/01/02 14:51:23  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
