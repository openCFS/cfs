#ifndef LIGHT_TOOLS_HH
#define LIGHT_TOOLS_HH

/** here he have tools which require not expensive includes, in contrast to ToolsFull.hh */

#include <cmath>
#include "General/defs.hh"
#include <def_use_cuda.hh>
#include "def_use_openmp.hh"

// C++ is so poor that there even is no really standard for pi - > why not????
#ifndef M_PI
  // avoid including boost for math constants
  #define M_PI 3.14159265358979323846 
#endif

namespace CoupledField 
{
//! Convert grad => rad
inline Double Grad2Rad( Double rad ) {
  return rad / 180.0 * M_PI;
}

//! Convert rad => grad
inline Double Rad2Grad( Double rad ) {
  return rad / M_PI * 180.0;
}

/** Compares if two doubles are close to each other */
inline bool close(Double d1, Double d2) { return std::abs(d1-d2) < 1e-6; }

/** clang complains about std::abs(c1 - c2) if i* is unsigned int but it would be used for templated Matrix::IsSymmetric()
 * @param eps not used but there for compatibelity*/
inline bool close(int i1, int i2, double eps = 0.0) { return i1 == i2; }
inline bool close(unsigned int i1, unsigned int i2, double eps = 0.0) { return i1 == i2; }

/** separate eps version to be faster! */
inline bool close(Double d1, Double d2, double eps) { return std::abs(d1-d2) < eps; }

/** Compared if two complex are close (if both the real and imaginary part are close) */
inline bool close(Complex c1, Complex c2) { return close(c1.real(), c2.real()) && close(c1.imag(), c2.imag()); }

inline bool close(Complex c1, Complex c2, double eps) { return close(c1.real(), c2.real(), eps) && close(c1.imag(), c2.imag(), eps); }

/** identifies numerical noise */
inline bool IsNoise(Double val) { return std::abs(val) < 1e-13; }

inline bool IsNoise(Complex val) { return IsNoise(val.real()) && IsNoise(val.imag()); }

inline bool IsNoise(int val) { return false; }

inline bool IsNoise(UInt val) { return false; }

/** http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// shortcuts for pow with int base - strange this is not in boost?!
inline double Pow2(double x) { return x * x; }
inline int Pow2(int x) { return x * x; }
inline unsigned int Pow2(unsigned int x) { return x * x; }
inline double Pow3(double x) { return x * x * x; }
inline int Pow3(int x) { return x * x * x; }
inline unsigned int Pow3(unsigned int x) { return x * x * x; }

/** Returns the sign of a value 
 * @return 0 if 0 or +/- 1 */ 
inline int Sign(int a) { return (a == 0) ? 0 : (a < 0 ? -1 : 1); } 
inline int Sign(Double a) { return (a == 0.0) ? 0 : (a < 0.0 ? -1 : 1); } 

/** Compares the sign for equality - note that we have three signs! -1, 0, 1 */ 
inline bool SameSign(Double a, Double b) { return (a < 0.0) == (b < 0.0); } 

inline bool UseOpenMP()
{
  #ifdef USE_OPENMP
    return true;
  #else
    return false;
  #endif
}

inline bool UseCuda()
{
  #ifdef USE_CUDA
    return true;
  #else
    return false;
  #endif
}

} // end of CoupledField

#endif // LIGHT_TOOLS_HH
