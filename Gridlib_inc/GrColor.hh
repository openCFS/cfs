/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRCOLOR_HH
#define GRCOLOR_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>
#include <iostream>

#include "GbTypes.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GrColor
{
public:
  // construction (components in range [0,1])
  INLINE GrColor();
  INLINE GrColor(T fR, T fG, T fB);
  INLINE GrColor(const GrColor<T>& rkColor);

  // access color C as C[0] = C.r, C[1] = C.g, C[2] = C.b
  INLINE T operator[] (int i) const;

  void clamp ();
  void scaleByMax ();

  INLINE GrColor<T>& operator= (const GrColor<T>& rkColor);
  INLINE GbBool operator== (const GrColor<T>& rkColor) const;
  INLINE GbBool operator!= (const GrColor<T>& rkColor) const;

  INLINE GrColor<T> operator+ (const GrColor<T>& rkColor) const;
  INLINE GrColor<T> operator- (const GrColor<T>& rkColor) const;
  INLINE GrColor<T> operator* (const GrColor<T>& rkColor) const;
  INLINE GrColor<T> operator/ (const GrColor<T>& rkColor) const;
  INLINE GrColor<T> operator- () const;
  INLINE GrColor<T>& operator+= (const GrColor<T>& rkColor);
  INLINE GrColor<T>& operator-= (const GrColor<T>& rkColor);
  INLINE GrColor<T>& operator*= (const GrColor<T>& rkColor);
  INLINE GrColor<T>& operator/= (const GrColor<T>& rkColor);

  friend GrColor<T> operator* (T fScalar, const GrColor<T>& rkColor);

  static const GrColor<T> BLACK; // = (0,0,0) 
  static const GrColor<T> WHITE; // = (1,1,1)

  friend std::ostream& operator<<(std::ostream& s, const GrColor<T>& v);

  T r_, g_, b_;
};

template <class T>
GrColor<T> 
operator* (T fScalar, const GrColor<T>& rkColor)
{
  return GrColor<T>(fScalar*rkColor[0],fScalar*rkColor[1],fScalar*rkColor[2]);
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GrColor<T>& v)
{
  s<<typeid(v).name()<<": color("<<v[0]<<","<<v[1]<<","<<v[2]<<")"<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GrColor.in"
#include "GrColor.T"

#else

#pragma instantiate GrColor<float>
#pragma instantiate GrColor<double>
#pragma instantiate GrColor<float> operator* (float fScalar, const GrColor<float>& rkColor)
#pragma instantiate GrColor<double> operator* (double fScalar, const GrColor<double>& rkColor)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GrColor<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GrColor<double>&)

#ifndef OUTLINE
#include "GrColor.in"
#endif

#endif  // g++


#endif // GRCOLOR_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.4  2001/06/19 16:30:20  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.3  2001/06/18 11:13:37  prkipfer
| made operator/ inline
|
| Revision 1.2  2001/06/15 09:02:30  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:29  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
