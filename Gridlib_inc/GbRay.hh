/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBRAY_HH
#define  GBRAY_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbVec3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbRay
{
public:
  //! Ray starts at origin P and has direction D
  GbRay() : origin_(), direction_() {}

  INLINE const GbVec3<T>& getOrigin() const { return origin_; }
  INLINE void setOrigin(const GbVec3<T>&o) { origin_ = o; }

  INLINE const GbVec3<T>& getDirection() const { return direction_; }
  INLINE void setDirection(const GbVec3<T>&d) { direction_ = d; }

  friend std::ostream& operator<<(std::ostream& s, const GbRay<T>& v);

private:
  GbVec3<T> origin_;  // P
  GbVec3<T> direction_;  // D
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbRay<T>& v)
{
  s<<typeid(v).name()<<": origin("<<v.getOrigin()[0]<<","<<v.getOrigin()[1]<<","<<v.getOrigin()[2]<<") ";
  s<<"direction("<<v.getDirection()[0]<<","<<v.getDirection()[1]<<","<<v.getDirection()[2]<<")"<<std::endl;
  return s;
}

#ifdef __GNUG__

//#include "GbRay.in"
//#include "GbRay.T"

#else

#pragma instantiate GbRay<float>
#pragma instantiate GbRay<double>
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbRay<float>& v)
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbRay<double>& v)

//#ifndef OUTLINE
//#include "GbRay.in"
//#endif

#endif  // g++

#endif  // GBRAY_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/01/02 14:56:32  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
