/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBLINE_HH
#define  GBLINE_HH

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
class GbLine
{
public:
  //! Line is L(t) = P+t*D for any real-valued t.  D is not necessarily unit length.
  GbLine() : origin_(), direction_() {}

  INLINE const GbVec3<T>& getOrigin() const { return origin_; }
  INLINE void setOrigin(const GbVec3<T>&o) { origin_ = o; }

  INLINE const GbVec3<T>& getDirection() const { return direction_; }
  INLINE void setDirection(const GbVec3<T>&d) { direction_ = d; }

  friend std::ostream& operator<<(std::ostream& s, const GbLine<T>& v);

private:
  GbVec3<T> origin_;  // P
  GbVec3<T> direction_;  // D
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbLine<T>& v)
{
  s<<typeid(v).name()<<": origin("<<v.getOrigin()[0]<<","<<v.getOrigin()[1]<<","<<v.getOrigin()[2]<<") ";
  s<<"direction("<<v.getDirection()[0]<<","<<v.getDirection()[1]<<","<<v.getDirection()[2]<<")"<<std::endl;
  return s;
}

#ifdef __GNUG__

//#include "GbLine.in"
//#include "GbLine.T"

#else

#pragma instantiate GbLine<float>
#pragma instantiate GbLine<double>
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbLine<float>& v)
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbLine<double>& v)

//#ifndef OUTLINE
//#include "GbLine.in"
//#endif

#endif  // g++

#endif  // GBLINE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/01/02 14:57:28  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
