/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBSEGMENT_HH
#define  GBSEGMENT_HH

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
class GbSegment
{
public:
  //! Segment has endpoints P and P+D
  GbSegment() : origin_(), direction_() {}

  INLINE const GbVec3<T>& getOrigin() const { return origin_; }
  INLINE void setOrigin(const GbVec3<T>&o) { origin_ = o; }

  INLINE const GbVec3<T>& getDirection() const { return direction_; }
  INLINE void setDirection(const GbVec3<T>&d) { direction_ = d; }

  friend std::ostream& operator<<(std::ostream& s, const GbSegment<T>& v);

private:
  GbVec3<T> origin_;  // P
  GbVec3<T> direction_;  // D
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbSegment<T>& v)
{
  s<<typeid(v).name()<<": start("<<v.getOrigin()[0]<<","<<v.getOrigin()[1]<<","<<v.getOrigin()[2]<<") ";
  GbVec3<T> endp = v.getOrigin() + v.getDirection();
  s<<"end("<<endp[0]<<","<<endp[1]<<","<<endp[2]<<")"<<std::endl;
  return s;
}

#ifdef __GNUG__

//#include "GbSegment.in"
//#include "GbSegment.T"

#else

#pragma instantiate GbSegment<float>
#pragma instantiate GbSegment<double>
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbSegment<float>& v)
#pragma instantiate std::ostream& operator<<(std::ostream& s, const GbSegment<double>& v)

//#ifndef OUTLINE
//#include "GbSegment.in"
//#endif

#endif  // g++

#endif  // GBSEGMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/01/02 14:55:56  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
