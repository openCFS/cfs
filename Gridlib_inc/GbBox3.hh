/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBBOX3_HH
#define  GBBOX3_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GbVec3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbBox3
{
public:
  GbBox3 ();

  INLINE GbVec3<T>& getCenter ();
  INLINE const GbVec3<T>& getCenter () const;
  INLINE void setCenter(const GbVec3<T>& c);

  INLINE GbVec3<T>& getAxis (int i);
  INLINE const GbVec3<T>& getAxis (int i) const;
  INLINE void setAxis(int i, const GbVec3<T>&a);
  INLINE GbVec3<T>* getAxes ();
  INLINE const GbVec3<T>* getAxes () const;
  INLINE void setAxes(const GbVec3<T>* a);

  INLINE T& getExtent (int i);
  INLINE const T& getExtent (int i) const;
  INLINE void setExtent(int i, const T& e);
  INLINE T* getExtents ();
  INLINE const T* getExtents () const;

  void computeVertices (GbVec3<T> akVertex[8]) const;

  //! This operator displays the values
  friend std::ostream& operator<<(std::ostream&, const GbBox3<T>&);

private:
  GbVec3<T> center_;
  GbVec3<T> axis_[3];
  T extent_[3];
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbBox3<T>& v)
{
  s<<typeid(v).name()<<": center("<<v.getCenter()[0]<<","<<v.getCenter()[1]<<","<<v.getCenter()[2]<<")"<<std::endl;
  s<<"\taxis 0: ("<<v.getAxis(0)[0]<<","<<v.getAxis(0)[1]<<","<<v.getAxis(0)[2]<<") = "<<v.getExtent(0)<<std::endl;
  s<<"\taxis 1: ("<<v.getAxis(1)[0]<<","<<v.getAxis(1)[1]<<","<<v.getAxis(1)[2]<<") = "<<v.getExtent(1)<<std::endl;
  s<<"\taxis 2: ("<<v.getAxis(2)[0]<<","<<v.getAxis(2)[1]<<","<<v.getAxis(2)[2]<<") = "<<v.getExtent(2)<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GbBox3.in"
#include "GbBox3.T"

#else

#pragma instantiate GbBox3<float>
#pragma instantiate GbBox3<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbBox3<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbBox3<double>&)

#ifndef OUTLINE
#include "GbBox3.in"
#endif

#endif  // g++

#endif  // GBBOX3_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/01/02 14:59:02  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
