/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTETRAHEDRONELEMENTBASE_HH
#define GOTETRAHEDRONELEMENTBASE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GbVec3.hh"
#include "GbMemPool.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"

template <class T> class GoGeometryElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoTetrahedronElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoTetrahedronElementBase();
  ~GoTetrahedronElementBase();

  //! Operations on vertices of the tetrahedra
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;

  //! Operations on edges of the tetrahedra
  INLINE void setEdge(int i, GoEdge<T> *e);
  INLINE GoEdge<T> *getEdge(int i) const;

  //! The neighboring elements
  INLINE void setNeighbour(int i, GoGeometryElement<T> *face);
  INLINE GoGeometryElement<T> *getNeighbour(int i) const;

  INLINE GoGeometryElement<T> *getChild(int i) const;
  INLINE void setChild(int i, GoGeometryElement<T> *face);

  //! Print memory pool statistics
  void poolStatistic() const;
  void shrink();

  //! The memory handling: we get all memory from a pool
  //! Note that these operators are automatic-inline
  static void *operator new(size_t size) { return pool_.alloc(size); }
  static void operator delete(void *p, size_t size) { pool_.free(p, size); }

  //! This operator pretty prints the vertex values
  friend std::ostream& operator<<(std::ostream&, const GoTetrahedronElementBase<T>&);

private:
  //! Inhibit assignment to a tetrahedra
  GoTetrahedronElementBase<T>& operator=(const GoTetrahedronElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[4];
    GoEdge<T>* edges[6];
    GoGeometryElement<T>* neighbours[4];
    GoGeometryElement<T>* children[8];
  } TetraMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<TetraMem> pool_;
  TetraMem data_;

};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoTetrahedronElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
//    s<<"\tn("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") "<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoTetrahedronElementBase.in"
#include "GoTetrahedronElementBase.T"

#else

#pragma instantiate GoTetrahedronElementBase<float>
#pragma instantiate GoTetrahedronElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoTetrahedronElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoTetrahedronElementBase<double>&)
#pragma instantiate GbMemPool<GoTetrahedronElementBase<float>::TetraMem>
#pragma instantiate GbMemPool<GoTetrahedronElementBase<double>::TetraMem>

#ifndef OUTLINE
#include "GoTetrahedronElementBase.in"
#endif

#endif // g++

#endif // GOTETRAHEDRONELEMENTBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.9  2002/03/18 10:00:26  prkipfer
| refactored element structure
|
| Revision 1.8  2001/09/12 09:20:04  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.7  2001/01/02 15:16:47  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.6  2000/11/07 17:27:39  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.5  2000/09/28 14:15:13  uflabsik
| errors corrected
|
| Revision 1.4  2000/09/03 14:05:42  uflabsik
| functionality added
|
| Revision 1.3  2000/08/28 09:25:22  prkipfer
| added localized subdivision mechanism
|
| Revision 1.2  2000/07/21 10:40:45  prkipfer
| dropped g++ support
|
| Revision 1.1  2000/07/03 16:15:03  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
