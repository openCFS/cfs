/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTRIANGLEELEMENTBASE_HH
#define GOTRIANGLEELEMENTBASE_HH

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
class GoTriangleElementBase 
{
public: 

  //! Constructor: all values initialized to zero, if not specified
  GoTriangleElementBase();
  ~GoTriangleElementBase();

  //! Operations on vertices of the triangle
  INLINE void setVertex(int i, GoVertex<T> *v);
  INLINE GoVertex<T> *getVertex(int i) const;

  //! Operations on edges of the triangle
  INLINE void setEdge(int i, GoEdge<T> *e);
  INLINE GoEdge<T> *getEdge(int i) const;

  //! The neighboring triangles
  INLINE void setNeighbour(int i, GoGeometryElement<T> *face);
  INLINE GoGeometryElement<T> *getNeighbour(int i) const;

  //! The parent and the child triangles
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
  friend std::ostream& operator<<(std::ostream&, const GoTriangleElementBase<T>&);

private:
  //! Inhibit assignment to a triangle
  GoTriangleElementBase<T>& operator=(const GoTriangleElementBase<T>& rhs);
  
  //! This is the memory layout of this class within the memory pool
  typedef struct {
    GoVertex<T>* vertices[3];
    GoEdge<T>* edges[3];
    GoGeometryElement<T>* neighbours[3];
    GoGeometryElement<T>* children[4];
  } TriangleMem;

  //! The memory pool and the link to the allocated memory space
  static GbMemPool<TriangleMem> pool_;
  TriangleMem data_;
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GoTriangleElementBase<T>& v)
{
//  s<<typeid(v).name()<<" id: "<<v.getId()<<std::endl;
//  s<<"\tn("<<v.getNormal()[0]<<","<<v.getNormal()[1]<<","<<v.getNormal()[2]<<") "<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GoTriangleElementBase.in"
#include "GoTriangleElementBase.T"

#else

#pragma instantiate GoTriangleElementBase<float>
#pragma instantiate GoTriangleElementBase<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoTriangleElementBase<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GoTriangleElementBase<double>&)
#pragma instantiate GbMemPool<GoTriangleElementBase<float>::TriangleMem>
#pragma instantiate GbMemPool<GoTriangleElementBase<double>::TriangleMem>

#ifndef OUTLINE
#include "GoTriangleElementBase.in"
#endif

#endif // g++

#endif // GOTRIANGLEELEMENTBASE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.12  2002/03/18 10:00:27  prkipfer
| refactored element structure
|
| Revision 1.11  2001/09/12 09:20:05  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.10  2001/01/02 15:16:48  prkipfer
| changed to use new classes and GbMath
|
| Revision 1.9  2000/11/07 17:27:39  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.8  2000/10/17 14:34:46  uflabsik
| new primitives added
|
| Revision 1.7  2000/09/03 14:05:43  uflabsik
| functionality added
|
| Revision 1.6  2000/08/28 09:25:22  prkipfer
| added localized subdivision mechanism
|
| Revision 1.5  2000/07/21 10:40:45  prkipfer
| dropped g++ support
|
| Revision 1.4  2000/07/20 13:32:38  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.3  2000/06/16 13:33:50  prkipfer
| further hierarchy fixes
|
| Revision 1.2  2000/06/14 15:39:14  prkipfer
| improved base classes and added funcstruct processing for mesh
|
| Revision 1.1.1.1  2000/06/08 16:24:44  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
