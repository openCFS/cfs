/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOGEOMETRYELEMENT_HH
#define GOGEOMETRYELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>

#include "GbTypes.hh"
#include "GbVec3.hh"

template <class T> class GoVertex;
template <class T> class GoEdge;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoGeometryElement

  Abstract base class for geometrical elements
 */
template <class T>
class GoGeometryElement
{
public: 

  //! Cannot be instantiated
  virtual ~GoGeometryElement() {}

  //! Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<T> *v) = 0;
  virtual GoVertex<T> *getVertex(int i) const = 0;
  virtual int findVertex(GoVertex<T> *v) const = 0;
  virtual GoVertex<T> *otherVertex(GoGeometryElement<T> *f) const = 0;
  virtual int getNumVertices() const = 0;

  //! Other faces related to this object
  virtual GoGeometryElement<T> *otherFace(GoVertex<T> *v) const = 0;

  //! The face normal
  virtual void computeNormal() = 0;
  virtual GbVec3<T> getNormal() const = 0;
  virtual void setNormal(const GbVec3<T>& n) = 0;

  //! The layout of the object
  virtual int getEdgeList(int **l) const = 0;
  virtual GbVec3<T> getOrigin() const = 0;
  virtual GbVec3<T> getEdgeVec(int i) const = 0;
  virtual void setEdge(int i, GoEdge<T> *e) = 0;
  virtual GoEdge<T> *getEdge(int i) const = 0;
  virtual GoEdge<T> *findEdge(GoVertex<T>* v0, GoVertex<T>* v1) const = 0;
  virtual int findEdge(GoEdge<T> *e) const = 0;
  virtual int getNumEdges() const = 0;

  //! Status flags indicating semantic defined by the
  //! object using this class
  virtual void setFlag(GbGeoStatusFlag f) = 0;
  virtual void delFlag(GbGeoStatusFlag f) = 0;
  virtual GbBool testFlag(GbGeoStatusFlag f) const = 0;

  //! The object can be split into sub-objects
//  virtual void setNumSplits(int s);
//  virtual int getNumSplits() const;

  //! Integer to identify the object
  //! Has no meaning to this class's implementation
  virtual void setId(int i) = 0;
  virtual int getId() const = 0;

  //! modification operations on the object
  virtual GbBool subdivide(const GbOracle &op) = 0;

  //! The scene part this object belongs to in the partition
  virtual void setPartition(int i) = 0;
  virtual int getPartition() const = 0;

  //! The neighboring face objects
  virtual int getNumNeighbours() const = 0;
  virtual void setNeighbour(int i, GoGeometryElement<T> *face) = 0;
  virtual GbBool trySetNeighbour(GoGeometryElement<T> *face) = 0;
  virtual GoGeometryElement<T> *getNeighbour(int i) const = 0;
  virtual int findNeighbour(GoGeometryElement<T> *face) = 0;

  //! Parents and children objects
  virtual GoGeometryElement<T> *getChild(int i) const = 0;
  virtual void setChild(int i, GoGeometryElement<T> *face) = 0;
  virtual int getNumChildren() const = 0;
  virtual void clearChildren() = 0;
  virtual GoGeometryElement<T> *getParent() const = 0;
  virtual void setParent(GoGeometryElement<T> *face) = 0;
  virtual void setLevel(int l) = 0;
  virtual int getLevel() const = 0;

  //! Print memory statistics and reduce memory space needed
  virtual void statistic() const = 0;
  virtual void shrink() = 0;

//  friend std::ostream& operator<<(std::ostream&, const GoGeometryElement<T>&);

};

#include "GbExternalAdapters.hh"

#endif // GOGEOMETRYELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.13  2002/03/18 09:58:55  prkipfer
| refactored element structure
|
| Revision 1.12  2001/09/12 09:28:42  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.11  2001/01/02 15:21:34  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.10  2000/11/07 17:27:26  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.9  2000/11/07 09:58:52  uflabsik
| docu added
|
| Revision 1.8  2000/10/17 14:33:00  uflabsik
| added new primitives
|
| Revision 1.7  2000/09/07 16:57:17  prkipfer
| added subdivision
|
| Revision 1.6  2000/09/03 14:06:14  uflabsik
| functionality added
|
| Revision 1.5  2000/08/28 09:25:21  prkipfer
| added localized subdivision mechanism
|
| Revision 1.4  2000/07/21 10:40:43  prkipfer
| dropped g++ support
|
| Revision 1.3  2000/06/16 13:33:43  prkipfer
| further hierarchy fixes
|
| Revision 1.2  2000/06/14 15:39:10  prkipfer
| improved base classes and added funcstruct processing for mesh
|
| Revision 1.1.1.1  2000/06/08 16:24:43  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
