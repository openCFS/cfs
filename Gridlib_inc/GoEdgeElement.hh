/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOEDGEELEMENT_HH
#define GOEDGEELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"

template <class T> class GoVertex;
template <class T> class GoGeometryElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


/*!
  class GoEdgeElement

  Abstract base class for edges
 */
template <class T>
class GoEdgeElement 
{
public: 
  //! Cannot be instantiated
  virtual ~GoEdgeElement() {}

  //! ... want some more ?
  virtual GoEdgeElement<T> *clone() const = 0;

  //! Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<T> *v) = 0;
  virtual GoVertex<T> *getVertex(int i) const = 0;
  virtual int findVertex(GoVertex<T> *v) const = 0;

  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbEdgeStatusFlag f) = 0;
  virtual void delFlag(GbEdgeStatusFlag f) = 0;
  virtual GbBool testFlag(GbEdgeStatusFlag f) const = 0;

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) = 0;
  virtual int getId() const = 0;

  //! The face object, this edge belongs to
  virtual void setElement(GoGeometryElement<T> *f) = 0;
  virtual GoGeometryElement<T> *getElement() const = 0;

  //! The neighbor valence of the edge
  virtual void setValence(int v) = 0;
  virtual int getValence() const = 0;

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const = 0;
  virtual void shrink() = 0;

//  friend std::ostream& operator<<(std::ostream&, const GoEdgeElement<T>&);
};


#include "GbExternalAdapters.hh"

#endif // GOEDGEELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 15:20:17  prkipfer
| introduced edge elements
|
|
+---------------------------------------------------------------------*/
