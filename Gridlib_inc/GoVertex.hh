/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOVERTEX_HH
#define GOVERTEX_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>

#include "GbTypes.hh"
#include "GbVec3.hh"

template <class T> class GoGeometryElement;
template <class T> class GoEdge;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoVertex

  Abstract base class for vertices
 */
template <class T>
class GoVertex
{
private:
  //! CS
  GbBool useElementsCounter;
  int elements;

public: 

  //! constructor zum probieren - CS
  GoVertex()
  {
    useElementsCounter = false;
    elements = 0;
  }

  //! CS
  void addElement() 
  {
    useElementsCounter = true;
    elements++;

    if (elements > 0)
      delFlag(DELETED_V);
  }

  //! CS
  void removeElement()
  {
    if (!useElementsCounter)
      return;

    elements--;
    if (elements <= 0)
      setFlag(DELETED_V);
  }
  //! CS
  int getNumElements()
  {return elements;}

  //! Cannot be instantiated
  virtual ~GoVertex() {}

  //! ... want some more ?
  virtual GoVertex<T> *clone() const = 0;

  //! The position in 3D space
  virtual void setPosition(T x, T y, T z) = 0;
  virtual void setPosition(const GbVec3<T>& v) = 0;
  virtual void getPosition(T& xx, T& yy, T& zz) const = 0;
  virtual GbVec3<T> getPosition() const = 0;

  //! A normal at this vertex in 3D space
  virtual void computeNormal(std::vector<GoGeometryElement<T> *>& S) = 0;
  virtual void setNormal(T x, T y, T z) = 0;
  virtual void setNormal(const GbVec3<T>& v) = 0;
  virtual void getNormal(T& xx, T& yy, T& zz) const = 0;
  virtual GbVec3<T> getNormal() const = 0;

  //! Integer to identify the vertex
  //! Has no meaning to this class's implementation
  virtual void setId(int i) = 0;
  virtual int getId() const = 0;

  //! Status flags of this vertex indicating semantic defined by the
  //! object using this class
  virtual void setFlag(GbVertexFlag f) = 0;
  virtual void delFlag(GbVertexFlag f) = 0;
  virtual GbBool testFlag(GbVertexFlag f) const = 0;

  //! The face object, this vertex belongs to
  virtual void setElement(GoGeometryElement<T> *f) = 0;
  virtual GoGeometryElement<T> *getElement() const = 0;

  //! The neighbor valence of the vertex
  virtual void setValence(int v) = 0;
  virtual int getValence() const = 0;

  //! Print memory pool statistics
  virtual void statistic() const = 0;
  virtual void shrink() = 0;

//  friend std::ostream& operator<<(std::ostream&, const GoVertex<T>&);

};

#include "GbExternalAdapters.hh"

#endif // GOVERTEX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.3  2001/09/12 09:28:43  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.2  2001/01/02 15:21:35  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.1  2000/11/07 17:27:27  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
