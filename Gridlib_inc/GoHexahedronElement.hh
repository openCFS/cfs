#ifndef GOHEXAHEDRONELEMENT_HH
#define GOHEXAHEDRONELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoHexahedronElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoHexahedronElement
 */
class GoHexahedronElement 
  : public GoDefaultElement
{
public: 
  GoHexahedronElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    {
      hexa_ = new GoHexahedronElementBase<float>; 
    }
  GoHexahedronElement(GoDefaultElementBase<float> *t) 
    : GoDefaultElement(t) 
    {
      hexa_ = new GoHexahedronElementBase<float>; 
    }

  virtual ~GoHexahedronElement() { delete hexa_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoHexahedronElementBase<float>,float> (hexa_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoHexahedronElementBase<float>,float> (hexa_,i);
  }
  virtual int getNumVertices() const { return 8; }

  // The face normal
  virtual void computeNormal() {
  }

  // The layout of the object
  virtual int getEdgeList(int **l) const { 
    static int table_hexa[24] = {0,1,0,3,0,4,1,5,4,5,1,2,5,6,4,7,3,2,7,6,3,7,2,6};
    *l = table_hexa;
    return 24;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoHexahedronElementBase<float>,float> (hexa_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoHexahedronElementBase<float>,float> (hexa_,i);
  }
  virtual int getNumEdges() const { return 12; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoHexahedronElementBase<float>,float> (hexa_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoHexahedronElementBase<float>,float> (hexa_,i);
  }
  virtual int getNumNeighbours() const { return 6; }
  
  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoHexahedronElementBase<float>,float> (hexa_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoHexahedronElementBase<float>,float> (hexa_,i,face);
  }
  virtual int getNumChildren() const { return 8; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoHexahedronElementBase<float>,float> (hexa_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoHexahedronElementBase<float>,float> (hexa_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoHexahedronElement&);

protected:
  GoHexahedronElementBase<float> *hexa_;
};

#endif // GOHEXAHEDRONELEMENT_HH

/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.8  2002/03/18 09:58:55  prkipfer
| refactored element structure
|
|
+---------------------------------------------------------------------*/
