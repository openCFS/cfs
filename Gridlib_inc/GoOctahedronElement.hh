#ifndef GOOCTAHEDRONELEMENT_HH
#define GOOCTAHEDRONELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoOctahedronElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoOctahedronElement
 */
class GoOctahedronElement 
  : public GoDefaultElement
{
public: 
  GoOctahedronElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    { 
      octa_ = new GoOctahedronElementBase<float>; 
    }
  GoOctahedronElement(GoDefaultElementBase<float> *t) 
    : GoDefaultElement(t) 
    {
      octa_ = new GoOctahedronElementBase<float>; 
    }

  virtual ~GoOctahedronElement() { delete octa_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoOctahedronElementBase<float>,float> (octa_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoOctahedronElementBase<float>,float> (octa_,i);
  }
  virtual int getNumVertices() const { return 6; }

  // The face normal
  virtual void computeNormal() {
  }

  // The layout of the object
  virtual int getEdgeList(int **l) const { 
    static int table_octa[24] = {0,2,0,3,0,4,0,5,1,2,1,3,1,4,1,5,3,2,4,3,5,4,2,5};
    *l = table_octa;
    return 24;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoOctahedronElementBase<float>,float> (octa_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoOctahedronElementBase<float>,float> (octa_,i);
  }
  virtual int getNumEdges() const { return 12; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoOctahedronElementBase<float>,float> (octa_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoOctahedronElementBase<float>,float> (octa_,i);
  }
  virtual int getNumNeighbours() const { return 8; }
  
  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoOctahedronElementBase<float>,float> (octa_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoOctahedronElementBase<float>,float> (octa_,i,face);
  }
  virtual int getNumChildren() const { return 10; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoOctahedronElementBase<float>,float> (octa_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoOctahedronElementBase<float>,float> (octa_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoOctahedronElement&);

protected:
  GoOctahedronElementBase<float> *octa_;
};

#endif // GOOCTAHEDRONELEMENT_HH

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
