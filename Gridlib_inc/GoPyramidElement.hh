#ifndef GOPYRAMIDELEMENT_HH
#define GOPYRAMIDELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoPyramidElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoPyramidElement
 */
class GoPyramidElement 
  : public GoDefaultElement
{
public: 
  GoPyramidElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    { 
      pyra_ = new GoPyramidElementBase<float>; 
    }
  GoPyramidElement(GoDefaultElementBase<float> *t) 
    : GoDefaultElement(t) 
    {
      pyra_ = new GoPyramidElementBase<float>; 
    }

  virtual ~GoPyramidElement() { delete pyra_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoPyramidElementBase<float>,float> (pyra_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoPyramidElementBase<float>,float> (pyra_,i);
  }
  virtual int getNumVertices() const { return 5; }

  // The face normal
  virtual void computeNormal() {
  }

  // The layout of the object
  virtual int getEdgeList(int **l) const { 
    static int table_pyra[16] = {0,1,0,2,0,3,0,4,2,1,3,2,4,3,1,4};
    *l = table_pyra;
    return 16;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoPyramidElementBase<float>,float> (pyra_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoPyramidElementBase<float>,float> (pyra_,i);
  }
  virtual int getNumEdges() const { return 8; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoPyramidElementBase<float>,float> (pyra_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoPyramidElementBase<float>,float> (pyra_,i);
  }
  virtual int getNumNeighbours() const { return 5; }
  
  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoPyramidElementBase<float>,float> (pyra_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoPyramidElementBase<float>,float> (pyra_,i,face);
  }
  virtual int getNumChildren() const { return 5; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoPyramidElementBase<float>,float> (pyra_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoPyramidElementBase<float>,float> (pyra_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoPyramidElement&);

protected:
  GoPyramidElementBase<float> *pyra_;
};

#endif // GOPYRAMIDELEMENT_HH

/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.8  2002/03/18 09:58:56  prkipfer
| refactored element structure
|
|
+---------------------------------------------------------------------*/
