#ifndef GOPRISMELEMENT_HH
#define GOPRISMELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoPrismElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoPrismElement
 */
class GoPrismElement 
  : public GoDefaultElement
{
public: 
  GoPrismElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    { 
      prism_ = new GoPrismElementBase<float>; 
    }
  GoPrismElement(GoDefaultElementBase<float> *t) 
    : GoDefaultElement(t) 
    {
      prism_ = new GoPrismElementBase<float>; 
    }

  virtual ~GoPrismElement() { delete prism_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoPrismElementBase<float>,float> (prism_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoPrismElementBase<float>,float> (prism_,i);
  }
  virtual int getNumVertices() const { return 6; }

  // The face normal
  virtual void computeNormal() {
  }

  // The layout of the object
  virtual int getEdgeList(int **l) const { 
    static int table_prism[18] = {0,1,0,3,0,5,1,5,1,2,5,4,3,4,3,2,2,4};
    *l = table_prism;
    return 18;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoPrismElementBase<float>,float> (prism_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoPrismElementBase<float>,float> (prism_,i);
  }
  virtual int getNumEdges() const { return 9; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoPrismElementBase<float>,float> (prism_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoPrismElementBase<float>,float> (prism_,i);
  }
  virtual int getNumNeighbours() const { return 5; }
  
  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoPrismElementBase<float>,float> (prism_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoPrismElementBase<float>,float> (prism_,i,face);
  }
  virtual int getNumChildren() const { return 2; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoPrismElementBase<float>,float> (prism_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoPrismElementBase<float>,float> (prism_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoPrismElement&);

protected:
  GoPrismElementBase<float> *prism_;
};

#endif // GOPRISMELEMENT_HH

/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.8  2002/03/18 09:58:55  prkipfer
| refactored element structure
|
|
+---------------------------------------------------------------------*/
