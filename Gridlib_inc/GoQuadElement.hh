/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOQUADELEMENT_HH
#define GOQUADELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoQuadElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoQuadElement
 */
class GoQuadElement 
  : public GoDefaultElement
{
public: 
  GoQuadElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    { 
      quad_ = new GoQuadElementBase<float>; 
    }
  GoQuadElement(GoDefaultElementBase<float> *t) 
    : GoDefaultElement(t) 
    {
      quad_ = new GoQuadElementBase<float>; 
    }

  virtual ~GoQuadElement() { delete quad_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoQuadElementBase<float>,float> (quad_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoQuadElementBase<float>,float> (quad_,i);
  }
  virtual int getNumVertices() const { return 4; }

  // The face normal
  virtual void computeNormal() {
    float x0,y0,z0, x1,y1, z1, x2,y2,z2, xn,yn,zn;
    float n;

//  debugmsg("computing face normal");
    getVertex(0)->getPosition(x0,y0,z0);
    getVertex(1)->getPosition(x1,y1,z1);
    getVertex(2)->getPosition(x2,y2,z2);
  
    xn = ((y1-y0)*(z2-z0) - (y2-y0)*(z1-z0));
    yn = ((z1-z0)*(x2-x0) - (z2-z0)*(x1-x0));
    zn = ((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0));
  
    n = GbMath<float>::Sqrt(xn*xn + yn*yn + zn*zn);
  
    if (n != 0.0f)
      setNormal(GbVec3<float>(xn/n,yn/n,zn/n));
    else
      setNormal(GbVec3<float>::ZERO);
  }

  // The layout of the object
  virtual int getEdgeList(int **l) const { 
    static int table_quad[8] = {0,1,0,3,1,2,3,2};
    *l = table_quad;
    return 8;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoQuadElementBase<float>,float> (quad_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoQuadElementBase<float>,float> (quad_,i);
  }
  virtual int getNumEdges() const { return 4; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoQuadElementBase<float>,float> (quad_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoQuadElementBase<float>,float> (quad_,i);
  }
  virtual int getNumNeighbours() const { return 4; }

  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoQuadElementBase<float>,float> (quad_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoQuadElementBase<float>,float> (quad_,i,face);
  }
  virtual int getNumChildren() const { return 4; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoQuadElementBase<float>,float> (quad_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoQuadElementBase<float>,float> (quad_);
  }


//  friend std::ostream& operator<<(std::ostream&, const GoQuadElement&);

protected:
  GoQuadElementBase<float> *quad_;
};

#endif // GOQUADELEMENT_HH

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
