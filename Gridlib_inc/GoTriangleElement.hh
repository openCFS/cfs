/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTRIANGLEELEMENT_HH
#define GOTRIANGLEELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoDefaultElement.hh"
#include "GoVertex.hh"
#include "GoEdge.hh"
#include "GoTriangleElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoTriangleElement
 */
class GoTriangleElement 
  : public GoDefaultElement
{
public: 
  GoTriangleElement(int i=-1, int p=0, int l=0, const GbVec3<float>* n=NULL)
    : GoDefaultElement(i,p,l,n)
    { 
      tri_ = new GoTriangleElementBase<float>; 
    }
  GoTriangleElement(GoDefaultElementBase<float> *t) 
    : GoDefaultElement(t) 
    {
      tri_ = new GoTriangleElementBase<float>; 
    }

  virtual ~GoTriangleElement() { delete tri_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoTriangleElementBase<float>,float> (tri_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoTriangleElementBase<float>,float> (tri_,i);
  }
  virtual int getNumVertices() const { return 3; }

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
    static int table_tri[6] = {0,1,0,2,1,2};
    *l = table_tri;
    return 6;
  }
  virtual void setEdge(int i, GoEdge<float> *e) {
    ::setEdge<GoTriangleElementBase<float>,float> (tri_, i, e);
  }
  virtual GoEdge<float> *getEdge(int i) const {
    return ::getEdge<GoTriangleElementBase<float>,float> (tri_,i);
  }
  virtual int getNumEdges() const { return 3; }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoTriangleElementBase<float>,float> (tri_,i,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoTriangleElementBase<float>,float> (tri_,i);
  }
  virtual int getNumNeighbours() const { return 3; }
  
  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoTriangleElementBase<float>,float> (tri_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoTriangleElementBase<float>,float> (tri_,i,face);
  }
  virtual int getNumChildren() const { return 4; }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultElementBase<float>,float> (this_);
    ::statistic<GoTriangleElementBase<float>,float> (tri_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultElementBase<float>,float> (this_);
    ::shrink<GoTriangleElementBase<float>,float> (tri_);
  }


//  friend std::ostream& operator<<(std::ostream&, const GoTriangleElement&);

protected:
  GoTriangleElementBase<float> *tri_;
};

#endif // GOTRIANGLEELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.14  2002/03/18 09:58:56  prkipfer
| refactored element structure
|
| Revision 1.13  2001/12/18 13:25:17  prkipfer
| changed edge enumeration to be consistent throughout the lib
|
| Revision 1.12  2001/09/12 09:28:43  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.11  2001/01/02 15:21:35  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.10  2000/11/07 17:27:27  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.9  2000/11/07 09:58:53  uflabsik
| docu added
|
| Revision 1.8  2000/10/17 14:33:03  uflabsik
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
