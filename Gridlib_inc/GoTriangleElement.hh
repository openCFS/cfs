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
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoTriangleElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


/*!
  class GoTriangleElement
 */
class GoTriangleElement 
  : public GoGeometryElement<float>
{
public: 
  GoTriangleElement() { this_ = new GoTriangleElementBase<float>; }
  GoTriangleElement(GoTriangleElementBase<float> *t) : this_(t) {}

  virtual ~GoTriangleElement() { delete this_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoTriangleElementBase<float>,float> (this_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoTriangleElementBase<float>,float> (this_,i);
  }
  virtual int findVertex(GoVertex<float> *v) const {
    return ::findVertex<GoTriangleElementBase<float>,float> (this_,v);
  }
  virtual GoVertex<float> *otherVertex(GoGeometryElement<float> *f) const {
    return ::otherVertex<GoTriangleElementBase<float>,float> (this_,f);
  }
  virtual int getNumVertices() const {
    return 3;
  }

  // Other faces related to this object
  virtual GoGeometryElement<float> *otherFace(GoVertex<float> *v) const {
    return ::otherFace<GoTriangleElementBase<float>,float> (this_,v);
  }

  // The face normal
  virtual void computeNormal() {
    ::computeNormal<GoTriangleElementBase<float>,float> (this_);
  }
  virtual GbVec3<float> getNormal() const {
    return ::getNormal<GoTriangleElementBase<float>,float> (this_);
  }
  virtual void setNormal(const GbVec3<float>& n) {
    ::setNormal<GoTriangleElementBase<float>,float> (this_,n);
  }

  // The layout of the object
  virtual GbVec3<float> getOrigin() const {
    return ::getOrigin<GoTriangleElementBase<float>,float> (this_);
  }
  virtual GbVec3<float> getEdge(int i) const {
    return ::getEdge<GoTriangleElementBase<float>,float> (this_,i);
  }

  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbGeoStatusFlag f) {
    ::setFlag<GoTriangleElementBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbGeoStatusFlag f) {
    ::delFlag<GoTriangleElementBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbGeoStatusFlag f) const {
    return ::testFlag<GoTriangleElementBase<float>,float> (this_,f);
  }

  // The object can be split into sub-objects
//  virtual void setNumSplits(int s);
//  virtual int getNumSplits() const;

  // modification operations on the object
  virtual GbBool subdivide(const GbOracle &op) {
    if (op())
      return ::subdivide<GoTriangleElementBase<float>,float> (this_);
    else
      return false;
  }

  // The scene part this object belongs to in the partition
  virtual void setPartition(int i) {
    ::setPartition<GoTriangleElementBase<float>,float> (this_,i);
  }
  virtual int getPartition() const {
    return ::getPartition<GoTriangleElementBase<float>,float> (this_);
  }

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoTriangleElementBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoTriangleElementBase<float>,float> (this_);
  }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoTriangleElementBase<float>,float> (this_,i,face);
  }
  virtual void setNeighbour(GoGeometryElement<float> *face) {
    ::setNeighbour<GoTriangleElementBase<float>,float> (this_,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoTriangleElementBase<float>,float> (this_,i);
  }
  
  virtual int findNeighbour(GoGeometryElement<float> *face) {
    return ::findNeighbour<GoTriangleElementBase<float>,float> (this_,face);
  }

  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoTriangleElementBase<float>,float> (this_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoTriangleElementBase<float>,float> (this_,i,face);
  }
  virtual GoGeometryElement<float> *getParent() const {
    return ::getParent<GoTriangleElementBase<float>,float> (this_);
  }
  virtual void setParent(GoGeometryElement<float> *face) {
    ::setParent<GoTriangleElementBase<float>,float> (this_,face);
  }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoTriangleElementBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoTriangleElementBase<float>,float> (this_);
  }


//  friend std::ostream& operator<<(std::ostream&, const GoTriangleElement&);

private:
  GoTriangleElementBase<float> *this_;
};

#endif // GOTRIANGLEELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
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
