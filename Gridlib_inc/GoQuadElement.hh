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
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoQuadElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoQuadElement
 */
class GoQuadElement 
  : public GoGeometryElement<float>
{
public: 
  GoQuadElement() { this_ = new GoQuadElementBase<float>; }
  GoQuadElement(GoQuadElementBase<float> *t) : this_(t) {}

  virtual ~GoQuadElement() { delete this_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoQuadElementBase<float>,float> (this_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoQuadElementBase<float>,float> (this_,i);
  }
  virtual int findVertex(GoVertex<float> *v) const {
    return ::findVertex<GoQuadElementBase<float>,float> (this_,v);
  }
  virtual GoVertex<float> *otherVertex(GoGeometryElement<float> *f) const {
    return ::otherVertex<GoQuadElementBase<float>,float> (this_,f);
  }
  virtual int getNumVertices() const {
    return 4;
  }

  // Other faces related to this object
  virtual GoGeometryElement<float> *otherFace(GoVertex<float> *v) const {
    return ::otherFace<GoQuadElementBase<float>,float> (this_,v);
  }

  // The face normal
  virtual void computeNormal() {
    ::computeNormal<GoQuadElementBase<float>,float> (this_);
  }
  virtual GbVec3<float> getNormal() const {
    return ::getNormal<GoQuadElementBase<float>,float> (this_);
  }
  virtual void setNormal(const GbVec3<float>& n) {
    ::setNormal<GoQuadElementBase<float>,float> (this_,n);
  }

  // The layout of the object
  virtual GbVec3<float> getOrigin() const {
//    return ::getOrigin<GoQuadElementBase<float>,float> (this_);
    return GbVec3<float>::ZERO;
  }
  virtual GbVec3<float> getEdge(int i) const {
//    return ::getEdge<GoQuadElementBase<float>,float> (this_,i);
    return GbVec3<float>::ZERO;
  }

  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbGeoStatusFlag f) {
    ::setFlag<GoQuadElementBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbGeoStatusFlag f) {
    ::delFlag<GoQuadElementBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbGeoStatusFlag f) const {
    return ::testFlag<GoQuadElementBase<float>,float> (this_,f);
  }

  // The object can be split into sub-objects
//  virtual void setNumSplits(int s);
//  virtual int getNumSplits() const;

  // modification operations on the object
  virtual GbBool subdivide(const GbOracle &op) {
    if (op())
      return ::subdivide<GoQuadElementBase<float>,float> (this_);
    else
      return false;
  }

  // The scene part this object belongs to in the partition
  virtual void setPartition(int i) {
    ::setPartition<GoQuadElementBase<float>,float> (this_,i);
  }
  virtual int getPartition() const {
    return ::getPartition<GoQuadElementBase<float>,float> (this_);
  }

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoQuadElementBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoQuadElementBase<float>,float> (this_);
  }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoQuadElementBase<float>,float> (this_,i,face);
  }
  virtual void setNeighbour(GoGeometryElement<float> *face) {
    ::setNeighbour<GoQuadElementBase<float>,float> (this_,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoQuadElementBase<float>,float> (this_,i);
  }
  
  virtual int findNeighbour(GoGeometryElement<float> *face) {
    return ::findNeighbour<GoQuadElementBase<float>,float> (this_,face);
  }

  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return ::getChild<GoQuadElementBase<float>,float> (this_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    ::setChild<GoQuadElementBase<float>,float> (this_,i,face);
  }
  virtual GoGeometryElement<float> *getParent() const {
    return ::getParent<GoQuadElementBase<float>,float> (this_);
  }
  virtual void setParent(GoGeometryElement<float> *face) {
    ::setParent<GoQuadElementBase<float>,float> (this_,face);
  }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoQuadElementBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoQuadElementBase<float>,float> (this_);
  }


//  friend std::ostream& operator<<(std::ostream&, const GoQuadElement&);

private:
  GoQuadElementBase<float> *this_;
};

#endif // GOQUADELEMENT_HH

