#ifndef GOPRISMELEMENT_HH
#define GOPRISMELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoPrismElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoPrismElement
 */
class GoPrismElement 
  : public GoGeometryElement<float>
{
public: 
  GoPrismElement() { this_ = new GoPrismElementBase<float>; }
  GoPrismElement(GoPrismElementBase<float> *t) : this_(t) {}

  virtual ~GoPrismElement() { delete this_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoPrismElementBase<float>,float> (this_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoPrismElementBase<float>,float> (this_,i);
  }
  virtual int findVertex(GoVertex<float> *v) const {
    return ::findVertex<GoPrismElementBase<float>,float> (this_,v);
  }
  virtual GoVertex<float> *otherVertex(GoGeometryElement<float> *f) const {
    return ::otherVertex<GoPrismElementBase<float>,float> (this_,f);
  }
  virtual int getNumVertices() const {
    return 6;
  }

  // Other faces related to this object
  virtual GoGeometryElement<float> *otherFace(GoVertex<float> *v) const {
    return ::otherFace<GoPrismElementBase<float>,float> (this_,v);
  }

  // The face normal
  virtual void computeNormal() {
    ::computeNormal<GoPrismElementBase<float>,float> (this_);
  }
  virtual GbVec3<float> getNormal() const {
    return ::getNormal<GoPrismElementBase<float>,float> (this_);
  }
  virtual void setNormal(const GbVec3<float>& n) {
    ::setNormal<GoPrismElementBase<float>,float> (this_,n);
  }

  // The layout of the object
  virtual GbVec3<float> getOrigin() const {
//    return ::getOrigin<GoPrismElementBase<float>,float> (this_);
    return GbVec3<float>::ZERO;
  }
  virtual GbVec3<float> getEdge(int i) const {
//    return ::getEdge<GoPrismElementBase<float>,float> (this_,i);
    return GbVec3<float>::ZERO;
  }

  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbGeoStatusFlag f) {
    ::setFlag<GoPrismElementBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbGeoStatusFlag f) {
    ::delFlag<GoPrismElementBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbGeoStatusFlag f) const {
    return ::testFlag<GoPrismElementBase<float>,float> (this_,f);
  }

  // The object can be split into sub-objects
//  virtual void setNumSplits(int s);
//  virtual int getNumSplits() const;

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoPrismElementBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoPrismElementBase<float>,float> (this_);
  }

  // modification operations on the object
  virtual GbBool subdivide(const GbOracle &op) {
    if (op())
      return ::subdivide<GoPrismElementBase<float>,float> (this_);
    else
      return false;
  }

  // The scene part this object belongs to in the partition
  virtual void setPartition(int i) {
    ::setPartition<GoPrismElementBase<float>,float> (this_,i);
  }
  virtual int getPartition() const {
    return ::getPartition<GoPrismElementBase<float>,float> (this_);
  }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoPrismElementBase<float>,float> (this_,i,face);
  }
  virtual void setNeighbour(GoGeometryElement<float> *face) {
    ::setNeighbour<GoPrismElementBase<float>,float> (this_,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoPrismElementBase<float>,float> (this_,i);
  }
  virtual int findNeighbour(GoGeometryElement<float> *face) {
    return ::findNeighbour<GoPrismElementBase<float>,float> (this_,face);
  }
  
  // Parents and children objects
  virtual GoGeometryElement<float> *getChild(int i) const {
    return NULL; //::getChild<GoTriangleElementBase<float>,float> (this_,i);
  }
  virtual void setChild(int i, GoGeometryElement<float> *face) {
    //::setChild<GoQuadElementBase<float>,float> (this_,i,face);
  }
  virtual GoGeometryElement<float> *getParent() const {
    return NULL; //::getParent<GoTetrahedronElementBase<float>,float> (this_);
  }
  virtual void setParent(GoGeometryElement<float> *face) {
    //::setParent<GoTertahedronElementBase<float>,float> (this_,face);
  }

  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoPrismElementBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoPrismElementBase<float>,float> (this_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoPrismElement&);

private:
  GoPrismElementBase<float> *this_;
};

#endif // GOPRISMELEMENT_HH

