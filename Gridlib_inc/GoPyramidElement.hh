#ifndef GOPYRAMIDELEMENT_HH
#define GOPYRAMIDELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoPyramidElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoPyramidElement
 */
class GoPyramidElement 
  : public GoGeometryElement<float>
{
public: 
  GoPyramidElement() { this_ = new GoPyramidElementBase<float>; }
  GoPyramidElement(GoPyramidElementBase<float> *t) : this_(t) {}

  virtual ~GoPyramidElement() { delete this_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoPyramidElementBase<float>,float> (this_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoPyramidElementBase<float>,float> (this_,i);
  }
  virtual int findVertex(GoVertex<float> *v) const {
    return ::findVertex<GoPyramidElementBase<float>,float> (this_,v);
  }
  virtual GoVertex<float> *otherVertex(GoGeometryElement<float> *f) const {
    return ::otherVertex<GoPyramidElementBase<float>,float> (this_,f);
  }
  virtual int getNumVertices() const {
    return 5;
  }

  // Other faces related to this object
  virtual GoGeometryElement<float> *otherFace(GoVertex<float> *v) const {
    return ::otherFace<GoPyramidElementBase<float>,float> (this_,v);
  }

  // The face normal
  virtual void computeNormal() {
    ::computeNormal<GoPyramidElementBase<float>,float> (this_);
  }
  virtual GbVec3<float> getNormal() const {
    return ::getNormal<GoPyramidElementBase<float>,float> (this_);
  }
  virtual void setNormal(const GbVec3<float>& n) {
    ::setNormal<GoPyramidElementBase<float>,float> (this_,n);
  }

  // The layout of the object
  virtual GbVec3<float> getOrigin() const {
//    return ::getOrigin<GoPyramidElementBase<float>,float> (this_);
    return GbVec3<float>::ZERO;
  }
  virtual GbVec3<float> getEdge(int i) const {
//    return ::getEdge<GoPyramidElementBase<float>,float> (this_,i);
    return GbVec3<float>::ZERO;
  }

  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbGeoStatusFlag f) {
    ::setFlag<GoPyramidElementBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbGeoStatusFlag f) {
    ::delFlag<GoPyramidElementBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbGeoStatusFlag f) const {
    return ::testFlag<GoPyramidElementBase<float>,float> (this_,f);
  }

  // The object can be split into sub-objects
//  virtual void setNumSplits(int s);
//  virtual int getNumSplits() const;

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoPyramidElementBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoPyramidElementBase<float>,float> (this_);
  }

  // modification operations on the object
  virtual GbBool subdivide(const GbOracle &op) {
    if (op())
      return ::subdivide<GoPyramidElementBase<float>,float> (this_);
    else
      return false;
  }

  // The scene part this object belongs to in the partition
  virtual void setPartition(int i) {
    ::setPartition<GoPyramidElementBase<float>,float> (this_,i);
  }
  virtual int getPartition() const {
    return ::getPartition<GoPyramidElementBase<float>,float> (this_);
  }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoPyramidElementBase<float>,float> (this_,i,face);
  }
  virtual void setNeighbour(GoGeometryElement<float> *face) {
    ::setNeighbour<GoPyramidElementBase<float>,float> (this_,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoPyramidElementBase<float>,float> (this_,i);
  }
  virtual int findNeighbour(GoGeometryElement<float> *face) {
    return ::findNeighbour<GoPyramidElementBase<float>,float> (this_,face);
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
    ::statistic<GoPyramidElementBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoPyramidElementBase<float>,float> (this_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoPyramidElement&);

private:
  GoPyramidElementBase<float> *this_;
};

#endif // GOPYRAMIDELEMENT_HH

