#ifndef GOHEXAHEDRONELEMENT_HH
#define GOHEXAHEDRONELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoHexahedronElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoHexahedronElement
 */
class GoHexahedronElement 
  : public GoGeometryElement<float>
{
public: 
  GoHexahedronElement() { this_ = new GoHexahedronElementBase<float>; }
  GoHexahedronElement(GoHexahedronElementBase<float> *t) : this_(t) {}

  virtual ~GoHexahedronElement() { delete this_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoHexahedronElementBase<float>,float> (this_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoHexahedronElementBase<float>,float> (this_,i);
  }
  virtual int findVertex(GoVertex<float> *v) const {
    return ::findVertex<GoHexahedronElementBase<float>,float> (this_,v);
  }
  virtual GoVertex<float> *otherVertex(GoGeometryElement<float> *f) const {
    return ::otherVertex<GoHexahedronElementBase<float>,float> (this_,f);
  }
  virtual int getNumVertices() const {
    return 8;
  }

  // Other faces related to this object
  virtual GoGeometryElement<float> *otherFace(GoVertex<float> *v) const {
    return ::otherFace<GoHexahedronElementBase<float>,float> (this_,v);
  }

  // The face normal
  virtual void computeNormal() {
    ::computeNormal<GoHexahedronElementBase<float>,float> (this_);
  }
  virtual GbVec3<float> getNormal() const {
    return ::getNormal<GoHexahedronElementBase<float>,float> (this_);
  }
  virtual void setNormal(const GbVec3<float>& n) {
    ::setNormal<GoHexahedronElementBase<float>,float> (this_,n);
  }

  // The layout of the object
  virtual GbVec3<float> getOrigin() const {
//    return ::getOrigin<GoHexahedronElementBase<float>,float> (this_);
    return GbVec3<float>::ZERO;
  }
  virtual GbVec3<float> getEdge(int i) const {
//    return ::getEdge<GoHexahedronElementBase<float>,float> (this_,i);
    return GbVec3<float>::ZERO;
  }

  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbGeoStatusFlag f) {
    ::setFlag<GoHexahedronElementBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbGeoStatusFlag f) {
    ::delFlag<GoHexahedronElementBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbGeoStatusFlag f) const {
    return ::testFlag<GoHexahedronElementBase<float>,float> (this_,f);
  }

  // The object can be split into sub-objects
//  virtual void setNumSplits(int s);
//  virtual int getNumSplits() const;

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoHexahedronElementBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoHexahedronElementBase<float>,float> (this_);
  }

  // modification operations on the object
  virtual GbBool subdivide(const GbOracle &op) {
    if (op())
      return ::subdivide<GoHexahedronElementBase<float>,float> (this_);
    else
      return false;
  }

  // The scene part this object belongs to in the partition
  virtual void setPartition(int i) {
    ::setPartition<GoHexahedronElementBase<float>,float> (this_,i);
  }
  virtual int getPartition() const {
    return ::getPartition<GoHexahedronElementBase<float>,float> (this_);
  }

  // The neighboring face objects
  virtual void setNeighbour(int i, GoGeometryElement<float> *face) {
    ::setNeighbour<GoHexahedronElementBase<float>,float> (this_,i,face);
  }
  virtual void setNeighbour(GoGeometryElement<float> *face) {
    ::setNeighbour<GoHexahedronElementBase<float>,float> (this_,face);
  }
  virtual GoGeometryElement<float> *getNeighbour(int i) const {
    return ::getNeighbour<GoHexahedronElementBase<float>,float> (this_,i);
  }
  virtual int findNeighbour(GoGeometryElement<float> *face) {
    return ::findNeighbour<GoHexahedronElementBase<float>,float> (this_,face);
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
    ::statistic<GoHexahedronElementBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoHexahedronElementBase<float>,float> (this_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoHexahedronElement&);

private:
  GoHexahedronElementBase<float> *this_;
};

#endif // GOHEXAHEDRONELEMENT_HH

