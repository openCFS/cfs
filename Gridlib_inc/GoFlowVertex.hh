/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOFLOWVERTEX_HH
#define GOFLOWVERTEX_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoDefaultVertexBase.hh"
#include "GoFlowVertexBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoFlowVertex

  Vertex class for fluid dynamics meshes
 */
class GoFlowVertex
  : public GoVertex<float>
{
public: 
  GoFlowVertex(int i = -1, float x = 0, float y = 0, float z = 0, 
	       float d = 0, float mx = 0, float my = 0, float mz = 0, float e = 0) {
    this_ = new GoDefaultVertexBase<float>(i, x, y, z);
    value_ = new GoFlowVertexBase<float>(d,mx,my,mz,e);
  }
  GoFlowVertex(GoDefaultVertexBase<float> *t) : this_(t) {
    value_ = new GoFlowVertexBase<float>;
  }

  virtual ~GoFlowVertex() { delete this_; delete value_; }

  //! ... want some more ?
  virtual GoVertex<float> *clone() const {
    return new GoFlowVertex;
  }

  //! The position in 3D space
  virtual void setPosition(float x, float y, float z) {
    ::setPosition<GoDefaultVertexBase<float>,float> (this_,x,y,z);
  }
  virtual void setPosition(const GbVec3<float>& v) {
    ::setPosition<GoDefaultVertexBase<float>,float> (this_,v);
  }
  virtual void getPosition(float& xx, float& yy, float& zz) const {
    ::getPosition<GoDefaultVertexBase<float>,float> (this_,xx,yy,zz);
  }
  virtual GbVec3<float> getPosition() const {
    return ::getPosition<GoDefaultVertexBase<float>,float> (this_);
  }

  //! A normal at this vertex in 3D space
  virtual void computeNormal(std::vector<GoGeometryElement<float> *>& S) {
    ::computeNormal<GoDefaultVertexBase<float>,float> (this_,S);
  }
  virtual void setNormal(float x, float y, float z) {
    ::setNormal<GoDefaultVertexBase<float>,float> (this_,x,y,z);
  }
  virtual void setNormal(const GbVec3<float>& v) {
    ::setNormal<GoDefaultVertexBase<float>,float> (this_,v);
  }
  virtual void getNormal(float& xx, float& yy, float& zz) const {
    ::getNormal<GoDefaultVertexBase<float>,float> (this_,xx,yy,zz);
  }
  virtual GbVec3<float> getNormal() const {
    return ::getNormal<GoDefaultVertexBase<float>,float> (this_);
  }

  //! Integer to identify the vertex
  //! Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoDefaultVertexBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoDefaultVertexBase<float>,float> (this_);
  }

  //! Status flags of this vertex indicating semantic defined by the
  //! object using this class
  virtual void setFlag(GbVertexFlag f) {
    ::setFlag<GoDefaultVertexBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbVertexFlag f) {
    ::delFlag<GoDefaultVertexBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbVertexFlag f) const {
    return ::testFlag<GoDefaultVertexBase<float>,float> (this_,f);
  }

  //! The face object, this vertex belongs to
  virtual void setElement(GoGeometryElement<float> *f) {
    ::setElement<GoDefaultVertexBase<float>,float> (this_,f);
  }
  virtual GoGeometryElement<float> *getElement() const {
    return ::getElement<GoDefaultVertexBase<float>,float> (this_);
  }

  //! The neighbor valence of the vertex
  virtual void setValence(int v) {
    ::setValence<GoDefaultVertexBase<float>,float> (this_,v);
  }
  virtual int getValence() const {
    return ::getValence<GoDefaultVertexBase<float>,float> (this_);
  }

  //! Print memory pool statistics
  virtual void statistic() const {
    ::statistic<GoDefaultVertexBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultVertexBase<float>,float> (this_);
  }

  //! Handle the scalar values
  virtual void setDensity(float v) {
    ::setDensity<GoFlowVertexBase<float>,float> (value_,v);
  }
  virtual float getDensity() const {
    return ::getDensity<GoFlowVertexBase<float>,float> (value_);
  }
  virtual void setEnergy(float v) {
    ::setEnergy<GoFlowVertexBase<float>,float> (value_,v);
  }
  virtual float getEnergy() const {
    return ::getEnergy<GoFlowVertexBase<float>,float> (value_);
  }

  //! Handle momentum vector
  virtual void setMomentum(const GbVec3<float>& v) {
    ::setMomentum<GoFlowVertexBase<float>,float> (value_,v);
  }
  virtual GbVec3<float> getMomentum() const {
    return ::getMomentum<GoFlowVertexBase<float>,float> (value_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoVertex<float>&);

protected:
  GoDefaultVertexBase<float> *this_;
  GoFlowVertexBase<float> *value_;
};

#endif // GOFLOWVERTEX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/02/13 11:03:17  prkipfer
| introduced boundary vertices
|
| Revision 1.2  2001/01/02 15:21:34  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.1  2000/11/07 17:27:26  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
