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

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoDefaultVertex.hh"
#include "GoFlowVertexBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoFlowVertex

  Vertex class for fluid dynamics meshes
 */
class GoFlowVertex
  : public GoDefaultVertex
{
public: 
  GoFlowVertex(int i = -1, float x = 0, float y = 0, float z = 0, float nx = 0, float ny = 0, float nz = 0, 
	       float d = 0, float mx = 0, float my = 0, float mz = 0, float e = 0) 
    : GoDefaultVertex(i,x,y,z,nx,ny,nz)
    {
      value_ = new GoFlowVertexBase<float>(d,mx,my,mz,e);
    }
  GoFlowVertex(GoDefaultVertexBase<float> *t) 
    : GoDefaultVertex(t) 
    {
      value_ = new GoFlowVertexBase<float>;
    }
  
  virtual ~GoFlowVertex() { delete value_; }

  //! ... want some more ?
  virtual GoVertex<float> *clone() const {
    return new GoFlowVertex;
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
  GoFlowVertexBase<float> *value_;
};

#endif // GOFLOWVERTEX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.4  2002/03/18 09:58:55  prkipfer
| refactored element structure
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
