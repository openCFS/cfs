/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOBOUNDARYFLOWVERTEX_HH
#define GOBOUNDARYFLOWVERTEX_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoFlowVertex.hh"
#include "GoBoundaryVertexBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoBoundaryFlowVertex

  Vertex class for fluid dynamics boundary meshes
 */
class GoBoundaryFlowVertex
  : public GoFlowVertex
{
public: 
  GoBoundaryFlowVertex(int i = -1, float x = 0, float y = 0, float z = 0, 
		       float d = 0, float mx = 0, float my = 0, float mz = 0, float e = 0)
    : GoFlowVertex(i,x,y,z,d,mx,my,mz,e)
  {
    boundary_ = new GoBoundaryVertexBase<float>();
  }

  virtual ~GoBoundaryFlowVertex() { delete boundary_; }

  //! ... want some more ?
  virtual GoVertex<float> *clone() const {
    return new GoBoundaryFlowVertex;
  }

  //! Link to corresponding data vertex in volume mesh
  virtual void setBoundaryVertex(GoVertex<float> *v) {
    ::setVertex<GoBoundaryVertexBase<float>,float> (boundary_,v);
  }
  virtual GoVertex<float> *getBoundaryVertex() const {
    return ::getVertex<GoBoundaryVertexBase<float>,float> (boundary_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoVertex<float>&);

protected:
  GoBoundaryVertexBase<float> *boundary_;
};

#endif // GOBOUNDARYFLOWVERTEX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/02/13 11:03:17  prkipfer
| introduced boundary vertices
|
|
+---------------------------------------------------------------------*/
