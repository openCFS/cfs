/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOSCALARVERTEX_HH
#define GOSCALARVERTEX_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoDefaultVertex.hh"
#include "GoScalarVertexBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoScalarVertex

  Vertex class for meshes with scalar values at the vertices
 */
class GoScalarVertex
  : public GoDefaultVertex
{
public: 
  GoScalarVertex(int i = -1, float x = 0, float y = 0, float z = 0, float nx = 0, float ny = 0, float nz = 0, float v = 0) 
    : GoDefaultVertex(i,x,y,z,nx,ny,nz) 
    {
      value_ = new GoScalarVertexBase<float>(v);
    }
  GoScalarVertex(GoDefaultVertexBase<float> *t) 
    : GoDefaultVertex(t) 
    {
      value_ = new GoScalarVertexBase<float>;
    }

  virtual ~GoScalarVertex() { delete value_; }

  //! ... want some more ?
  virtual GoVertex<float> *clone() const {
    return new GoScalarVertex;
  }

  //! Handle the scalar value
  virtual void setValue(float v) {
    ::setValue<GoScalarVertexBase<float>,float> (value_,v);
  }
  virtual float getValue() const {
    return ::getValue<GoScalarVertexBase<float>,float> (value_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoVertex<float>&);

protected:
  GoScalarVertexBase<float> *value_;
};

#endif // GOSCALARVERTEX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.4  2002/03/18 09:58:56  prkipfer
| refactored element structure
|
| Revision 1.3  2001/02/13 11:03:16  prkipfer
| introduced boundary vertices
|
| Revision 1.2  2001/01/02 15:21:35  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.1  2000/11/07 17:27:27  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
