/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GODEFAULTVERTEX_HH
#define GODEFAULTVERTEX_HH

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

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

/*!
  class GoDefaultVertex

  Average useful vertex class
 */
class GoDefaultVertex
  : public GoVertex<float>
{
public: 
  GoDefaultVertex(int i = -1, float x = 0, float y = 0, float z = 0, float nx = 0, float ny = 0, float nz = 0) {
    this_ = new GoDefaultVertexBase<float>(i, x, y, z, nx, ny ,nz); 
  }
  GoDefaultVertex(GoDefaultVertexBase<float> *t) : this_(t) {}

  //! ... want some more ?
  virtual GoVertex<float> *clone() const {
    return new GoDefaultVertex;
  }

  virtual ~GoDefaultVertex() { delete this_; }

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
  // compute normal at this vertex by averaging the normals
  // of neighboring face objects
  virtual void computeNormal(std::vector<GoGeometryElement<float> *>& S) {
    GbVec3<float> n(0.0f, 0.0f, 0.0f);

    if (S.empty()) {
      if (getElement()) {
	setNormal(getElement()->getNormal());
	warningmsg("no neighbors ? possible mesh consistency problem or mesh not closed ?");
      }
      else {
	warningmsg("no face associated ? possible mesh setup problem ?");
      }
      return;
    }
  
    for (std::vector<GoGeometryElement<float> *>::iterator sIter = S.begin(); 
	 sIter != S.end(); 
	 ++sIter)
      
      n += (*sIter)->getNormal();

    n.normalize();
    setNormal(n);
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

//  friend std::ostream& operator<<(std::ostream&, const GoVertex<float>&);

protected:
  GoDefaultVertexBase<float> *this_;
};

#endif // GODEFAULTVERTEX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.4  2002/03/18 09:58:55  prkipfer
| refactored element structure
|
| Revision 1.3  2001/02/13 11:03:16  prkipfer
| introduced boundary vertices
|
| Revision 1.2  2001/01/02 15:21:34  prkipfer
| introduced cloning and support for new base classes
|
| Revision 1.1  2000/11/07 17:27:25  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
