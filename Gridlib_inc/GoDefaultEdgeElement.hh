/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GODEFAULTEDGEELEMENT_HH
#define GODEFAULTEDGEELEMENT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"
#include "GoVertex.hh"
#include "GoEdgeElement.hh"
#include "GoDefaultEdgeElementBase.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


/*!
  class GoEdgeElement

  Average useful edge class
 */
class GoDefaultEdgeElement 
  : public GoEdgeElement<float>
{
public: 
  GoDefaultEdgeElement(int i = -1) { 
    this_ = new GoDefaultEdgeElementBase<float>(i); 
  }
  GoDefaultEdgeElement(GoDefaultEdgeElementBase<float> *t) : this_(t) {}

  //! ... want some more ?
  virtual GoEdgeElement<float> *clone() const {
    return new GoDefaultEdgeElement;
  }

  virtual ~GoDefaultEdgeElement() { delete this_; }

  // Operations on vertices of the geometry object
  virtual void setVertex(int i, GoVertex<float> *v) {
    ::setVertex<GoDefaultEdgeElementBase<float>,float> (this_,i,v);
  }
  virtual GoVertex<float> *getVertex(int i) const {
    return ::getVertex<GoDefaultEdgeElementBase<float>,float> (this_,i);
  }
  virtual int findVertex(GoVertex<float> *v) const {
    return ::findVertex<GoDefaultEdgeElementBase<float>,float> (this_,v);
  }
  
  // Status flags indicating semantic defined by the
  // object using this class
  virtual void setFlag(GbEdgeStatusFlag f) {
    ::setFlag<GoDefaultEdgeElementBase<float>,float> (this_,f);
  }
  virtual void delFlag(GbEdgeStatusFlag f) {
    ::delFlag<GoDefaultEdgeElementBase<float>,float> (this_,f);
  }
  virtual GbBool testFlag(GbEdgeStatusFlag f) const {
    return ::testFlag<GoDefaultEdgeElementBase<float>,float> (this_,f);
  }

  // Integer to identify the object
  // Has no meaning to this class's implementation
  virtual void setId(int i) {
    ::setId<GoDefaultEdgeElementBase<float>,float> (this_,i);
  }
  virtual int getId() const {
    return ::getId<GoDefaultEdgeElementBase<float>,float> (this_);
  }

  //! The face object, this edge belongs to
  virtual void setElement(GoGeometryElement<T> *f) {
    ::setElement<GoDefaultEdgeElementBase<float>, float> (this_, f);
  }
  virtual GoGeometryElement<T> *getElement() const {
    return ::getElement<GoDefaultEdgeElementBase<float>, float> (this_);
  }

  //! The neighbor valence of the edge
  virtual void setValence(int v) {
    ::setValence<GoDefaultEdgeElementBase<float>, float> (this_, v);
  }
  virtual int getValence() const {
    return ::getValence<GoDefaultEdgeElementBase<float>, float> (this_);
  }
  
  // Get memory statistics and reduce memory space needed
  virtual void statistic() const {
    ::statistic<GoDefaultEdgeElementBase<float>,float> (this_);
  }
  virtual void shrink() {
    ::shrink<GoDefaultEdgeElementBase<float>,float> (this_);
  }

//  friend std::ostream& operator<<(std::ostream&, const GoDefaultEdgeElement&);

private:
  GoDefaultEdgeElementBase<float> *this_;
};

#endif // GOEDGEELEMENT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 15:20:16  prkipfer
| introduced edge elements
|
|
+---------------------------------------------------------------------*/
