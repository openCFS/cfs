/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBSUBDIVIDEUNIFORM_HH
#define GBSUBDIVIDEUNIFORM_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GoMesh.hh"
#include "GoVertex.hh"
#include "GoGeometryElement.hh"
#include "GoTriangleElement.hh"
#include "GoQuadElement.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class UniformOracle 
  : public GbOracle 
{
public:
//  UniformOracle() { debugmsg("constructor"); }
//  UniformOracle(const UniformOracle& t) { debugmsg("copy constructor"); }
//  virtual ~UniformOracle() { debugmsg("destructor"); }
  virtual GbBool operator()() const {
    // always true -> subdivide everything
    return true;
  }
  virtual GbBool operator()(int l, int id) const {
    // always true -> subdivide everything
    return true;
  }
};

class GbSubdivideUniformFunctor
{
public:
  GbSubdivideUniformFunctor(GoMesh &m, int l=0, 
			    GbSubdivisionStyle st=PRIMAL_STYLE, 
			    GbSubdivisionScheme sch=LINEAR_SCHEME);

  INLINE void operator()(GoGeometryElement<float> *t) {
    if (typeid(*t)==typeid(GoTriangleElement))
#ifdef HAS_DYN_CAST
      subdiv(dynamic_cast<GoTriangleElement*>(t));
#else
      subdiv(static_cast<GoTriangleElement*>(t));
#endif
    else if (typeid(*t)==typeid(GoQuadElement))
#ifdef HAS_DYN_CAST
      subdiv(dynamic_cast<GoQuadElement*>(t));
#else
      subdiv(static_cast<GoQuadElement*>(t));
#endif
    else
      errormsg("cannot subdivide element type "<<typeid(*t).name());
  }

protected:
  void subdiv(GoTriangleElement *t);
  void subdiv(GoQuadElement *t);

  GbVec3<float> getNewPosition(GoTriangleElement *t, int j);
  GbVec3<float> getNewPosition(GoQuadElement *t, int j);

  const UniformOracle oracle_;
  GoMesh &mesh_;
  int currentLevel_;
  GbSubdivisionStyle style_;
  GbSubdivisionScheme scheme_;
};

class GbSubdivideUniform 
{
public:
  GbSubdivideUniform(GoMesh &m, int l=0, 
		     GbSubdivisionStyle st=PRIMAL_STYLE, 
		     GbSubdivisionScheme sch=LINEAR_SCHEME);
  ~GbSubdivideUniform();

protected:
  GoMesh &mesh_;
  int currentLevel_;
  GbSubdivisionStyle style_;
  GbSubdivisionScheme scheme_;

  INLINE void doSqrt3Style(GoTriangleElement* elem);

private:
  // assignment to functor is forbidden
  GbSubdivideUniform& operator=(GbSubdivideUniform& rhs);

  // do not copy
  GbSubdivideUniform(const GbSubdivideUniform& t);
};

#ifndef OUTLINE
#include "GbSubdivideUniform.in"
#endif

#endif // GBSUBDIVIDEUNIFORM_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.5  2001/12/11 12:44:59  prkipfer
| fixes for KCC compiler on new PCs
|
| Revision 1.4  2001/09/12 11:53:02  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.3  2001/01/02 14:27:56  prkipfer
| introduced layered subdivision functors
|
| Revision 1.2  2000/10/17 14:36:22  uflabsik
| new primitives added
|
| Revision 1.1  2000/09/07 16:55:50  prkipfer
| added subdivision
|
|
+---------------------------------------------------------------------*/
