/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOTRIANGLEMESH_HH
#define GOTRIANGLEMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GoSurfaceMesh.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! class GoTriangleMesh
/*!

 */
class GoTriangleMesh 
  : public GoSurfaceMesh
{
public:
  //! Constructor
  GoTriangleMesh();
  virtual ~GoTriangleMesh();

  //! Query topology
  virtual GoVertex<float> *         nextVertex(GoVertex<float> *v);
  virtual GoGeometryElement<float> *nextElement(GoVertex<float> *v);
//  virtual GoEdge<float> *           nextEdge(GoVertex<float> *v);

  //! Mesh smoothing methods
//  virtual void smoothUmbrella();
//  virtual void smoothCurvatureFlow();

  //! This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoTriangleMesh&);
};


#ifndef OUTLINE
#include "GoTriangleMesh.in"
#endif

#endif // GOTRIANGLEMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.11  2001/09/12 11:53:06  prkipfer
| introduced adaptive tet subdivision
|
| Revision 1.10  2000/12/12 10:03:50  prkipfer
| moved consistency check to GoMesh
|
| Revision 1.9  2000/11/07 17:27:56  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.8  2000/10/31 16:29:07  uflabsik
| new version
|
| Revision 1.7  2000/10/23 09:35:20  prkipfer
| moved bounding box calculation to GoMesh
|
| Revision 1.6  2000/08/22 17:15:24  prkipfer
| added possibility to query mesh contents
|
| Revision 1.5  2000/07/28 15:25:28  uflabsik
| functionality added
|
| Revision 1.4  2000/07/20 13:32:41  prkipfer
| complete rework: now KCC is the standard compiler for Linux
|
| Revision 1.3  2000/06/16 13:32:22  prkipfer
| started mesh class hierarchy
|
| Revision 1.2  2000/06/14 15:39:18  prkipfer
| improved base classes and added funcstruct processing for mesh
|
| Revision 1.1.1.1  2000/06/08 16:24:43  prkipfer
| Imported source tree to start the project
|
|
+---------------------------------------------------------------------*/
