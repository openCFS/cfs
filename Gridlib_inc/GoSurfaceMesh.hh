/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOSURFACEMESH_HH
#define GOSURFACEMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <typeinfo>
#include <algorithm>

#include "GoVertex.hh"
#include "GoGeometryElement.hh"
#include "GoMesh.hh"
#include "GbMeshFunctions.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


//! class GoSurfaceMesh
/*!

 */
class GoSurfaceMesh 
  : public GoMesh
{
public:
  //! Constructor
  GoSurfaceMesh();
  virtual ~GoSurfaceMesh();

  //! Query topology
  virtual GoVertex<float> *nextVertex(GoVertex<float> *v);
  virtual GoGeometryElement<float>   *nextElement(GoVertex<float> *v);

  //! Actions on the mesh
  virtual void setupNeighbours();
  virtual void switchNormals();

  //! Mesh smoothing methods
  virtual void smoothUmbrella(int level=0);
  virtual void smoothCurvatureFlow(int level=0);

  virtual std::vector<GoGeometryElement<float> *> * getAllNeighbourElements(GoVertex<float>* vertex);


  //! This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoSurfaceMesh&);

private:
  struct TableEntry {
    int i;
    GoGeometryElement<float> *f;
    TableEntry *next;
  };
  
  void searchTable(int, int, GoGeometryElement<float> *, TableEntry **);
};


#ifndef OUTLINE
#include "GoSurfaceMesh.in"
#endif

#endif // GOSURFACEMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.7  2001/05/03 08:05:51  uflabsik
| functionality added
|
| Revision 1.6  2001/02/13 11:12:55  prkipfer
| fixed wrong RCS id's
|
| Revision 1.5  2000/12/12 10:03:49  prkipfer
| moved consistency check to GoMesh
|
| Revision 1.4  2000/11/07 17:27:55  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.3  2000/10/31 16:29:07  uflabsik
| new version
|
| Revision 1.2  2000/10/23 09:35:20  prkipfer
| moved bounding box calculation to GoMesh
|
| Revision 1.1  2000/09/26 09:26:39  uflabsik
| initital version
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
