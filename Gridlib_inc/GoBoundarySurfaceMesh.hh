/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOBOUNDARYSURFACEMESH_HH
#define GOBOUNDARYSURFACEMESH_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>
#include <iostream>
#include <typeinfo>
#include <algorithm>

#include "GoSurfaceMesh.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


//! class GoBoundarySurfaceMesh
/*!
  This is a helper class that links vertices instead of copying.
 */
class GoBoundarySurfaceMesh 
  : public GoSurfaceMesh
{
public:
  //! Constructor
  GoBoundarySurfaceMesh();
  virtual ~GoBoundarySurfaceMesh();

  //! This operator pretty prints info about the mesh
//  friend std::ostream& operator<<(std::ostream&, const GoBoundarySurfaceMesh&);
};


//  #ifndef OUTLINE
//  #include "GoBoundarySurfaceMesh.in"
//  #endif

#endif // GOBOUNDARYSURFACEMESH_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2001/02/13 11:17:03  prkipfer
| introduced boundary surface mesh
|
|
+---------------------------------------------------------------------*/
