/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBGAUSSPOINTSFIT_HH
#define  GBGAUSSPOINTSFIT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GbVec3.hh"
#include "GbBox3.hh"
#include "GoMesh.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

// Original Source: Magic Software, Inc.
// FREE SOURCE CODE

/*! Fit points with a Gaussian distribution.  The center is the mean of the
    points, the axes are the eigenvectors of the covariance matrix, and the
    extents are the eigenvalues of the covariance matrix and are returned in
    increasing order.
*/
GbBox3<float> 
GbGaussPointsFit (const GoMesh* mesh, int level=0);

/*! This function allows selection of valid
    vertices from a pool.  The return value is 'true' if and only if at least
    one vertex was valid.
*/
GbBool 
GbGaussPointsFit (const GoMesh* mesh, GbVertexFlag flag, GbBox3<float>& box, int level=0);


//#ifndef OUTLINE
//#include "GbGaussPointsFit.in"
//#endif

#endif  // GBGAUSSPOINTSFIT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/09/12 09:15:35  prkipfer
| made GaussPointsFit work on a specific mesh level
|
| Revision 1.1  2001/01/02 15:04:01  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
