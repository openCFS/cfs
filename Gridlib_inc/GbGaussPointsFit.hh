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
GbGaussPointsFit (const GoMesh* mesh);

/*! This function allows selection of valid
    vertices from a pool.  The return value is 'true' if and only if at least
    one vertex was valid.
*/
GbBool 
GbGaussPointsFit (const GoMesh* mesh, GbVertexFlag flag, GbBox3<float>& box);


//#ifndef OUTLINE
//#include "GbGaussPointsFit.in"
//#endif

#endif  // GBGAUSSPOINTSFIT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 15:04:01  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
