/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBCONTAINERBOX_HH
#define  GBCONTAINERBOX_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>
//#include <cstdio>

#include "GbMath.hh"
#include "GbBox3.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

//! Test for containment.  Let X = C + y0*U0 + y1*U1 + y2*U2 where C is the
//! box center and U0, U1, U2 are the orthonormal axes of the box.  X is in
//! the box if |y_i| <= E_i for all i where E_i are the extents of the box.
//! If an epsilon e is supplied, the tests are |y_i| <= E_i + e.
GbBool 
GbInBox (const GbVec3<float>& rkPoint, const GbBox3<float>& rkBox,
	 float fEpsilon = 0.0);

//! Compute union of two boxes
GbBox3<float> 
GbMergedBox (const GbBox3<float>& rkBox0, const GbBox3<float>& rkBox1);

//! Compute minimum volume oriented box containing the specified points
GbBox3<float>
GbMinVolumeBox (int iQuantity, const GbVec3<float>* akPoint);

//#ifndef OUTLINE
//#include "GbContainerBox.in"
//#endif

#endif  // GBCONTAINERBOX_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.3  2001/03/20 09:36:45  prkipfer
| added in-box-test
|
| Revision 1.2  2001/01/09 11:12:53  prkipfer
| new feature: compute minimum volume box
|
| Revision 1.1  2001/01/02 14:59:55  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
