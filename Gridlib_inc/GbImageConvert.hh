/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBIMAGECONVERT_HH
#define GBIMAGECONVERT_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

void GbImageConvert (int uiQuantity, int uiSrcRTTI, void* pvSrcData, int uiTrgRTTI, void* pvTrgData);


//#ifndef OUTLINE
//#include "GbImageConvert.in"
//#endif

#endif // GBIMAGECONVERT_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:25:00  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:14  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
