/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBTIMAGE3D_HH
#define  GBTIMAGE3D_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>

#include "GbTypes.hh"
#include "GbTImage.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


template <class T>
class GbTImage3D 
  : public GbTImage<T>
{
public:
  // Construction and destruction.  GbTImage3D accepts responsibility for
  // deleting the input data array.
  GbTImage3D (int uiXBound, int uiYBound, int uiZBound, T* atData = NULL);
  GbTImage3D (const GbTImage3D& rkImage);
  GbTImage3D (const char* acFilename);

  // data access
  INLINE T& operator() (int uiX, int uiY, int uiZ) const;

  // conversion between 3D coordinates and 1D indexing
  INLINE int getIndex (int uiX, int uiY, int uiZ) const;
  INLINE void getCoordinates (int uiIndex, int& ruiX, int& ruiY, int& ruiZ) const;

  // assignment
  INLINE GbTImage3D& operator= (const GbTImage3D& rkImage);
  INLINE GbTImage3D& operator= (T tValue);
};

#ifdef __GNUG__

#include "GbTImage3D.in"
#include "GbTImage3D.T"

#else

#ifndef OUTLINE
#include "GbTImage3D.in"
#endif

#endif  // g++

#endif  // GBTIMAGE3D_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:27:31  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:21  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
