/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBTIMAGE2D_HH
#define  GBTIMAGE2D_HH

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
class GbTImage2D 
  : public GbTImage<T>
{
public:
  // Construction and destruction.  GbTImage2D accepts responsibility for
  // deleting the input data array.
  GbTImage2D (int uiXBound, int uiYBound, T* atData = NULL);
  GbTImage2D (const GbTImage2D& rkImage);
  GbTImage2D (const char* acFilename);

  // data access
  INLINE T& operator() (int uiX, int uiY) const;

  // conversion between 2D coordinates and 1D indexing
  INLINE int getIndex (int uiX, int uiY) const;
  INLINE void getCoordinates (int uiIndex, int& ruiX, int& ruiY) const;

  // assignment
  INLINE GbTImage2D& operator= (const GbTImage2D& rkImage);
  INLINE GbTImage2D& operator= (T tValue);
};

#ifdef __GNUG__

#include "GbTImage2D.in"
#include "GbTImage2D.T"

#else

#ifndef OUTLINE
#include "GbTImage2D.in"
#endif

#endif  // g++

#endif  // GBTIMAGE2D_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:27:30  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:19  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
