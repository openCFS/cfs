/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBTIMAGE_HH
#define  GBTIMAGE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <fstream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbLattice.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


/*! The class T is intended to be a wrapper of native data (int, float, char*,
    etc.).  The code uses memcpy, memcmp, and memset on the array of T values.
    Class T must have the following member functions:
    T::T ()
    T& T::operator= (T)
    static const char* getRTTI ()
    The static member function returns a string that is used for streaming.
*/

template <class T>
class GbTImage 
  : public GbLattice
{
public:
  GbTImage (int uiDimensions, int* auiBound, T* atData = NULL);
  GbTImage (const GbTImage& rkImage);
  GbTImage (const char* acFilename);
  virtual ~GbTImage ();

  // data access
  INLINE T* getData () const;
  INLINE T& operator[] (int uiIndex) const;

  // assignment
  INLINE GbTImage& operator= (const GbTImage& rkImage);
  INLINE GbTImage& operator= (T tValue);

  // comparison
  INLINE GbBool operator== (const GbTImage& rkImage) const;
  INLINE GbBool operator!= (const GbTImage& rkImage) const;

  // streaming
  GbBool load (const char* acFilename);
  GbBool save (const char* acFilename) const;

protected:
  // for deferred creation of bounds
  GbTImage (int uiDimensions);
  void setData (T* atData);

  T* data_;
};

#ifdef __GNUG__

#include "GbTImage.in"
#include "GbTImage.T"

#else

#ifndef OUTLINE
#include "GbTImage.in"
#endif

#endif  // g++

#endif  // GBTIMAGE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:27:31  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:18  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
