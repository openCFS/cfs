/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GBLATTICE_HH
#define GBLATTICE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include <fstream>

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GbLattice
{
public:
  // abstract base class for GbTImage
  virtual ~GbLattice ();

  // data acesss
  INLINE int getDimensions () const;
  INLINE const int* getBounds () const;
  INLINE int getBound (int uiI) const;
  INLINE int getQuantity () const;
  INLINE const int* getOffsets () const;
  INLINE int getOffset (int uiI) const;

  // assignment
  GbLattice& operator= (const GbLattice& rkLattice);

  // comparisons
  GbBool operator== (const GbLattice& rkLattice) const;
  GbBool operator!= (const GbLattice& rkLattice) const;

  // Conversions between n-dim and 1-dim structures.  The coordinate arrays
  // must have the same number of elements as the dimensions of the lattice.
  int getIndex (const int* auiCoord) const;
  void getCoordinates (int uiIndex, int* auiCoord) const;

  // streaming
  GbBool load (std::ifstream& rkIStr);
  GbBool save (std::ofstream& rkOStr) const;

  static GbBool loadRaw (const char* acFilename, int& ruiDimensions,
			 int*& rauiBound, int& ruiQuantity,
			 int& ruiRTTI, int& ruiSizeOf, char*& racData);
  static GbBool loadVoxels (const char* acFilename, int& ruiDimensions,
			    int*& rauiBound, int& ruiQuantity,
			    int& ruiRTTI, int& ruiSizeOf, char*& racData);

protected:
  // Construction.  GbLattice accepts responsibility for deleting the
  // bound array.
  GbLattice (int uiDimensions, int* auiBound);
  GbLattice (const GbLattice& rkLattice);
  GbLattice ();

  // For deferred creation of bounds.  GbLattice accepts responsibility
  // for deleting the bound array.
  GbLattice (int uiDimensions);
  void setBounds (int* auiBound);
  void computeQuantityAndOffsets ();

  int dimensions_;
  int* iBound_;
  int iSize_;
  int* offset_;

  // streaming
  static const char* magicHeader_;
};

#ifndef OUTLINE
#include "GbLattice.in"
#endif

#endif // GBLATTICE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/06/15 09:26:30  prkipfer
| removed useless signedness
|
| Revision 1.2  2001/04/23 10:02:40  prkipfer
| minor bug fixes
|
| Revision 1.1  2001/03/20 09:50:17  prkipfer
| introduced image handling tool
|
|
+---------------------------------------------------------------------*/
