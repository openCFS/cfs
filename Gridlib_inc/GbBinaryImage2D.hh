/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBBINARYIMAGE2D_HH
#define  GBBINARYIMAGE2D_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <vector>
#include "GbMath.hh"
#include "GbTList.hh"
#include "GbImages.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/


class GbBinaryImage2D 
  : public GbImageInt2D
{
public:
  // Construction and destruction.  GbBinaryImage2D accepts responsibility for
  // deleting the input data array.
  GbBinaryImage2D (int uiXBound, int uiYBound, GbEint* atData = 0);
  GbBinaryImage2D (const GbBinaryImage2D& rkImage);
  GbBinaryImage2D (const char* acFilename);

  // Extract boundaries from blobs in a binary image.  The input image is
  // treated as binary in the sense that pixels are either zero or not zero.
  // The blobs correspond to the non-zero values.  The boundaries are
  // treated as 8-connected (neighbors can occur in any of 8 directions from
  // the given pixel).  The blob interiors are treated as 4-connected.
  //
  // For sake of naming, let the returned array be called Boundary.  If the
  // image is identically zero, there are no blobs, so Boundary = null.  If
  // there are blobs, the boundaries are packed into Boundary as follows.
  //
  // Boundary[0] is the number of blobs (connected components).  Boundary[1]
  // is the number of pixels in the first boundary.  Boundary[2] through
  // Boundary[Boundary[1]+1] are the pixel locations stored as unsigned
  // integers (see below).  The next boundary (if any) is stored in the
  // same format:  number of pixel locations followed by the array of pixel
  // locations.
  //
  // The pixel locations are stored as unsigned integers.  If I is the index
  // and (x,y) is the pixel location, then I = x + XDIM*y where the image is
  // XDIM-by-YDIM.  You can use rkImage.GetCoordinates(I,x,y) to extract the
  // (x,y) values for a given I.
  //
  // The application is responsible for deleting the returned array.

  int* getBoundaries () const;


  // Compute the connected components of a binary image.  The components in
  // the returned image are labeled with positive integer values.  If the
  // image is identically zero, then the components image is identically
  // zero and the returned quantity is zero.

  void getComponents (int& riQuantity, GbImageInt2D& rkComponents) const;


  // Compute the L1 distance transform.  Given a pixel (x,y), the neighbors
  // (x+1,y), (x-1,y), (x,y+1), and (x,y-1) are 1 unit of distance from
  // (x,y).  The neighbors (x+1,y+1), (x+1,y-1), (x-1,y+1), and (x-1,y-1)
  // are 2 units of distance from (x,y).

  void getL1Distance (int& riMaxDistance, GbImageInt2D& rkTransform) const;

  // Compute the L2 distance transform (Euclidean distance transform).  The
  // distances are exact as long as they are smaller than 100 (see the
  // comments in the source code).
  void getL2Distance (float& rfMaxDistance, GbImageDouble2D& rkTransform) const;

  // Compute a skeleton of the image.  Pixels are trimmed from outside to
  // inside using L1 distance.  The maximum distance from the L1 transform
  // is returned in case the application needs this.  Connectivity and
  // topology of the original blobs are preserved.
  void getSkeleton (int& riMaxDistance, GbImageInt2D& rkSkeleton) const;

protected:
  // helper for boundary extraction
  typedef std::vector<int> BoundaryList;
  BoundaryList* extractBoundary (int uiX0, int uiY0,
				 GbImageInt2D& rkTemp) const;

  // helper for component labeling
  void addToAssociative (int iI0, int iI1, int* aiAssoc) const;

  // helper for L1 distance (no embedding in larger image, done in-place)
  void getL1DistanceZeroBoundary (int& riMaxDistance, GbImageInt2D& rkTemp) const;

  // helpers for L2 distance
  void L2Initialize (GbImageInt2D& rkXNear, GbImageInt2D& rkYNear,
		     GbImageInt2D& rkDist) const;

  void L2Check (int iX, int iY, int iDx, int iDy, GbImageInt2D& rkXNear,
		GbImageInt2D& rkYNear, GbImageInt2D& rkDist) const;

  void L2XpYp (GbImageInt2D& rkXNear, GbImageInt2D& rkYNear,
	       GbImageInt2D& rkDist) const;

  void L2XpYm (GbImageInt2D& rkXNear, GbImageInt2D& rkYNear,
	       GbImageInt2D& rkDist) const;

  void L2XmYp (GbImageInt2D& rkXNear, GbImageInt2D& rkYNear,
	       GbImageInt2D& rkDist) const;

  void L2XmYm (GbImageInt2D& rkXNear, GbImageInt2D& rkYNear,
	       GbImageInt2D& rkDist) const;

  void L2Finalize (const GbImageInt2D& rkDist, float& rfMaxDistance,
		   GbImageDouble2D& rkTransform) const;

  // helpers for thinning
  void trimContour (int iMaxDistance, const GbImageInt2D& rkTemp) const;
  void trimSkeleton (int iMaxDistance, const GbImageInt2D& rkTemp) const;
  int convertToByte (int uiX, int uiY, const GbImageInt2D& rkTemp) const;
  int pixelType4 (int uiX, int uiY, const GbImageInt2D& rkTemp) const;
  GbBool isBranch (int uiX, int uiY, const GbImageInt2D& rkTemp) const;
  GbBool notNeeded (int uiX, int uiY, const GbImageInt2D& rkTemp) const;
};

// #ifndef OUTLINE
// #include "GbBinaryImage2D.in"
// #endif

#endif  // GBBINARYIMAGE2D_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/06/19 16:27:16  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.2  2001/06/15 09:24:16  prkipfer
| introduced L2 distance, switched to STL containers, added compiler hints and removed useless signedness
|
| Revision 1.1  2001/04/23 10:02:22  prkipfer
| introduced morphological image operations and TList class
|
|
+---------------------------------------------------------------------*/
