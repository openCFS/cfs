/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBDISTANCE_HH
#define  GBDISTANCE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbMath.hh"
#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GbBox3.hh"
#include "GbFrustum.hh"
#include "GbRay.hh"
#include "GbSegment.hh"
#include "GbLine.hh"
#include "GoTriangleElement.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbDistance
{
public:
  // squared distance measurements
  static T squareDistance (const GbVec3<T>& rkPoint, const GbBox3<T>& rkBox,
			   T* pfBParam0=NULL, T* pfBParam1=NULL, T* pfBParam2=NULL);

  static T squareDistance (const GbVec3<T>& rkPoint, const GoTriangleElement& rkTri,
			   T* pfSParam=NULL, T* pfTParam=NULL);

  static T squareDistance (const GbVec3<T>& rkPoint, const GbRay<T>& rkRay,
			   T* pfParam = NULL);

  static T squareDistance (const GbVec3<T>& rkPoint, const GbSegment<T>& rkSegment, 
			   T* pfParam = NULL);

  static T squareDistance (const GbVec3<T>& rkPoint, const GbLine<T>& rkLine,
			   T* pfParam = NULL);

  static T squareDistance (const GbLine<T>& rkLine0, const GbLine<T>& rkLine1,
			   T* pfLinP0 = NULL, T* pfLinP1 = NULL);

  static T squareDistance (const GbLine<T>& rkLine, const GbRay<T>& rkRay,
			   T* pfLinP = NULL, T* pfRayP = NULL);

  static T squareDistance (const GbLine<T>& rkLine, const GbSegment<T>& rkSeg,
			   T* pfLinP = NULL, T* pfSegP = NULL);

  static T squareDistance (const GbLine<T>& rkLine, const GbBox3<T>& rkBox,
			   T* pfLParam = NULL, T* pfBParam0 = NULL, T* pfBParam1 = NULL, T* pfBParam2 = NULL);

  static T squareDistance (const GbLine<T>& rkLine, const GoTriangleElement& rkTri,
			   T* pfLinP = NULL, T* pfTriP0 = NULL, T* pfTriP1 = NULL);

  static T squareDistance (const GoTriangleElement& rkTri0, const GoTriangleElement& rkTri1, 
			   T* pfTri0P0=NULL, T* pfTri0P1=NULL, T* pfTri1P0=NULL, T* pfTri1P1=NULL);

  static T squareDistance (const GbSegment<T>& rkSeg, const GoTriangleElement& rkTri,
			   T* pfSegP=NULL, T* pfTriP0=NULL, T* pfTriP1=NULL);

  static T squareDistance (const GbSegment<T>& rkSeg0, const GbSegment<T>& rkSeg1,
			   T* pfSegP0=NULL, T* pfSegP1=NULL);

  static T squareDistance (const GbSegment<T>& rkSeg, const GbBox3<T>& rkBox,
 			   T* pfLParam = NULL, T* pfBParam0 = NULL, T* pfBParam1 = NULL, T* pfBParam2 = NULL);

  static T squareDistance (const GbRay<T>& rkRay, const GbBox3<T>& rkBox,
 			   T* pfLParam = NULL, T* pfBParam0 = NULL, T* pfBParam1 = NULL, T* pfBParam2 = NULL);

  static T squareDistance (const GbRay<T>& rkRay0, const GbRay<T>& rkRay1,
			   T* pfRayP0=NULL, T* pfRayP1=NULL);

  static T squareDistance (const GbRay<T>& rkRay, const GbSegment<T>& rkSeg,
			   T* pfRayP=NULL, T* pfSegP=NULL);

  static T squareDistance (const GbRay<T>& rkRay, const GoTriangleElement& rkTri,
			   T* pfRayP=NULL, T* pfTriP0=NULL, T* pfTriP1=NULL);

  static T squareDistance (const GbVec3<T>& rkPoint, const GbFrustum<T>& rkFrustum,
			   GbVec3<T>* pkClosest = NULL);


  // distance measurements
  static T distance (const GbVec3<T>& rkPoint, const GbBox3<T>& rkBox,
		     T* pfBParam0=NULL, T* pfBParam1=NULL, T* pfBParam2=NULL);

  static T distance (const GbVec3<T>& rkPoint, const GoTriangleElement& rkTri,
		     T* pfSParam=NULL, T* pfTParam=NULL);

  static T distance (const GbVec3<T>& rkPoint, const GbRay<T>& rkRay,
		     T* pfParam = NULL);

  static T distance (const GbVec3<T>& rkPoint, const GbSegment<T>& rkSegment,
		     T* pfParam = NULL);

  static T distance (const GbVec3<T>& rkPoint, const GbLine<T>& rkLine,
		     T* pfParam = NULL);

  static T distance (const GbLine<T>& rkLine0, const GbLine<T>& rkLine1,
		     T* pfLinP0 = NULL, T* pfLinP1 = NULL);

  static T distance (const GbLine<T>& rkLine, const GbRay<T>& rkRay,
		     T* pfLinP = NULL, T* pfRayP = NULL);

  static T distance (const GbLine<T>& rkLine, const GbSegment<T>& rkSeg,
		     T* pfLinP = NULL, T* pfSegP = NULL);

  static T distance (const GbLine<T>& rkLine, const GbBox3<T>& rkBox,
		     T* pfLParam = NULL, T* pfBParam0 = NULL, T* pfBParam1 = NULL, T* pfBParam2 = NULL);

  static T distance (const GbLine<T>& rkLine, const GoTriangleElement& rkTri,
		     T* pfLinP = NULL, T* pfTriP0 = NULL, T* pfTriP1 = NULL);

  static T distance (const GoTriangleElement& rkTri0, const GoTriangleElement& rkTri1,
		     T* pfTri0P0=NULL, T* pfTri0P1=NULL, T* pfTri1P0=NULL, T* pfTri1P1=NULL);

  static T distance (const GbSegment<T>& rkSeg, const GbBox3<T>& rkBox,
 		     T* pfLParam = NULL, T* pfBParam0 = NULL, T* pfBParam1 = NULL, T* pfBParam2 = NULL);

  static T distance (const GbSegment<T>& rkSeg, const GoTriangleElement& rkTri,
		     T* pfSegP=NULL, T* pfTriP0=NULL, T* pfTriP1=NULL);

  static T distance (const GbSegment<T>& rkSeg0, const GbSegment<T>& rkSeg1,
		     T* pfSegP0=NULL, T* pfSegP1=NULL);

  static T distance (const GbRay<T>& rkRay, const GbBox3<T>& rkBox,
 		     T* pfLParam = NULL, T* pfBParam0 = NULL, T* pfBParam1 = NULL, T* pfBParam2 = NULL);

  static T distance (const GbRay<T>& rkRay0, const GbRay<T>& rkRay1,
		     T* pfRayP0=NULL, T* pfRayP1=NULL);

  static T distance (const GbRay<T>& rkRay, const GbSegment<T>& rkSeg,
		     T* pfRayP=NULL, T* pfSegP=NULL);

  static T distance (const GbRay<T>& rkRay, const GoTriangleElement& rkTri,
		     T* pfRayP=NULL, T* pfTriP0=NULL, T* pfTriP1=NULL);

  //! Distance from point to plane: Dot(N,X-A) = 0. The plane normal N must be unit length.
  static T distance (const GbVec3<T>& rkPoint, const GbVec3<T>& rkNormal,
		     const GbVec3<T>& rkOrigin, GbVec3<T>* pkClosest = NULL);

  //! Distance from point to plane: Dot(N,X) = c. The plane normal N must be unit length.
  static T distance (const GbVec3<T>& rkPoint, const GbVec3<T>& rkNormal,
		     T fConstant, GbVec3<T>* pkClosest = NULL);

  static T distance (const GbVec3<T>& rkPoint, const GbFrustum<T>& rkFrustum,
		     GbVec3<T>* pkClosest = NULL);

  static const T TOLERANCE;

private:
  static void case000 (GbVec3<T>& rkPnt, const GbBox3<T>& rkBox, T& rfSqrDistance);
  static void case00 (int i0, int i1, int i2, 
		      GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, const GbBox3<T>& rkBox, 
		      T* pfLParam, T& rfSqrDistance);
  static void case0 (int i0, int i1, int i2, 
		     GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, const GbBox3<T>& rkBox, 
		     T* pfLParam, T& rfSqrDistance);
  static void caseNoZeros (GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, const GbBox3<T>& rkBox, 
			   T* pfLParam, T& rfSqrDistance);
  static void face (int i0, int i1, int i2, 
		    GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, 
		    const GbBox3<T>& rkBox, const GbVec3<T>& rkPmE,
		    T* pfLParam, T& rfSqrDistance);
};

#ifdef __GNUG__

//#include "GbDistance.in"
#include "GbDistance.T"

#else

#pragma instantiate GbDistance<float>
//#pragma instantiate GbDistance<double>

//#ifndef OUTLINE
//#include "GbDistance.in"
//#endif

#endif  // g++

#endif  // GBDISTANCE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.3  2001/12/18 13:26:30  prkipfer
| introduced sphere object
|
| Revision 1.2  2001/06/18 10:43:03  prkipfer
| introduced Frustum object
|
| Revision 1.1  2001/01/02 15:07:14  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
