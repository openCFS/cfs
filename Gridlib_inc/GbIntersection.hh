/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBINTERSECTION_HH
#define  GBINTERSECTION_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GbVec3.hh"
#include "GbBox3.hh"
#include "GbSegment.hh"
#include "GbRay.hh"
#include "GbLine.hh"
#include "GbPlane.hh"
#include "GbFrustum.hh"
#include "GbSphere.hh"
#include "GbDistance.hh"

class GoTriangleElement;
class GoTetrahedronElement;

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbIntersection
{
public:
  // return value is 'true' if and only if objects intersect
  static GbBool test (const GbSegment<T>& rkSegment, const GbBox3<T>& rkBox);

  static GbBool test (const GbRay<T>& rkRay, const GbBox3<T>& rkBox);

  static GbBool test (const GbSegment<T>& rkSegment, const GoTriangleElement& rkTriangle);

  static GbBool test (const GbRay<T>& rkRay, const GoTriangleElement& rkTriangle);

  static GbBool test (const GoTriangleElement& rkTriangle, const GbBox3<T>& rkBox);

  static GbBool test (const GbLine<T>& rkLine, const GbBox3<T>& rkBox);

  static GbBool test (const GbLine<T>& rkLine, const GoTriangleElement& rkTriangle);

  static GbBool test (const GbPlane<T>& rkPlane, const GbBox3<T>& rkBox);

  static GbBool test (const GbPlane<T>& rkPlane0, const GbPlane<T>& rkPlane1);

  static GbBool test (const GbSegment<T>& rkSegment, const GbSphere<T>& rkSphere);

  static GbBool test (const GbRay<T>& rkRay, const GbSphere<T>& rkSphere);

  static GbBool test (const GbLine<T>& rkLine, const GbSphere<T>& rkSphere);

  /*! The boolean bUnitNormal is a hint about whether or not the plane normal
      is unit length.  If it is not, the length must be calculated by these
      routines.  For batch calls, the plane normal should be unitized in advance
      to avoid the expensive length calculation. */
  static GbBool test (const GbPlane<T>& rkPlane, const GbSphere<T>& rkSphere,
		      GbBool bUnitNormal);

  static GbBool test (const GbSphere<T>& rkSphere, const GbFrustum<T>& rkFrustum);

  // Test for intersection of static spheres.
  static GbBool test (const GbSphere<T>& rkS0, const GbSphere<T>& rkS1);

  /*! Determine if triangle transversely intersects sphere.  Return value is
    'true' if and only if they intersect.  The Function does not indicate an
    intersection if one or more vertices are the only points of intersection
    or if the triangle is tangent to the sphere. */
  static GbBool test (const GoTriangleElement& rkT, const GbSphere<T>& rkS);

  // This method considers boxes to be stationary
  static GbBool test (const GbBox3<T>& rkBox0, const GbBox3<T>& rkBox1);

  static GbBool test(const GoTriangleElement& triU, const GoTriangleElement& triV);

  static GbBool test(const GbBox3<T>& rkBox, const GbFrustum<T>& rkFrustum);

  /*! Clipping of a linear component 'origin'+t*'direction' against an
      axis-aligned box [-e0,e0]x[-e1,e1]x[-e2,e2] where 'extent'=(e0,e1,e2).
      The values of t0 and t1 must be set by the caller.  If the component is a
      segment, set t0 = 0 and t1 = 1.  If the component is a ray, set t0 = 0 and
      t1 = INFINITY.  If the component is a line, set t0 = -INFINITY and
      t1 = INFINITY.  The values are (possibly) modified by the clipper.
  */
  static GbBool find (const GbVec3<T>& rkOrigin, const GbVec3<T>& rkDirection, 
		      const T afExtent[3], T& rfT0, T& rfT1);


  static GbBool find (const GbSegment<T>& rkSegment, const GbBox3<T>& rkBox,
		      int& riQuantity, GbVec3<T> akPoint[2]);

  static GbBool find (const GbRay<T>& rkRay, const GbBox3<T>& rkBox,
		      int& riQuantity, GbVec3<T> akPoint[2]);

  static GbBool find (const GbSegment<T>& rkSegment, const GoTriangleElement& rkTriangle, 
		      GbVec3<T>& rkPoint);

  static GbBool find (const GbRay<T>& rkRay, const GoTriangleElement& rkTriangle, 
		      GbVec3<T>& rkPoint);
  
  static GbBool find (const GbLine<T>& rkLine, const GbBox3<T>& rkBox,
		      int& riQuantity, GbVec3<T> akPoint[2]);

  static GbBool find (const GbLine<T>& rkLine, const GoTriangleElement& rkTriangle, 
		      GbVec3<T>& rkPoint);

  static GbBool find (const GbPlane<T>& rkPlane0, const GbPlane<T>& rkPlane1,
		      GbLine<T>& rkLine);

  static GbBool find (const GbSegment<T>& rkSegment, const GbSphere<T>& rkSphere, 
		      int& riQuantity, GbVec3<T> akPoint[2]);

  static GbBool find (const GbRay<T>& rkRay, const GbSphere<T>& rkSphere,
		      int& riQuantity, GbVec3<T> akPoint[2]);
  
  static GbBool find (const GbLine<T>& rkLine, const GbSphere<T>& rkSphere,
		      int& riQuantity, GbVec3<T> akPoint[2]);

  /*! Find the intersection of static spheres.  The circle of intersection is
      X(t) = C + R*(cos(t)*U + sin(t)*V) for 0 <= t < 2*pi. */
  static GbBool find (const GbSphere<T>& rkS0, const GbSphere<T>& rkS1,
		      GbVec3<T>& rkU, GbVec3<T>& rkV, GbVec3<T>& rkC, T& rfR);

  /*! If there is a intersection, eventually generated addition tets are returned in rkIntr */
  static GbBool find (const GoTetrahedronElement& rkT0, const GoTetrahedronElement& rkT1,
		      std::vector<GoTetrahedronElement*>& rkIntr);
  
  
  static const T TOLERANCE;
  
private:
  INLINE static GbBool clip(T fDenom, T fNumer, T& rfT0, T& rfT1);

  INLINE static void projection(const GbVec3<T>& rkD, const GbVec3<T> akV[3],
				T& rfMin, T& rfMax);
  INLINE static void projectTriangle (const GbVec3<T>& rkD, const GoTriangleElement& rkTriangle,
				      T& rfMin, T& rfMax);
  INLINE static void projectBox (const GbVec3<T>& rkD, const GbBox3<T>& rkBox, 
				 T& rfMin, T& rfMax);

  INLINE static void splitAndDecompose (GoTetrahedronElement* kTetra, const GbPlane<T>& rkPlane,
					std::vector<GoTetrahedronElement*>& rkInside);

  static GoTetrahedronElement* makeatet(GoTetrahedronElement* kTetra,
					const GbVec3<T>& v0, const GbVec3<T>& v1, const GbVec3<T>& v2, const GbVec3<T>& v3);
};

#ifdef __GNUG__

#include "GbIntersection.in"
#include "GbIntersection.T"

#else

#pragma instantiate GbIntersection<float>
//#pragma instantiate GbIntersection<double>

#ifndef OUTLINE
#include "GbIntersection.in"
#endif

#endif  // g++

#endif  // GBINTERSECTION_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.6  2002/03/18 10:00:16  prkipfer
| added tet-tet intersection
|
| Revision 1.5  2001/12/18 13:26:30  prkipfer
| introduced sphere object
|
| Revision 1.4  2001/09/12 09:16:13  prkipfer
| added tri-box test
|
| Revision 1.3  2001/06/15 08:30:32  prkipfer
| made helper methods private
|
| Revision 1.2  2001/03/20 09:36:10  prkipfer
| added tri-tri test
|
| Revision 1.1  2001/01/03 15:30:08  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
+---------------------------------------------------------------------*/
