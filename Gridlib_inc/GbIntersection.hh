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
#include "GoTriangleElement.hh"
#include "GbDistance.hh"

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

  static GbBool test (const GbLine<T>& rkLine, const GbBox3<T>& rkBox);

  static GbBool test (const GbLine<T>& rkLine, const GoTriangleElement& rkTriangle);

  static GbBool test (const GbPlane<T>& rkPlane, const GbBox3<T>& rkBox);

  static GbBool test (const GbPlane<T>& rkPlane0, const GbPlane<T>& rkPlane1);

  // This method considers boxes to be stationary
  static GbBool test (const GbBox3<T>& rkBox0, const GbBox3<T>& rkBox1);

  static GbBool test(const GoTriangleElement& triU, const GoTriangleElement& triV);

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

  static const T TOLERANCE;

private:
  INLINE static GbBool clip(T fDenom, T fNumer, T& rfT0, T& rfT1);

  INLINE static void projection(const GbVec3<T>& rkD, const GbVec3<T> akV[3],
				T& rfMin, T& rfMax);
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
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
