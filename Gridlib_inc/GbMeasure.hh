/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBMEASURE_HH
#define  GBMEASURE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbMath.hh"
#include "GbTypes.hh"
#include "GbVec3.hh"
#include "GoGeometryElement.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbMeasure
{
public:
  typedef struct {
    T x, y, z;         // pixel coordinates
  } Point3;

  typedef struct {
    int x, y, z;       // pixel coordinates 
    int nx, ny, nz;    // normal vector
  } Point;

  static T surface(const GoGeometryElement<T>* el);
  static T volume(const GoGeometryElement<T>* el);
  static void volume(int n, Point* p, T dx, T dy, T dz, T& surfaceArea, T& volume);

  static T polyArea(int n, const Point3* p);
  static T polyVolume(int n, const Point3* p, int polys, const int* connectivity);

protected:
  static int tetraconn_[];
  static int hexaconn_[];
  static int pyraconn_[];
  static int prismconn_[];
  static int octaconn_[];
};

#ifdef __GNUG__

//#include "GbMeasure.in"
#include "GbMeasure.T"

#else

#pragma instantiate GbMeasure<float>
#pragma instantiate GbMeasure<double>

//#ifndef OUTLINE
//#include "GbMeasure.in"
//#endif

#endif  // g++

#endif  // GBMEASURE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/05/03 11:06:03  prkipfer
| made connectivity data static, added missing type cases
|
| Revision 1.1  2001/04/23 10:07:47  prkipfer
| added measurements class
|
|
+---------------------------------------------------------------------*/
