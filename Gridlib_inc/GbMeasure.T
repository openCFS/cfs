// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/08/16 16:53:19  prkipfer
| improved type safety for template parameter
|
| Revision 1.2  2001/05/03 11:06:02  prkipfer
| made connectivity data static, added missing type cases
|
| Revision 1.1  2001/04/23 10:07:47  prkipfer
| added measurements class
|
|
|
+---------------------------------------------------------------------*/

#include <typeinfo>

#include "GoTriangleElement.hh"
#include "GoQuadElement.hh"
#include "GoTetrahedronElement.hh"
#include "GoHexahedronElement.hh"
#include "GoPrismElement.hh"
#include "GoPyramidElement.hh"
#include "GoOctahedronElement.hh"

template <class T>
int
GbMeasure<T>::tetraconn_[16] = { 3, 0,1,3,
				 3, 1,2,3,
				 3, 2,0,3,
				 3, 0,2,1 };
template <class T>
int
GbMeasure<T>::hexaconn_[30] = { 4,  3,0,4,7,
				4,  2,3,7,6,
			        4,  1,2,6,5,
			        4,  0,1,5,4,
			        4,  4,5,6,7,
			        4,  3,2,1,0 };
template <class T>
int
GbMeasure<T>::pyraconn_[21] = { 3, 0,1,4,
				3, 1,2,4,
				3, 2,3,4,
			        3, 3,0,4,
			        4, 3,2,1,0 };
template <class T>
int
GbMeasure<T>::prismconn_[23] = { 4, 1,2,5,4, 
				 4, 2,0,3,5, 
				 4, 0,1,4,3, 
				 3, 2,1,0,
				 3, 3,4,5 };
template <class T>
int
GbMeasure<T>::octaconn_[32] = {3, 0,1,4,
			       3, 1,2,4,
			       3, 2,3,4,
			       3, 3,0,4,
			       3, 1,0,5,
			       3, 2,1,5,
			       3, 3,2,5,
			       3, 0,3,5 };

template <class T>
T
GbMeasure<T>::surface(const GoGeometryElement<T>* el)
{
  T retval=T(0.0);
  Point3 *p=NULL;

  if ( (typeid(*el) == typeid(GoTriangleElement)) ||
       (typeid(*el) == typeid(GoQuadElement)) ) {
    int n=el->getNumVertices();
    p = new Point3[n];
    for (int i=0; i<n; ++i)
      el->getVertex(i)->getPosition().getVec3(p[i].x,p[i].y,p[i].z);
    retval = GbMeasure<T>::polyArea(n,p);
  }
  else if (typeid(*el) == typeid(GoTetrahedronElement)) {
    int count=0;
    for (int t=0; t<4; ++t) {
      int n=tetraconn_[count++];
      p = new Point3[n];
      for (int tt=0; tt<n; ++tt)
	el->getVertex(tetraconn_[count++])->getPosition().getVec3(p[tt].x,p[tt].y,p[tt].z);
      retval += GbMeasure<T>::polyArea(n,p);
    }
  }
  else if (typeid(*el) == typeid(GoHexahedronElement)) {
    int count=0;
    for (int t=0; t<6; ++t) {
      int n=hexaconn_[count++];
      p = new Point3[n];
      for (int tt=0; tt<n; ++tt)
	el->getVertex(hexaconn_[count++])->getPosition().getVec3(p[tt].x,p[tt].y,p[tt].z);
      retval += GbMeasure<T>::polyArea(n,p);
    }
  }
  else if (typeid(*el) == typeid(GoPrismElement)) {
    int count=0;
    for (int t=0; t<5; ++t) {
      int n=prismconn_[count++];
      p = new Point3[n];
      for (int tt=0; tt<n; ++tt)
	el->getVertex(prismconn_[count++])->getPosition().getVec3(p[tt].x,p[tt].y,p[tt].z);
      retval += GbMeasure<T>::polyArea(n,p);
    }
  }
  else if (typeid(*el) == typeid(GoPyramidElement)) {
    int count=0;
    for (int t=0; t<5; ++t) {
      int n=pyraconn_[count++];
      p = new Point3[n];
      for (int tt=0; tt<n; ++tt)
	el->getVertex(pyraconn_[count++])->getPosition().getVec3(p[tt].x,p[tt].y,p[tt].z);
      retval += GbMeasure<T>::polyArea(n,p);
    }
  }
  else if (typeid(*el) == typeid(GoOctahedronElement)) {
    int count=0;
    for (int t=0; t<8; ++t) {
      int n=octaconn_[count++];
      p = new Point3[n];
      for (int tt=0; tt<n; ++tt)
	el->getVertex(octaconn_[count++])->getPosition().getVec3(p[tt].x,p[tt].y,p[tt].z);
      retval += GbMeasure<T>::polyArea(n,p);
    }
  }
  else {
    std::cerr<<"cannot compute surface for element type "<<typeid(*el).name()<<std::endl;
  }

  delete[] p;
  return retval;
}

template <class T>
T
GbMeasure<T>::volume(const GoGeometryElement<T>* el)
{
  T retval=T(0.0);
  if ( (typeid(*el) != typeid(GoTriangleElement)) &&
       (typeid(*el) != typeid(GoQuadElement)) ) {
    // its a 3D object
    Point3 *p=NULL;
    int n=el->getNumVertices();
    p = new Point3[n];
    for (int i=0; i<n; ++i)
      el->getVertex(i)->getPosition().getVec3(p[i].x,p[i].y,p[i].z);

    int poly=0;
    int *conn=NULL;

    if (typeid(*el) == typeid(GoTetrahedronElement)) {
      poly=4;
      conn=tetraconn_;
    }
    else if (typeid(*el) == typeid(GoHexahedronElement)) {
      poly=6;
      conn=hexaconn_;
    }
    else if (typeid(*el) == typeid(GoPrismElement)) {
      poly=5;
      conn=prismconn_;
    }
    else if (typeid(*el) == typeid(GoPyramidElement)) {
      poly=5;
      conn=pyraconn_;
    }
    else if (typeid(*el) == typeid(GoOctahedronElement)) {
      poly=8;
      conn=octaconn_;
    }
    else {
      std::cerr<<"cannot compute volume for element type "<<typeid(*el).name()<<std::endl;
    }

    if(conn) 
      retval = polyVolume(n,p,poly,conn);

    delete[] p;
  }
  //  else volume of tri/quad is zero
  
  return retval;
}

// Area of a polygon in space.  The points must lie in a plane Dot(N,X) = c
// and be ordered within that plane (clockwise or counterclockwise).
template <class T>
T
GbMeasure<T>::polyArea(int n, const Point3* p)
{
  Point3 edge0 = { p[1].x-p[0].x, p[1].y-p[0].y, p[1].z-p[0].z };
  Point3 edge1 = { p[2].x-p[0].x, p[2].y-p[0].y, p[2].z-p[0].z };
  Point3 nor = {
    edge0.y*edge1.z - edge0.z*edge1.y,
    edge0.z*edge1.x - edge0.x*edge1.z,
    edge0.x*edge1.y - edge0.y*edge1.x
  };
  T invlen = T(1.0)/GbMath<T>::Sqrt(nor.x*nor.x+nor.y*nor.y+nor.z*nor.z);
  nor.x *= invlen;
  nor.y *= invlen;
  nor.z *= invlen;

  Point3 csum = { T(0.0), T(0.0), T(0.0) };

  for (int i = 0, j = 1; i <= n-2; i++, j++) {
    csum.x += p[i].y*p[j].z-p[i].z*p[j].y;
    csum.y += p[i].z*p[j].x-p[i].x*p[j].z;
    csum.z += p[i].x*p[j].y-p[i].y*p[j].x;
  }
  csum.x += p[n-1].y*p[0].z-p[n-1].z*p[0].y;
  csum.y += p[n-1].z*p[0].x-p[n-1].x*p[0].z;
  csum.z += p[n-1].x*p[0].y-p[n-1].y*p[0].x;

  return T(0.5)*GbMath<T>::Abs(nor.x*csum.x+nor.y*csum.y+nor.z*csum.z);
}

// Volume of a polyhedron in space.  The polyhederon must not have triangular
// faces.  The vertices of the each face must be ordered (clockwise or
// counterclockwise) in the plane containing the face.  The orderings for all
// the faces must be consistent (that is, all orderings are clockwise or all
// orderings are counterclockwise).
// Points are p[0] through p[numPoints-1].  The connectivity list has 'polys'
// blocks, each block representing a polygonal face.  A block starts with the
// number M of the vertices of that polygon.  The next M entries are indices
// into the point list (each entry is between 0 and n-1).  The vertices for
// the polygon must be ordered (clockwise or counterclockwise) in the plane
// of that polygon.  All orderings must be consistent (all polygons clockwise
// or all polygons counterclockwise).
template <class T>
T
GbMeasure<T>::polyVolume(int n, const Point3* p, int polys, const int* connectivity)
{
  static const T oneSixth = T(1.0)/T(6.0);
  T sum = T(0.0);

  for (int i = 0, j = 0; i < polys; i++) {
    // number of vertices for polygonal face i
    int numVerts = connectivity[j++];

    // indices of three vertices for polygon (assert: noncollinear verts)
    int k0 = connectivity[j];
    int k1 = connectivity[j+1];
    int k2 = connectivity[j+2];

    // compute normal for polygonal face
    Point3 edge0 = { p[k1].x-p[k0].x, p[k1].y-p[k0].y, p[k1].z-p[k0].z };
    Point3 edge1 = { p[k2].x-p[k0].x, p[k2].y-p[k0].y, p[k2].z-p[k0].z };
    Point3 nor = {
      edge0.y*edge1.z - edge0.z*edge1.y,
      edge0.z*edge1.x - edge0.x*edge1.z,
      edge0.x*edge1.y - edge0.y*edge1.x
    };
    T invlen = T(1.0)/GbMath<T>::Sqrt(nor.x*nor.x+nor.y*nor.y+nor.z*nor.z);
    nor.x *= invlen;
    nor.y *= invlen;
    nor.z *= invlen;

    // compute 2*area of polygonal face
    Point3 csum = { T(0.0), T(0.0), T(0.0) };

    for (int k = 0, j0 = j, j1 = j+1; k <= numVerts-2; k++, j0++, j1++) {
      k0 = connectivity[j0];
      k1 = connectivity[j1];
      csum.x += p[k0].y*p[k1].z-p[k0].z*p[k1].y;
      csum.y += p[k0].z*p[k1].x-p[k0].x*p[k1].z;
      csum.z += p[k0].x*p[k1].y-p[k0].y*p[k1].x;
    }
    k0 = connectivity[j+numVerts-1];
    k1 = connectivity[j];
    csum.x += p[k0].y*p[k1].z-p[k0].z*p[k1].y;
    csum.y += p[k0].z*p[k1].x-p[k0].x*p[k1].z;
    csum.z += p[k0].x*p[k1].y-p[k0].y*p[k1].x;

    T twoArea = GbMath<T>::Abs(nor.x*csum.x+nor.y*csum.y+nor.z*csum.z);
    T mult = GbMath<T>::Abs(p[k1].x*nor.x+p[k1].y*nor.y+p[k1].z*nor.z);
    T incr = mult*twoArea;
    sum += incr;

    j += numVerts;
  }

  return oneSixth*GbMath<T>::Abs(sum);
}

// This code is an implementation of an algorithm developed in
//
//   On gray scale image measurements:  II. Surface area and volume,
//   Dave Eberly, Jack Lancaster, and Majeid Alyassin,
//   CVGIP: Graphical Models and Image Processing, vol. 53, no. 6 (1991),
//   pp. 550-562.
//
// To use the method you need a point list which contains boundary
// points and corresponding normal vectors for a simply-connected object
// (no holes).  If you have an object with holes, you can always
// separate the disjoint surfaces making up the boundary and run the
// algorithm separately on these point lists (and add/subtract 
// appropriately).
//
// The voxel spacings are (dx,dy,dz) and are typically measured in
// millimeters or centimeters.  The values dx and dy are the spacings
// between pixels in a slice, and dz is the spacing between slices.
//
// Each voxel has location (x,y,z) in the lattice of voxels, where the
// coordinates are the array indices (integer-type).  To compute the
// corresponding normal vector (nx,ny,nz) you can use centralized
// differences:  nx = I(x+1,y,z)-I(x-1,y,z), ny = I(x,y+1,z)-I(x,y-1,z),
// and  nz = I(x,y,z+1)-I(x,y,z-1)  where I(x,y,z) are the gray scale
// voxel values.
template <class T>
void 
GbMeasure<T>::volume(int n, Point* p, T dx, T dy, T dz, T& surfaceArea, T& volume)
{
  const T epsilon = std::numeric_limits<T>::epsilon(); // 1e-06
  const T dxdy = dx*dy;
  const T dxdz = dx*dz;
  const T dydz = dy*dz;

  // volume elements
  T volumeX = T(0.0), volumeY = T(0.0), volumeZ = T(0.0);

  // initialize surface area
  surfaceArea = T(0.0);

  // bounding box of data points
  T xmax = -std::numeric_limits<T>::max(), xmin = std::numeric_limits<T>::max();
  T ymax = -std::numeric_limits<T>::max(), ymin = std::numeric_limits<T>::max();
  T zmax = -std::numeric_limits<T>::max(), zmin = std::numeric_limits<T>::max();
  for (int i = 0; i < n; i++) {
    // construct boundary point (bx,by,bz) and unit normal (nx,ny,nz)
    T bx = dx*T(p[i].x);
    T by = dy*T(p[i].y);
    T bz = dz*T(p[i].z);
    T nx = dx*T(p[i].nx);
    T ny = dy*T(p[i].ny);
    T nz = dz*T(p[i].nz);
    T length = GbMath<T>::Sqrt(nx*nx+ny*ny+nz*nz);
    if ( length > epsilon ) {
      T invLength = T(1.0)/length;
      nx *= invLength;
      ny *= invLength;
      nz *= invLength;
    }
    else {
      // Zero normal vector indicates constant density region, hard
      // to determine contribution to volume.
      continue;
    }

    // update bounding box of object
    if ( bx < xmin ) xmin = bx;
    else if ( bx > xmax ) xmax = bx;
    if ( by < ymin ) ymin = by;
    else if ( by > ymax ) ymax = by;
    if ( bz < zmin ) zmin = bz;
    else if ( bz > zmax ) zmax = bz;

    // compute surface area contribution
    T ax = GbMath<T>::Abs(nx), ay = GbMath<T>::Abs(ny), az = GbMath<T>::Abs(nz);
    T dS;
    if ( ax > ay ) {
      if ( ax > az ) {
	// |nx| is maximum
	dS = dydz/ax;
      }
      else {
	// |nz| is maximum
	dS = dxdy/az;
      }
    }
    else {
      if ( ay > az ) {
	// |ny| is maximum
	dS = dxdz/ay;
      }
      else {
	// |nz| is maximum
	dS = dxdy/az;
      }
    }

    // update running sums for surface area and volume
    surfaceArea += dS;
    volumeX += bx*nx*dS;
    volumeY += by*ny*dS;
    volumeZ += bz*nz*dS;
  }

  // extent of object
  T ex = xmax-xmin, ey = ymax-ymin, ez = zmax-zmin;

  // shape parameters
  T invDenom = T(1.0)/(ex*ey+ex*ez+ey*ez);
  T sx = ey*ez*invDenom;
  T sy = ex*ez*invDenom;
  T sz = ex*ey*invDenom;

  volume = sx*volumeX + sy*volumeY + sz*volumeZ;

  // adjust for inwards pointing normals
  if ( volume < T(0.0) )
    volume = -volume;
}
