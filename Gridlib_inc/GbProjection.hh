/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBPROJECTION_HH
#define  GBPROJECTION_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>
#include <vector>

#include "GbMath.hh"
#include "GbTypes.hh"
#include "GbVec2.hh"
#include "GbMatrix2.hh"
#include "GbVec3.hh"
#include "GbPlane.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbProjection
{
public:

  /*!
    Homogeneous mapping of quadrilateral <p00,p10,p11,p01> to square [0,1]^2.
    The quadrilateral points are ordered counterclockwise and map onto the
    corners (0,0), (1,0), (1,1), and (0,1), respectively.
  */
  class HmQuadToSquare
  {
  public:
    HmQuadToSquare (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
		    const GbVec2<T>& rkP11, const GbVec2<T>& rkP01);
    
    GbVec2<T> transform (const GbVec2<T>& rkP);
    
  protected:
    GbVec2<T> kT_, kG_, kD_;
    GbMatrix2<T> kM_;
  };

  /*!
    Homogeneous mapping of square [0,1]^2 to quadrilateral <p00,p10,p11,p01>.
    The quadrilateral points are ordered counterclockwise and map onto the
    corners (0,0), (1,0), (1,1), and (0,1), respectively.
  */
  class HmSquareToQuad
  {
  public:
    HmSquareToQuad (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
		    const GbVec2<T>& rkP11, const GbVec2<T>& rkP01);

    GbVec2<T> transform (const GbVec2<T>& rkP);
    
  protected:
    GbVec2<T> kT_, kG_, kD_;
    GbMatrix2<T> kM_;
  };

  /*!
    Bilinear mapping of quadrilateral <p00,p10,p11,p01> to square [0,1]^2.
    The quadrilateral points are ordered counterclockwise and map onto the
    corners (0,0), (1,0), (1,1), and (0,1), respectively.
  */
  class BiQuadToSquare
  {
  public:
    BiQuadToSquare (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
		    const GbVec2<T>& rkP11, const GbVec2<T>& rkP01);

    GbVec2<T> transform (const GbVec2<T>& rkP);
    
  protected:
    GbVec2<T> kA_, kB_, kC_, kD_;
    T BCdet_, CDdet_;
  };

  /*!
    Bilinear mapping of square [0,1]^2 to quadrilateral <p00,p10,p11,p01>.
    The quadrilateral points are ordered counterclockwise and map onto the
    corners (0,0), (1,0), (1,1), and (0,1), respectively.
  */
  class BiSquareToQuad
  {
  public:
    BiSquareToQuad (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
		    const GbVec2<T>& rkP11, const GbVec2<T>& rkP01);

    GbVec2<T> transform (const GbVec2<T>& rkP);
    
  protected:
    GbVec2<T> S00_, S01_, S10_, S11_;
  };

  /*!
    Compute terminator and silhouette of convex polyhedron.  The "silhouette"
    is the convex hull of the planar projection of the polyhedron vertices.
    The "terminator" is the set of polyhedron vertices that map to the
    silhouette.
  */
  class ConvexPolyhedron
  {
  public:
    // forward declares for vertex-edge-triangle table
    class Edge;
    class Triangle;

    class Vertex
    {
    public:
      // The application is responsible for deleting the Vertex points if
      // they were dynamically allocated.  Of course, it is possible the
      // pointer is into a vertex pool that persists across application
      // life time, so ConvexPolyhedron is not going to make any
      // assumptions about this.
      GbVec3<T>* point_;
      int nbrEdges_;
      Edge** edge_;
    };

    class Edge
    {
    public:
      Vertex* vertex_[2];
      Triangle* triangle_[2];
    };

    class Triangle
    {
    public:
      //! plane of triangle (outward pointing normal)
      GbPlane<T> plane_;

      //! pointers to triangle vertices (counterclockwise order)
      Vertex* vertex_[3];

      //! triangle edges: E0 = <V0,V1>, E1 = <V1,V2>, E2 = <V2,V0>
      Edge* edge_[3];

      //! adjacent triangles:  adjacent[i] shares edge[i]
      Triangle* adjacent_[3];

      //! Temporary storage for signed distance from eye point to plane of
      //! triangle (initialized to MAXREAL).
      T distance_;
    };

    /*!
      Construction and destruction.  ConvexPolyhedron accepts
      responsibility for deleting the input arrays, but not for deleting the
      Vertex points if they were dynamically allocated.  The application is
      responsible for the ensuring that the polyhedron is convex.
    */
    ConvexPolyhedron (int iVertexQuantity, Vertex* akVertex,
		      int iEdgeQuantity, Edge* akEdge, int iTriangleQuantity,
		      Triangle* akTriangle);

    ~ConvexPolyhedron ();

    INLINE int getVertexQuantity () const;
    INLINE Vertex& getVertex (int i);
    INLINE Vertex* getVertices ();

    INLINE int getEdgeQuantity () const;
    INLINE Edge& getEdge (int i);
    INLINE Edge* getEdges ();

    INLINE int getTriangleQuantity () const;
    INLINE Triangle& getTriangle (int i);
    INLINE Triangle* getTriangles ();

    /*!
      The eye point must be outside the polyhedron.  The output is the
      terminator, an ordered list of vertices forming a simple closed
      polyline that separates the visible from invisible faces of the
      polyhedron.
    */
    void computeTerminator (const GbVec3<T>& rkEye, std::vector<GbVec3<T>*>& rkTerminator);

    /*!
      If projection plane is Dot(N,X) = c where N is unit length, then the
      application must ensure that Dot(N,eye) > c.  That is, the eye point is
      on the side of the plane to which N points.  The application must also
      specify two vectors U and V in the projection plane so that {U,V,N} is
      a right-handed and orthonormal set (the matrix [U V N] is orthonormal
      with determinant 1).  The origin of the plane is computed internally as
      the closest point to the eye point (an orthogonal pyramid for the
      perspective projection).  If all vertices P on the terminator satisfy
      Dot(N,P) < Dot(N,eye), then the polyhedron is completely visible (in
      the sense of perspective projection onto the viewing plane).  In this
      case the silhouette is computed by projecting the terminator points
      onto the viewing plane.  The return value of the function is 'true'
      when this happens.  However, if at least one terminator point P
      satisfies Dot(N,P) >= Dot(N,eye), then the silhouette is unbounded in
      the view plane.  It is not computed and the function returns 'false'.
      A silhouette point (x,y) is extracted from the point Q that is the
      intersection of the ray whose origin is the eye point and that contains
      a terminator point, Q = K+x*U+y*V+z*N where K is the origin of the
      plane.
    */
    GbBool computeSilhouette (const GbVec3<T>& rkEye, const GbPlane<T>& rkPlane,
			      const GbVec3<T>& rkU, const GbVec3<T>& rkV,
			      std::vector<GbVec2<T> >& rkSilhouette);

    GbBool computeSilhouette (std::vector<GbVec3<T>*>& rkTerminator,
			      const GbVec3<T>& rkEye, const GbPlane<T>& rkPlane, const GbVec3<T>& rkU,
			      const GbVec3<T>& rkV, std::vector<GbVec2<T> >& rkSilhouette);

  protected:
    T getDistance (const GbVec3<T>& rkEye, Triangle* pkTriangle);
    GbBool isNegativeProduct (T fDist0, T fDist1);

    int nbrVertices_;
    Vertex* vertex_;

    int nbrEdges_;
    Edge* edge_;

    int nbrTriangles_;
    Triangle* triangle_;
  };

protected:
  static const T TOLERANCE;

};


#ifdef __GNUG__

#include "GbProjection.in"
#include "GbProjection.T"

#else

#pragma instantiate GbProjection<float>
#pragma instantiate GbProjection<double>

#ifndef OUTLINE
#include "GbProjection.in"
#endif

#endif  // g++

#endif  // GBPROJECTION_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.1  2001/06/26 08:31:24  prkipfer
| added 2D Projections
|
|
+---------------------------------------------------------------------*/
