// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/06/26 08:31:24  prkipfer
| added 2D Projections
|
|
|
+---------------------------------------------------------------------*/

template <class T>
const T
GbProjection<T>::TOLERANCE = std::numeric_limits<T>::epsilon(); // max. 1e-6

template <class T>
GbProjection<T>::HmQuadToSquare::HmQuadToSquare (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
						 const GbVec2<T>& rkP11, const GbVec2<T>& rkP01)
{
  // translate to origin
  kT_ = rkP00;
  GbVec2<T> kQ10 = rkP10 - rkP00;
  GbVec2<T> kQ11 = rkP11 - rkP00;
  GbVec2<T> kQ01 = rkP01 - rkP00;

  GbMatrix2<T> kInvM(kQ10[0],kQ01[0],kQ10[1],kQ01[1]);
  GbBool bInverted = kInvM.inverse(kM_,GbProjection<T>::TOLERANCE);
  assert( bInverted );  // else quad degenerate

  // compute where p11 is mapped to
  GbVec2<T> kCorner = kM_*kQ11;  // = (a,b)

  // Compute homogeneous transform of quadrilateral
  // {(0,0),(1,0),(a,b),(0,1)} to square {(0,0),(1,0),(1,1),(0,1)}
  kG_[0] = (kCorner[1] - 1.0)/kCorner[0];
  kG_[1] = (kCorner[0] - 1.0)/kCorner[1];
  kD_[0] = 1.0 + kG_[0];
  kD_[1] = 1.0 + kG_[1];
}

template <class T>
GbVec2<T> 
GbProjection<T>::HmQuadToSquare::transform (const GbVec2<T>& rkP)
{
  GbVec2<T> kProd = kM_*(rkP - kT_);
  T fInvDenom = 1.0/(1.0 + kG_.dot(kProd));
  GbVec2<T> kResult = fInvDenom*kProd;
  kResult[0] *= kD_[0];
  kResult[1] *= kD_[1];
  return kResult;
}

template <class T>
GbProjection<T>::HmSquareToQuad::HmSquareToQuad (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
						 const GbVec2<T>& rkP11, const GbVec2<T>& rkP01)
{
  // translate to origin
  kT_ = rkP00;
  kM_[0][0] = rkP10[0] - rkP00[0];
  kM_[0][1] = rkP01[0] - rkP00[0];
  kM_[1][0] = rkP10[1] - rkP00[1];
  kM_[1][1] = rkP01[1] - rkP00[1];

  GbMatrix2<T> kInvM;
  GbBool bInverted = kM_.inverse(kInvM,GbProjection<T>::TOLERANCE);
  assert( bInverted );  // else quad degenerate

  // find point which is mapped to p11
  GbVec2<T> kCorner = kInvM*(rkP11-rkP00);  // = (a,b)

  // compute homogeneous transform of square {(0,0),(1,0),(1,1),(0,1)} to
  // quadrilateral {(0,0),(1,0),(a,b),(0,1)}
  T fInvDenom = 1.0/(kCorner[0] + kCorner[1] - 1.0);
  kG_[0] = fInvDenom*(1.0 - kCorner[1]);
  kG_[1] = fInvDenom*(1.0 - kCorner[0]);
  kD_[0] = fInvDenom*kCorner[0];
  kD_[1] = fInvDenom*kCorner[1];
}

template <class T>
GbVec2<T> 
GbProjection<T>::HmSquareToQuad::transform (const GbVec2<T>& rkP)
{
  T fInvDenom = 1.0/(1.0 + kG_.dot(rkP));
  GbVec2<T> kResult(kD_[0]*rkP[0],kD_[1]*rkP[1]);
  GbVec2<T> kProd = kM_*kResult;
  kResult[0] = fInvDenom*kProd[0] + kT_[0];
  kResult[1] = fInvDenom*kProd[1] + kT_[1];
  return kResult;
}

template <class T>
GbProjection<T>::BiQuadToSquare::BiQuadToSquare (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
						 const GbVec2<T>& rkP11, const GbVec2<T>& rkP01)
{
  kA_ = rkP00;
  kB_ = rkP10 - rkP00;
  kC_ = rkP01 - rkP00;
  kD_ = rkP11 + rkP00 - rkP10 - rkP01;
  BCdet_ = kB_[0]*kC_[1] - kB_[1]*kC_[0];
  assert( GbMath<T>::Abs(BCdet_) >= GbProjection<T>::TOLERANCE );  // else quad degenerate
  CDdet_ = kC_[1]*kD_[0] - kC_[0]*kD_[1];
}

template <class T>
GbVec2<T> 
GbProjection<T>::BiQuadToSquare::transform (const GbVec2<T>& rkP)
{
  GbVec2<T> kDiff = kA_ - rkP;
  T fABdet = kDiff[1]*kB_[0]-kDiff[0]*kB_[1];
  T fADdet = kDiff[1]*kD_[0]-kDiff[0]*kD_[1];

  T fA = CDdet_;
  T fB = fADdet + BCdet_;
  T fC = fABdet;

  GbVec2<T> kResult;

  if ( GbMath<T>::Abs(fA) >= GbProjection<T>::TOLERANCE ) {
    // t-equation is quadratic
    T fDiscr = GbMath<T>::Sqrt(GbMath<T>::Abs(fB*fB-4.0*fA*fC));
    kResult[0] = (-fB + fDiscr)/(2.0*fA);
    if ( kResult[0] < 0.0 || kResult[0] > 1.0 ) {
      kResult[0] = (-fB - fDiscr)/(2.0*fA);
      if ( kResult[0] < 0.0 || kResult[0] > 1.0 ) {
	// point p not inside quadrilateral, return invalid result
	return GbVec2<T>(-1.0,-1.0);
      }
    }
  }
  else {
    // t-equation is linear
    kResult[0] = -fC/fB;
  }

  kResult[1] = -(kDiff[0] + kResult[0]*kC_[0])/(kB_[0] + kResult[0]*kD_[0]);

  return kResult;
}

template <class T>
GbProjection<T>::BiSquareToQuad::BiSquareToQuad (const GbVec2<T>& rkP00, const GbVec2<T>& rkP10,
						 const GbVec2<T>& rkP11, const GbVec2<T>& rkP01)
{
  S00_ = rkP00;
  S10_ = rkP10;
  S11_ = rkP11;
  S01_ = rkP01;
}

template <class T>
GbVec2<T> 
GbProjection<T>::BiSquareToQuad::transform (const GbVec2<T>& rkP)
{
  GbVec2<T> kOmP(1.0-rkP[0],1.0-rkP[1]);
  GbVec2<T> kResult;
  kResult[0] = kOmP[1]*(kOmP[0]*S00_[0] + rkP[0]*S10_[0]) +
               rkP[1] *(kOmP[0]*S01_[0] + rkP[0]*S11_[0]);
  kResult[1] = kOmP[1]*(kOmP[0]*S00_[1] + rkP[0]*S10_[1]) +
               rkP[1] *(kOmP[0]*S01_[1] + rkP[0]*S11_[1]);
  return kResult;
}

template <class T>
GbProjection<T>::ConvexPolyhedron::ConvexPolyhedron (int iVertexQuantity, Vertex* akVertex,
						     int iEdgeQuantity, Edge* akEdge, int iTriangleQuantity,
						     Triangle* akTriangle)
{
  nbrVertices_ = iVertexQuantity;
  vertex_ = akVertex;
  nbrEdges_ = iEdgeQuantity;
  edge_ = akEdge;
  nbrTriangles_ = iTriangleQuantity;
  triangle_ = akTriangle;
}

template <class T>
GbProjection<T>::ConvexPolyhedron::~ConvexPolyhedron ()
{
  for (int i = 0; i < nbrVertices_; i++)
    delete[] vertex_[i].edge_;

  delete[] vertex_;
  delete[] edge_;
  delete[] triangle_;
}

template <class T>
T 
GbProjection<T>::ConvexPolyhedron::getDistance (const GbVec3<T>& rkEye, Triangle* pkTriangle)
{
  // Signed distance from eye to plane of triangle.  When distance is
  // positive, triangle is visible from eye (front-facing).  When distance
  // is negative, triangle is not visible from eye (back-facing).  When
  // distance is zero, triangle is visible "on-edge" from eye.

  if ( pkTriangle->distance_ == std::numeric_limits<T>::max() ) {
    pkTriangle->distance_ = pkTriangle->plane_.distanceTo(rkEye);

    if ( GbMath<T>::Abs(pkTriangle->distance_) < GbProjection<T>::TOLERANCE )
      pkTriangle->distance_ = 0.0;
  }

  return pkTriangle->distance_;
}

template <class T>
GbBool 
GbProjection<T>::ConvexPolyhedron::isNegativeProduct (T fDist0, T fDist1)
{
  return ( fDist0 != 0.0 ? (fDist0*fDist1 <= 0.0) : (fDist1 != 0.0) );
}

template <class T>
void 
GbProjection<T>::ConvexPolyhedron::computeTerminator (const GbVec3<T>& rkEye,
						      std::vector<GbVec3<T>*>& rkTerminator)
{
  int i;
  for (i = 0; i < nbrTriangles_; i++)
    triangle_[i].distance_ = std::numeric_limits<T>::max();

  // Start a search for a front-facing triangle that has an adjacent
  // back-facing triangle or for a back-facing triangle that has an
  // adjacent front-facing triangle.
  Triangle* pkTriangle = triangle_;
  T fTriDist = getDistance(rkEye,pkTriangle);
  Edge* pkFirstEdge = 0;
  int iIter;
  for (iIter = 0; iIter < nbrTriangles_; iIter++) {
    // Check adjacent neighbors for edge of terminator.  Such an
    // edge occurs if the signed distance changes sign.
    int iMinIndex = -1;
    T fMinAbsDist = std::numeric_limits<T>::max();
    T afAdjDist[3];
    for (i = 0; i < 3; i++) {
      afAdjDist[i] = getDistance(rkEye,pkTriangle->adjacent_[i]);
      if ( isNegativeProduct(fTriDist,afAdjDist[i]) ) {
	pkFirstEdge = pkTriangle->edge_[i];
	break;
      }

      T fAbsDist = GbMath<T>::Abs(afAdjDist[i]);
      if ( fAbsDist < fMinAbsDist ) {
	fMinAbsDist = fAbsDist;
	iMinIndex = i;
      }
    }
    if ( i < 3 )
      break;

    // First edge not found during this iteration.  Move to adjacent
    // triangle whose distance is smallest of all adjacent triangles.
    pkTriangle = pkTriangle->adjacent_[iMinIndex];
    fTriDist = afAdjDist[iMinIndex];
  }
  assert( iIter < nbrTriangles_ );

  rkTerminator.push_back(pkFirstEdge->vertex_[0]->point_);
  rkTerminator.push_back(pkFirstEdge->vertex_[1]->point_);

  // walk along the terminator
  Vertex* pkFirstVertex = pkFirstEdge->vertex_[0];
  Vertex* pkVertex = pkFirstEdge->vertex_[1];
  Edge* pkEdge = pkFirstEdge;
  for (iIter = 0; iIter < nbrEdges_; iIter++) {
    // search all edges sharing the vertex for another terminator edge
    for (i = 0; i < pkVertex->nbrEdges_; i++) {
      Edge* pkNextEdge = pkVertex->edge_[i];
      if ( pkNextEdge == pkEdge )
	continue;

      T fDist0 = getDistance(rkEye,pkNextEdge->triangle_[0]);
      T fDist1 = getDistance(rkEye,pkNextEdge->triangle_[1]);
      if ( isNegativeProduct(fDist0,fDist1) ) {
	if ( pkNextEdge->vertex_[0] == pkVertex ) {
	  pkVertex = pkNextEdge->vertex_[1];
	  rkTerminator.push_back(pkVertex->point_);
	  if ( pkVertex == pkFirstVertex )
	    return;
	}
	else {
	  pkVertex = pkNextEdge->vertex_[0];
	  rkTerminator.push_back(pkVertex->point_);
	  if ( pkVertex == pkFirstVertex )
	    return;
	}

	pkEdge = pkNextEdge;
	break;
      }
    }
    assert( i < pkVertex->nbrEdges_ );
  }
  assert( iIter < nbrEdges_ );
}

template <class T>
GbBool 
GbProjection<T>::ConvexPolyhedron::computeSilhouette (const GbVec3<T>& rkEye,
						      const GbPlane<T>& rkPlane, const GbVec3<T>& rkU, const GbVec3<T>& rkV,
						      std::vector<GbVec2<T> >& rkSilhouette)
{
  std::vector<GbVec3<T>*> kTerminator;
  computeTerminator(rkEye,kTerminator);
  return computeSilhouette(kTerminator,rkEye,rkPlane,rkU,rkV,rkSilhouette);
}

template <class T>
GbBool 
GbProjection<T>::ConvexPolyhedron::computeSilhouette (std::vector<GbVec3<T>*>& rkTerminator,
						      const GbVec3<T>& rkEye, const GbPlane<T>& rkPlane, const GbVec3<T>& rkU,
						      const GbVec3<T>& rkV, std::vector<GbVec2<T> >& rkSilhouette)
{
  T fEDist = rkPlane.distanceTo(rkEye);  // assert:  fEDist > 0

  // closest planar point to E is K = E-dist*N
  GbVec3<T> kClosest = rkEye - fEDist*rkPlane.getNormal();

  // project polyhedron points onto plane
  for (int i = 0; i < rkTerminator.size(); i++) {
    GbVec3<T>* pkPoint = rkTerminator[i];

    T fVDist = rkPlane.distanceTo(*pkPoint);
    if ( fVDist >= fEDist ) {
      // cannot project vertex onto plane
      return false;
    }

    // compute projected point Q
    T fRatio = fEDist/(fEDist-fVDist);
    GbVec3<T> kProjected(rkEye[0] + fRatio*((*pkPoint)[0] - rkEye[0]),
			 rkEye[1] + fRatio*((*pkPoint)[1] - rkEye[1]),
			 rkEye[2] + fRatio*((*pkPoint)[2] - rkEye[2]));

    // compute (x,y) so that Q = K+x*U+y*V+z*N
    GbVec3<T> kDiff = kProjected - kClosest;
    rkSilhouette.push_back(GbVec2<T>(rkU.dot(kDiff),rkV.dot(kDiff)));
  }

  return true;
}
