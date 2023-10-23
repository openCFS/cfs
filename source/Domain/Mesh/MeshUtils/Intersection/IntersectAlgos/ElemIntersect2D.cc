// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PlaneIntersect.cc
 *       \brief    This code is mainly taken from the Grid class.
 *                 TODO replace intersections in Grid class
 *                 With this to avoid the doubleing of code...
 *
 *       \date     Apr 13, 2016
 *       \author   ahueppe
 */
//================================================================================================

#include "ElemIntersect2D.hh"

#ifdef USE_CGAL

#pragma GCC diagnostic push
// for gcc7 and CGAL 4.9.1
#pragma GCC diagnostic ignored "-Wuninitialized"

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2 CGALPoint2;
typedef CGAL::Polygon_2<Kernel> CGALPolygon2;
typedef CGAL::Polygon_with_holes_2<Kernel> CGALPolygonWithHoles2;

#pragma GCC diagnostic pop

#endif

namespace CoupledField{


void ElemIntersect2D::SetTElem(UInt tNum){
  lastIntersection_.Clear(true);

  this->tElemNum_ = tNum;
}

bool ElemIntersect2D::Intersect(UInt sNum){
  sElemNum_ = sNum;
  const Elem * trgElem = tGrid_->GetElem(tElemNum_);
  const Elem * srcElem = sGrid_->GetElem(sElemNum_);
  lastVolume_ = IntersectPolys( trgElem , srcElem);

  if(lastVolume_> 0){
    return true;
  }else{
    return false;
  }
}

void ElemIntersect2D::GetIntersectionElems(StdVector<IntersectionElem*>& interElems){
  EXCEPTION("Not Implemented")
  //create elements from polygon, see mortrar interface for details.
}

void ElemIntersect2D::GetVolumeAndCenters(StdVector<ElemIntersect::VolCenterInfo>& infos){
  infos.Resize(1);
  infos[0].volume = lastVolume_;
  //in case of 2D grids we compute 3D intersections and just take the
  //the first two coordinates
  if(is2dTrg_){
    Vector<Double> tmpCenter;
    PolyCentroid(lastIntersection_,tmpCenter);
    infos[0].center.Resize(2);
    infos[0].center[0] = tmpCenter[0];
    infos[0].center[1] = tmpCenter[1];
  }else{
    PolyCentroid(lastIntersection_,infos[0].center);
  }
  infos[0].targetElemNum = tElemNum_;
  infos[0].sourceElemNum = sElemNum_;
}



Double ElemIntersect2D::IntersectPolys( const Elem * trgElem , const Elem *srcElem){
  //Generate two polyhedra from the elements
  tPoly_.Resize(trgElem->connect.GetSize());
  for(UInt i=0;i<trgElem->connect.GetSize(); ++i){
    tGrid_->GetNodeCoordinate(tPoly_[i],trgElem->connect[i],false);
  }
  sPoly_.Resize(srcElem->connect.GetSize());
  for(UInt i=0;i<srcElem->connect.GetSize(); ++i){
    sGrid_->GetNodeCoordinate(sPoly_[i],srcElem->connect[i],false);
  }
  if(sPoly_[0].GetSize()==2){
    for(UInt i=0;i<srcElem->connect.GetSize(); ++i){
      sPoly_[i].Push_back(0.0);
    }
  }
  if(tPoly_[0].GetSize()==2){
    is2dTrg_ = true;
    for(UInt i=0;i<trgElem->connect.GetSize(); ++i){
      tPoly_[i].Push_back(0.0);
    }
  }else{
    is2dTrg_ = false;
  }
  return CutPolysCGAL(tPoly_,sPoly_,false,lastIntersection_);
}

#ifdef USE_CGAL
Double ElemIntersect2D::CutPolysCGAL(StdVector<Vector<Double> > &p1,
                                   StdVector<Vector<Double> > &p2, const bool coplanar,
                                   StdVector<Vector<Double> > &r) {
    UInt i;
   Vector<Double>& temp1 = temp1_;
   Vector<Double>& temp2 = temp2_;

   temp1.Resize(3);
   temp2.Resize(3);
 #ifdef CHECK_INDEX
   // check that we have actually polygons
   if ((p1.GetSize() < 3) || (p2.GetSize() < 3)) {
     EXCEPTION("A polygon must consist of 3 points at least");
     return 0;
   }
 #endif

   // compute surface normal of p2
   Vector<Double>& n = e_;
   temp1 = p2[1] - p2[0];
   temp2 = p2[2] - p2[0];
   temp1.CrossProduct(temp2, n);
   n.Normalize();

   // if interface is not coplanar then project p1 onto p2
   Double scale;

   // project each point of p1
   for (i = 0; i < p1.GetSize(); ++i) {
     temp1 = p1[i] - p2[0];
     scale = n.Inner(temp1);
     p1[i] -= n * scale;
   }


   // Now both polygons have the same normal and we can proceed
   // with the rotation of the polygons, in order to make them
   // lying parallel to the XY-plane. After that, CGAL 2D procedures can
   // be used to calculate intersections of the polygons.
   Vector<Double>& ez = c1_;
   Vector<Double>& v = c2_;
   ez.Resize(3);
   v.Resize(3);
   ez[2] = 1.0;
   n.CrossProduct(ez, v);
   //quick fix, if v has zero length, i.e. ez == n we set continue, setting v to zero
   if(v.NormL2()!=0)
         v.Normalize();

   Double ca = n[2]; // cos(n ^ ez) = n_z/|n| = n_z
   Double sa = sqrt(n[0] * n[0] + n[1] * n[1]); // sin(n ^ ez) = sqrt(n_x^2 + n_y^2)/|n| = sqrt(n_x^2 + n_y^2)
   Double ci = 1 - ca; // 1 - cos(n ^ ez)
   // rotation matrix
   Matrix<Double>& rMat = rMat_;
   rMat.Resize(3, 3);
   rMat[0][0] = ca + v[0] * v[0] * ci;
   rMat[1][0] = v[0] * v[1] * ci + v[2] * sa;
   rMat[2][0] = v[0] * v[2] * ci - v[1] * sa;
   rMat[0][1] = v[0] * v[1] * ci - v[2] * sa;
   rMat[1][1] = ca + v[1] * v[1] * ci;
   rMat[2][1] = v[1] * v[2] * ci + v[0] * sa;
   rMat[0][2] = v[0] * v[2] * ci + v[1] * sa;
   rMat[1][2] = v[1] * v[2] * ci - v[0] * sa;
   rMat[2][2] = ca + v[2] * v[2] * ci;
   // Perform rotations and create CGAL 2D polygons.
   // We do not change the initial polygons, because, for
   // example, p1 is used further to create a projected master element.
   StdVector<Vector<Double> >& p1Rot = p1Rot_;
   StdVector<Vector<Double> >& p2Rot = p2Rot_;
   p1Rot.Resize(p1.GetSize());
   p2Rot.Resize(p2.GetSize());
   CGALPolygon2 plgn1, plgn2;
   for (i = 0; i < p1.GetSize(); i++) {
     p1Rot[i] = rMat * p1[i];
     plgn1.push_back(CGALPoint2(p1Rot[i][0], p1Rot[i][1]));
   }
   for (i = 0; i < p2.GetSize(); i++) {
     p2Rot[i] = rMat * p2[i];
     plgn2.push_back(CGALPoint2(p2Rot[i][0], p2Rot[i][1]));
   }

   // store the z-coordinate of the rotated polygons
   Double zCoord = p1Rot[0][2];

   if (plgn1.is_clockwise_oriented())
     plgn1.reverse_orientation();
   if (plgn2.is_clockwise_oriented())
     plgn2.reverse_orientation();

   std::vector<CGALPolygonWithHoles2> intrsLst;
   CGAL::intersection(plgn1, plgn2, std::back_inserter(intrsLst));

   // If the intersection list has the size = 0, the polygons don't
   // intersect with each other. If the size is > 1, then something
   // must be wrong, because the intersection two convex sets is a convex set.
   if (intrsLst.size() != 1)
     return 0;

   // If the relative area of the resulting polygon is too little,
   // we omit this intersection.
   if ((intrsLst[0].outer_boundary().area() < 1.0e-5 * plgn1.area())
       || (intrsLst[0].outer_boundary().area() < 1.0e-5 * plgn2.area()))
     return 0;

   // We need to transform the vertices of the CGAL-polygon, gained
   // from the intersection procedure, back to a CFS-Vector. However,
   // sometimes the polygon can contain a negligibly short edge which
   // must be omitted. Therefore, we iterate over the edges, check their
   // lengths, and take the starting points of those having considerable
   // lengths. The reason why such intersections take place is still unknown.
   Matrix<Double>& rMatTrans = rMatTrans_;
   rMat.Transpose(rMatTrans);

   for (CGALPolygon2::Edge_const_iterator edgeIt =
       intrsLst[0].outer_boundary().edges_begin();
       edgeIt != intrsLst[0].outer_boundary().edges_end(); ++edgeIt) {
     if (edgeIt->squared_length() < 1.0e-5 * intrsLst[0].outer_boundary().area())
       continue;

     // temp1 - a vertex of the intersection polygon; temp2 - a vertex rotated back
     temp1[0] = CGAL::to_double(edgeIt->vertex(0).x());
     temp1[1] = CGAL::to_double(edgeIt->vertex(0).y());
     // The polygon is parallel to the Oxy-plane, so we can choose any z-coordinate
     // for its vertices. The only condition - it must be the same for all vertices.
     temp1[2] = zCoord;
     // Perform the back-rotation, so that the polygon be parallel to the master/slave plane.
     // It is required in order to get correct results transforming Local-to-Global and back
     // within SurfaceMortarABInt::CalcElementMatrix method.
     temp2 = rMatTrans * temp1;
     r.Push_back(temp2);
   }
   //assume volume preserving projection
   return CGAL::to_double(intrsLst[0].outer_boundary().area());
 }
#else
Double ElemIntersect2D::CutPolysCGAL(StdVector<Vector<Double> > &p1,
                                   StdVector<Vector<Double> > &p2, const bool coplanar,
                                   StdVector<Vector<Double> > &r) {
  EXCEPTION("CGAL SUPPORT is needed for 2D element intersetion.")
  return 0;
}
#endif


void ElemIntersect2D::PolyCentroid(const StdVector< Vector<Double> > &p,
                                  Vector<Double> &c)
{
  UInt i, n = p.GetSize();
  // set c to 0
  c.Resize(3);
  c.Init(0.0);

  // compute center of gravity
  for (i = 0; i < n; ++i){
    c[0] += p[i][0];
    c[1] += p[i][1];
    c[2] += p[i][2];
  }
  for(i=0;i<3;i++)
    c[i] /= (Double) n;
}

/*bool ElemIntersect2D::IntersectRects( const Elem * trgElem , const Elem *srcElem){
  Vector<Double> c0, c1, c2, c3, d0, d1, d2, d3;
  Vector<Double> diffS, diffX, diffY, diffX2;
  Vector<Double> s, t;
  StdVector<UInt> connect2;
  StdVector<Vector<Double> > coord2;
  Double distX, distY, distX2, facX, facY, r;
  UInt nodeNr;
  // Introduce a tolerance to account for roundoff errors during the calculation of
  // normed new x base vector.
  Double tol_r;

  s.Resize(4);
  t.Resize(4);
  connect2.Resize(4);
  coord2.Resize(4);

  // The meaning of the points c0, c1, c3, d0 and d2
  // is as follows:
  //                x------------------x d2
  //                |                  |
  //                |                  |
  // c3 x-----------+----------x       |
  //    |           |          |       |
  //    |           |          |       |
  //    |        d0 x----------+-------x d1
  //    |                      |
  //    |                      |
  // c0 x----------------------x c1


  // Get coordinates of the endpoints
  tGrid_->GetNodeCoordinate(c0, trgElem->connect[0], false);
  tGrid_->GetNodeCoordinate(c1, trgElem->connect[1], false);
  tGrid_->GetNodeCoordinate(c2, trgElem->connect[2], false);
  tGrid_->GetNodeCoordinate(c3, trgElem->connect[3], false);
  sGrid_->GetNodeCoordinate(d0, srcElem->connect[0], false);
  sGrid_->GetNodeCoordinate(d1, srcElem->connect[1], false);
  sGrid_->GetNodeCoordinate(d2, srcElem->connect[2], false);
  sGrid_->GetNodeCoordinate(d3, srcElem->connect[3], false);

  // Compute and normalize vector from c0 to c1.
  // This becomes the new x-unit vector.
  diffX = c1 - c0;
  distX = diffX.NormL2();
  facX = 1.0 / distX;
  diffX *= facX;

  // Compute and normalize vector from c0 to c3
  // This becomes the new y-unit vector.
  diffY = c3 - c0;
  distY = diffY.NormL2();
  facY = 1.0 / distY;
  diffY *= facY;

  // Now compute vector from c0 to d0 and project
  // the result onto the new x- and y-axis.
  diffS = d0 - c0;
  diffS.Inner(diffX, s[0]);
  diffS.Inner(diffY, s[1]);
  s[0] *= facX;
  s[1] *= facY;

  // Now compute vector from c0 to d2 and project
  // the result onto the new x- and y-axis.
  diffS = d2 - c0;
  diffS.Inner(diffX, s[2]);
  diffS.Inner(diffY, s[3]);
  s[2] *= facX;
  s[3] *= facY;

  // Determine the orientation of the second rectangle
  // to make sure that the edges which connect c0 and c1
  // are parallel to the edges which connect d0 and d1.
  diffX2 = d1 - d0;
  distX2 = diffX2.NormL2();
  diffX.Inner(diffX2, r);

  // Set the tolerance for determining if the edges
  // mentioned in the last comment are parallel.
  tol_r = distX2 < distX ? distX2 / 10 : distX / 10;

  // Bring the x- and y-coordinates of the intersection
  // into an order, where the smaller coordinates come
  // first.
  if(s[2] < s[0])
  {
    t[0] = s[2];
    t[2] = s[0];

    if(s[3] < s[1])
    {
      t[1] = s[3];
      t[3] = s[1];
      connect2[0] = srcElem->connect[2];
      connect2[2] = srcElem->connect[0];
      coord2[0] = d2;
      coord2[2] = d0;
      if (fabs(r) < tol_r) {
        connect2[1] = srcElem->connect[1];
        connect2[3] = srcElem->connect[3];
        coord2[1] = d1;
        coord2[3] = d3;
      }
      else {
        connect2[1] = srcElem->connect[3];
        connect2[3] = srcElem->connect[1];
        coord2[1] = d3;
        coord2[3] = d1;
      }
    }
    else
    {
      t[1] = s[1];
      t[3] = s[3];
      connect2[1] = srcElem->connect[0];
      connect2[3] = srcElem->connect[2];
      coord2[1] = d0;
      coord2[3] = d2;
      if (fabs(r) < tol_r) {
        connect2[0] = srcElem->connect[3];
        connect2[2] = srcElem->connect[1];
        coord2[0] = d3;
        coord2[2] = d1;
      }
      else {
        connect2[0] = srcElem->connect[1];
        connect2[2] = srcElem->connect[3];
        coord2[0] = d1;
        coord2[2] = d3;
      }
    }

  }
  else
  {
    t[0] = s[0];
    t[2] = s[2];

    if(s[3] < s[1])
    {
      t[1] = s[3];
      t[3] = s[1];
      connect2[1] = srcElem->connect[2];
      connect2[3] = srcElem->connect[0];
      coord2[1] = d2;
      coord2[3] = d0;
      if (fabs(r) < tol_r) {
        connect2[0] = srcElem->connect[1];
        connect2[2] = srcElem->connect[3];
        coord2[0] = d1;
        coord2[2] = d3;
      }
      else {
        connect2[0] = srcElem->connect[3];
        connect2[2] = srcElem->connect[1];
        coord2[0] = d3;
        coord2[2] = d1;
      }
    }
    else
    {
      t[1] = s[1];
      t[3] = s[3];
      connect2[0] = srcElem->connect[0];
      connect2[2] = srcElem->connect[2];
      coord2[0] = d0;
      coord2[2] = d2;
      if (fabs(r) < tol_r) {
        connect2[1] = srcElem->connect[3];
        connect2[3] = srcElem->connect[1];
        coord2[1] = d3;
        coord2[3] = d1;
      }
      else {
        connect2[1] = srcElem->connect[1];
        connect2[3] = srcElem->connect[3];
        coord2[1] = d1;
        coord2[3] = d3;
      }
    }
  }

  // Check if an intersection between rectangle1
  // and rectangle2 exists.
  if(t[0] >= 1.0)
    return false;

  if(t[2] <= 0.0)
    return false;

  if(t[1] >= 1.0)
    return false;

  if(t[3] <= 0.0)
    return false;

  IntersectionElem* genElem = new IntersectionElem();
  genElem->eNum1 = trgElem;
  genElem->eNum1 = srcElem;
  genElem->connect.Resize(4);
  genElem->nodeCoords.Resize(4);

  diffX *= distX;
  diffY *= distY;

  // If an intersection actually exist, we eventually
  // have to compute the intersection points.
  // There exist 16 different cases how two axiparallel
  // rectangles can intersect each other.

  Vector<Double> tmp;

  if(t[0] <= 0)
  {
    if(t[2] >= 1)
    {
      if(t[1] <= 0)
      {
        genElem->connect[0] = trgElem->connect[0];
        genElem->connect[1] = trgElem->connect[1];
        genElem->nodeCoords[0] = c0;
        genElem->nodeCoords[1] = c1;
        if(t[3] >= 1)
        {
          genElem->connect[2] = trgElem->connect[2];
          genElem->connect[3] = trgElem->connect[3];
          genElem->nodeCoords[2] = c2;
          genElem->nodeCoords[3] = c3;
        }
        else
        {
          tmp = c0 + diffX     + diffY*t[3];
          //do not add to grid
          //ptGrid_->AddNode(tmp, nodeNr);
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          tmp = c0 + diffX*0.0 + diffY*t[3];
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*0.0 + diffY*t[1];
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;

          tmp = c0 + diffX     + diffY*t[1];
          genElem->connect[1] = 0;
          genElem->connect[2] = trgElem->connect[2];
          genElem->connect[3] = trgElem->connect[3];
          genElem->nodeCoords[1] = tmp;
          genElem->nodeCoords[2] = c2;
          genElem->nodeCoords[3] = c3;
        }
        else
        {
          tmp = c0 + diffX*0.0 + diffY*t[1];
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;
          tmp = c0 + diffX     + diffY*t[1];
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;
          tmp = c0 + diffX     + diffY*t[3];
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          tmp = c0 + diffX*0.0 + diffY*t[3];
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }
      }
    }
    else
    {
      if(t[1] <= 0)
      {
        if(t[3] >= 1)
        {
          genElem->connect[0] = trgElem->connect[0];
          genElem->nodeCoords[0] = c0;
          tmp = c0 + diffX*t[2] + diffY*0.0;
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;
          tmp = c0 + diffX*t[2] + diffY;
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          genElem->connect[3] = trgElem->connect[3];
          genElem->nodeCoords[3] = c3;
        }
        else
        {
          genElem->connect[0] = trgElem->connect[0];
          genElem->nodeCoords[0] = c0;
          tmp = c0 + diffX*t[2] + diffY*0.0;
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;

          genElem->connect[2] = connect2[2];
          genElem->nodeCoords[2] = coord2[2];

          tmp = c0 + diffX*0.0  + diffY*t[3];
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*0.0  + diffY*t[1];
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;
          genElem->connect[1] = connect2[1];
          tmp = c0 + diffX*t[2] + diffY;
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          genElem->connect[3] = trgElem->connect[3];
          genElem->nodeCoords[3] = c3;
        }
        else
        {
          tmp = c0 + diffX*0.0  + diffY*t[1];
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;
          genElem->connect[1] = connect2[1];
          genElem->connect[2] = connect2[2];
          genElem->nodeCoords[1] = coord2[1];
          genElem->nodeCoords[2] = coord2[2];
          tmp = c0 + diffX*0.0  + diffY*t[3];
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }
      }
    }
  }
  else
  {
    if(t[2] >= 1)
    {
      if(t[1] <= 0)
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;
          genElem->connect[1] = trgElem->connect[1];
          genElem->connect[2] = trgElem->connect[2];
          genElem->nodeCoords[1] = c1;
          genElem->nodeCoords[2] = c2;
          tmp = c0 + diffX*t[0] + diffY;
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }
        else
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;

          genElem->connect[1] = trgElem->connect[1];
          genElem->nodeCoords[1] = c1;
          tmp = c0 + diffX      + diffY*t[3];
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;

          genElem->connect[3] = connect2[3];
          genElem->nodeCoords[2] = coord2[3];
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          genElem->connect[0] = connect2[0];
          genElem->nodeCoords[0] = coord2[0];
          tmp = c0 + diffX      + diffY*t[1];
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;

          genElem->connect[2] = trgElem->connect[2];
          genElem->nodeCoords[2] = c2;
          tmp = c0 + diffX*t[0] + diffY;
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }
        else
        {
          genElem->connect[0] = connect2[0];
          genElem->nodeCoords[0] = coord2[0];
          tmp = c0 + diffX      + diffY*t[1];
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;
          tmp = c0 + diffX      + diffY*t[3];
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          genElem->connect[3] = connect2[3];
          genElem->nodeCoords[3] = coord2[3];
        }
      }
    }
    else
    {
      if(t[1] <= 0)
      {
        if(t[3] >= 1)
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;
          tmp = c0 + diffX*t[2] + diffY*0.0;
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;
          tmp = c0 + diffX*t[2] + diffY;
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          tmp = c0 + diffX*t[0] + diffY;
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }
        else
        {
          tmp = c0 + diffX*t[0] + diffY*0.0;
          genElem->connect[0] = 0;
          genElem->nodeCoords[0] = tmp;
          tmp = c0 + diffX*t[2] + diffY*0.0;
          genElem->connect[1] = 0;
          genElem->nodeCoords[1] = tmp;
          genElem->connect[2] = connect2[2];
          genElem->connect[3] = connect2[3];
          genElem->nodeCoords[2] = coord2[2];
          genElem->nodeCoords[3] = coord2[3];
        }

      }
      else
      {
        if(t[3] >= 1)
        {
          genElem->connect[0] = connect2[0];
          genElem->connect[1] = connect2[1];
          genElem->nodeCoords[0] = coord2[0];
          genElem->nodeCoords[1] = coord2[1];
          tmp = c0 + diffX*t[2] + diffY;
          genElem->connect[2] = 0;
          genElem->nodeCoords[2] = tmp;
          tmp = c0 + diffX*t[0] + diffY;
          genElem->connect[3] = 0;
          genElem->nodeCoords[3] = tmp;
        }
        else
        {
          genElem->connect[0] = connect2[0];
          genElem->connect[1] = connect2[1];
          genElem->connect[2] = connect2[2];
          genElem->connect[3] = connect2[3];
          genElem->nodeCoords[0] = coord2[0];
          genElem->nodeCoords[1] = coord2[1];
          genElem->nodeCoords[2] = coord2[2];
          genElem->nodeCoords[3] = coord2[3];
        }
      }
    }
  }

  genElem->type = Elem::ET_QUAD4;
  lastIntersections_.Push_back(genElem);

  return true;
}
*/

}
