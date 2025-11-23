// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AnalyticalFields.cc
 *       \brief    <Description>
 *
 *       \date     Dec 1, 2015
 *       \author   ahueppe
 */
//================================================================================================


//#include "AnalyticalFields.hh"

namespace CFSDat{



template<>
void CylinderVortex<3>::ComputeVortexVelocity(CF::Vector<Double>& v, const CF::Vector<Double>& point ){
  dist.Init(0.0);
  pd.Init(0.0);
  tmp.Init(0.0);
  Double dot,dsq,pdLen,scalProj;
  Double radius_sq = radius*radius;
  Double lenSq;

  v.Resize(3);
  v.Init();

  dist = this->pt2 - this->pt1; // translate so pt1 is origin.  Make vector from
  lenSq = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];

  pd = point - this->pt1;

  dot = pd * dist;
  if( dot > 0.0f && dot < lenSq ){
    pdLen = (pd[0]*pd[0] + pd[1]*pd[1] + pd[2]*pd[2]);
    scalProj = dot*dot/lenSq;
    dsq = pdLen - scalProj;
    if( dsq <= radius_sq ){
      //find way along axis
      dist *= scalProj/dot;
      tmp = pd - dist;
      dist.CrossProduct(tmp,v);
      //may be unnecessary
      v.Normalize();
      v *=  ((this->rpm/60) * 2 * M_PI * sqrt(dsq));
    }
  }
}

template<>
void CylinderVortex<2>::ComputeVortexVelocity(CF::Vector<Double>& v, const CF::Vector<Double>& point){

  //obtain center
  CF::Vector<Double> dist(2);
  dist = v - pt1;
  v.Resize(2);
  v.Init();
  Double rad = dist.NormL2();
  if(rad <= radius){
    v[0] = -(rpm/60) * 2 * M_PI *  (v[1] - pt1[1]);
    v[1] =  (rpm/60) * 2 * M_PI *  (v[0] - pt1[0]);
  }

}


//    //found somewhere but seemd ok
//    Vector<Double> dist;
//    Vector<Double> pd;
//    Double dot,dsq;
//    Double radius_sq = radius*radius;
//    Double lenSq;
//
//
//    dist = pt1 - pt2; // translate so pt1 is origin.  Make vector from
//    lenSq = dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2];
//
//    pd = v - pt1;
//
//    dot = pd * dist;
//    if( dot < 0.0f || dot > lenSq ){
//       return false;
//    }else{
//      dsq = (pd[0]*pd[0] + pd[1]*pd[1] + pd[2]*pd[2]) - (dot*dot/lenSq);
//
//      if( dsq > radius_sq )
//        return false;
//      else
//        return true;
//    }
//
//
//  pd = v[0] - p1x;   // vector from pt1 to test point.
//  pdy = v[1] - p1y;
//  pdz = v[2] - p1z;
//
//  // Dot the d and pd vectors to see if point lies behind the
//  // cylinder cap at pt1.x, pt1.y, pt1.z
//
//  dot = pdx * dx + pdy * dy + pdz * dz;
//
//  // If dot is less than zero the point is behind the pt1 cap.
//  // If greater than the cylinder axis line segment length squared
//  // then the point is outside the other end cap at pt2.
//  if( dot < 0.0f || dot > lenSq ){
//    return false;
//  }else{
//    // Point lies within the parallel caps, so find
//    // distance squared from point to line, using the fact that sin^2 + cos^2 = 1
//    // the dot = cos() * |d||pd|, and cross*cross = sin^2 * |d|^2 * |pd|^2
//    // Carefull: '*' means mult for scalars and dotproduct for vectors
//    // In short, where dist is pt distance to cyl axis:
//    // dist = sin( pd to d ) * |pd|
//    // distsq = dsq = (1 - cos^2( pd to d)) * |pd|^2
//    // dsq = ( 1 - (pd * d)^2 / (|pd|^2 * |d|^2) ) * |pd|^2
//    // dsq = pd * pd - dot * dot / lengthsq
//    //  where lengthsq is d*d or |d|^2 that is passed into this function
//
//    // distance squared to the cylinder axis:
//    dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot/lenSq;
//
//    if( dsq > radius_sq )
//      return false;
//    else
//      return true;
//  }


template struct CylinderVortex<2>;
template struct CylinderVortex<3>;
}
