
// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CutTetra.cc
 *       \brief    <Description>
 *
 *       \date     Dec 11, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include <Domain/Mesh/MeshUtils/Intersection/IntersectAlgos/CoordTetra.hh>

namespace CoupledField{

void CoordTetra::ComputeClipPlanes(){
  // Define edge vectors
  typedef Vector<Double> vec;

  vec e10 = points[1] - points[0];
  vec e20 = points[2] - points[0];
  vec e30 = points[3] - points[0];
  vec e21 = points[2] - points[1];
  vec e31 = points[3] - points[1];

  // Cross-products
  e20.CrossProduct(e10,clipPlanes_[0].first);
  e10.CrossProduct(e30,clipPlanes_[1].first);
  e30.CrossProduct(e20,clipPlanes_[2].first);
  e21.CrossProduct(e31,clipPlanes_[3].first);

  // Normalize
  clipPlanes_[0].first /= clipPlanes_[0].first.NormL2() + EPS;
  clipPlanes_[1].first /= clipPlanes_[1].first.NormL2() + EPS;
  clipPlanes_[2].first /= clipPlanes_[2].first.NormL2() + EPS;
  clipPlanes_[3].first /= clipPlanes_[3].first.NormL2() + EPS;

  // Compute magnitude of clipping tetrahedron
  tetSize_ = (1.0 / 6.0) * (e10 * clipPlanes_[3].first);
  if (tetSize_ < 0.0)
  {
    // Reverse normal directions
    clipPlanes_[0].first = -clipPlanes_[0].first;
    clipPlanes_[1].first = -clipPlanes_[1].first;
    clipPlanes_[2].first = -clipPlanes_[2].first;
    clipPlanes_[3].first = -clipPlanes_[3].first;
    // Reverse sign
    tetSize_ = std::abs(tetSize_);
  }
  // Determine plane constants
  clipPlanes_[0].second = (points[0] * clipPlanes_[0].first);
  clipPlanes_[1].second = (points[1] * clipPlanes_[1].first);
  clipPlanes_[2].second = (points[2] * clipPlanes_[2].first);
  clipPlanes_[3].second = (points[3] * clipPlanes_[3].first);
}


}



