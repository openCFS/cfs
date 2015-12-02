// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AnalyticalFields.hh
 *       \brief    <Description>
 *
 *       \date     Dec 1, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_UTILS_ANALYTICALFIELDS_HH_
#define SOURCE_CFSDAT_UTILS_ANALYTICALFIELDS_HH_

#include "cfsdat/Utils/Defines.hh"
#include "MatVec/Vector.hh"

namespace CFSDat{

template<unsigned int t>
struct CylinderVortex{

  CylinderVortex(){
    pt1.Resize(t);
    pt2.Resize(t);
    radius = 0;
    rpm = 0;
  }

  CylinderVortex(CF::Vector<Double> p1, CF::Vector<Double> p2, Double rad, Double rotPerMinute){
      pt1 = p1;
      pt2 = p2;
      radius = rad;
      rpm = rotPerMinute;
  }

  inline void ComputeVortexVelocity(CF::Vector<Double>& v, const CF::Vector<Double>& point){
    EXCEPTION("Undefined dimension")
  }

  CF::Vector<Double> pt1;
  CF::Vector<Double> pt2;
  Double         radius;
  Double         rpm;
};

}


#endif /* SOURCE_CFSDAT_UTILS_ANALYTICALFIELDS_HH_ */
