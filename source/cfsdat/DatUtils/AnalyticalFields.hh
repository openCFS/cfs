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

#include "cfsdat/DatUtils/Defines.hh"
#include "MatVec/Vector.hh"

namespace CFSDat{

template<unsigned int t>
struct CylinderVortex{

  CylinderVortex() : 
    pt1(t),pt2(t),dist(t),pd(t),tmp(t)
  {
    this->pt1.Init(0);
    this->pt2.Init(0);
    this->radius = 0;
    this->rpm = 0;
  }

  CylinderVortex(CF::Vector<Double> p1, CF::Vector<Double> p2, Double rad, Double rotPerMinute) : 
  pt1(t),pt2(t),dist(t),pd(t),tmp(t) 
  {
      this->pt1 = p1;
      this->pt2 = p2;
      this->radius = rad;
      this->rpm = rotPerMinute;
  }

  void ComputeVortexVelocity(CF::Vector<Double>& v, const CF::Vector<Double>& point){
    EXCEPTION("Undefined dimension : " << t)
  }

  CF::Vector<Double> pt1;
  CF::Vector<Double> pt2;
  CF::Vector<Double> dist,pd,tmp;
  Double         radius;
  Double         rpm;
};

}

#include "AnalyticalFields.cc"
//template<>
//void CylinderVortex<2>::ComputeVortexVelocity(CF::Vector<Double>& v, const CF::Vector<Double>& point);
//
//template<>
//void CylinderVortex<3>::ComputeVortexVelocity(CF::Vector<Double>& v, const CF::Vector<Double>& point);

#endif /* SOURCE_CFSDAT_UTILS_ANALYTICALFIELDS_HH_ */
