// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     LineIntersect.hh
 *       \brief    <Description>
 *
 *       \date     Jan 6, 2016
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_LINEINTERSECT_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_LINEINTERSECT_HH_


#include "SurfElemIntersect.hh"

namespace CoupledField{

class LineIntersect : public SurfElemIntersect{

public:
  LineIntersect(){

  }

  virtual ~LineIntersect(){

  }

  virtual void Intersect(const SurfElem* e1, const SurfElem* e2,StdVector<NcSurfElem*> iElems) = 0;

};

}

#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_LINEINTERSECT_HH_ */
