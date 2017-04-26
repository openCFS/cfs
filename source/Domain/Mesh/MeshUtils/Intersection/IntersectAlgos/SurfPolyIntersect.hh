// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SurfPolyIntersect.hh
 *       \brief    <Description>
 *
 *       \date     Jan 6, 2016
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_SURFPOLYINTERSECT_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_SURFPOLYINTERSECT_HH_


#include "SurfElemIntersect.hh"

namespace CoupledField{

class SurfPolyIntersect : public SurfElemIntersect{

public:
  SurfPolyIntersect(){

  }

  virtual ~SurfPolyIntersect(){

  }

  virtual void Intersect(const SurfElem* e1, const SurfElem* e2,StdVector<NcSurfElem*> iElems) = 0;

};

}

#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_SURFPOLYINTERSECT_HH_ */
