// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SurfaceGridIntersection.hh
 *       \brief    <Description>
 *
 *       \date     Dec 6, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_SURFACEGRIDINTERSECTION_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_SURFACEGRIDINTERSECTION_HH_

#include "Domain/ElemMapping/SurfElem.hh"

namespace CoupledField{

class SurfaceGridIntersection{

public:
  SurfaceGridIntersection(Grid* target, Grid* source, StdVector<RegionIdType> targetRegions, StdVector<RegionIdType> sourceRegions);

  ~SurfaceGridIntersection();

  StdVector<UInt> GetIntersectingElements(UInt tENum);

  void GetElemIntersection(UInt tEnum, StdVector<NcSurfElem*> iElements);

};

}
#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_SURFACEGRIDINTERSECTION_HH_ */
