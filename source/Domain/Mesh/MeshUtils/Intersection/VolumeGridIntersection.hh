// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     VolumeGridIntersection.hh
 *       \brief    <Description>
 *
 *       \date     Dec 6, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_VOLUMEGRIDINTERSECTION_HH_
#define SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_VOLUMEGRIDINTERSECTION_HH_

#include "Domain/ElemMapping/Elem.hh"


#include "IntersectAlgos/TriaIntersect.hh"
#include "IntersectAlgos/ElemIntersect2D.hh"
#include <unordered_map>

namespace CoupledField{

//! Class providing interface for volume grid intersection

//! Can provide the user with the information about all intersection centroids and volumes
//! or provides a real list with intersection elements for further computations
template<class INTER>
class VolumeGridIntersection{

public:

  //!constructor
  //!\param(in) target pointer to base grid
  //!\param(in) source pointer to source grid
  //!\param(in) targetRegions set of active regions in the target grid
  //!\param(in) sourceRegions set of active regions in the source grid
  //!\param(in) tolerance geometric tolerance for intersection algorithms
  VolumeGridIntersection(Grid* target, Grid* source, std::set<std::string> targetRegions, std::set<std::string> sourceRegions, Double tolerance = 1e-13);

  //! destructor
  ~VolumeGridIntersection(){

  }

  //! retuns interseciton result as a list of intersection elements
  //! usually all tetra or triangles with information about original elements
  StdVector<IntersectionElem*>  GetElemIntersection();

  //! returns list of intersection centroids and volumes
  StdVector<ElemIntersect::VolCenterInfo> GetVolCenterInfo();

private:

  //! determines possible intersection candidates by emplying
  //! an element bounding box intersection
  void PreComputeCandidates();

  //! vector of intersection candidates
  std::vector< std::pair<UInt, UInt > >  elemCandidates_;

  //! pointer to base/target grid
  Grid* tGrid_;

  //! pointer to operating/source grid
  Grid* sGrid_;

  //! set of acive target regions
  std::set<std::string> tRegs_;

  //! set of acive source regions
  std::set<std::string> sRegs_;

  //! geimetric tolerance of algorithms
  Double tol_;
};

}

#endif /* SOURCE_DOMAIN_MESH_MESHUTILS_INTERSECTION_VOLUMEGRIDINTERSECTION_HH_ */
