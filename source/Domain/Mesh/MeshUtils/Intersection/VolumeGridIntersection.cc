// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     VolumeGridIntersection.cc
 *       \brief    <Description>
 *
 *       \date     Dec 15, 2015
 *       \author   ahueppe
 */
//================================================================================================

#include "VolumeGridIntersection.hh"
#include "Domain/Mesh/MeshUtils/EntityAssociation.hh"
#include "Domain/Mesh/Grid.hh"

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
const UInt numThreads=omp_get_num_threads();
#endif

#ifndef USE_OPENMP
const UInt numThreads=1;
#endif

namespace CoupledField{

template<class INTER>
VolumeGridIntersection<INTER>::VolumeGridIntersection(Grid* target, Grid* source,
                                               std::set<std::string> targetRegions, std::set<std::string> sourceRegions,
                                               Double tolerance)
                       :tGrid_(target),sGrid_(source),
                        tRegs_(targetRegions),sRegs_(sourceRegions), tol_(tolerance),
                        interAlgo_(INTER(tGrid_,sGrid_)){

  if(tGrid_->GetDim() != sGrid_->GetDim()){
    //we would need additional implementation e.g. to map a 3D surface region to a 2D Grid.
    EXCEPTION("Grids need to have same dimensions at the moment. Projections 3D->2D are not supported right now.");
  }

  //check for dimensional compatibility
  std::set<std::string>::iterator tIter = tRegs_.begin();
  std::set<std::string>::iterator sIter = sRegs_.begin();
  UInt gDim = target->GetDim();
  for(;tIter != tRegs_.end(); ++tIter){
    UInt tDim = tGrid_->GetEntityDim(*tIter);
    if(gDim != tDim){
      EXCEPTION("Incompatible region dimensions in Volume intersection");
    }
  }
  for(;sIter != sRegs_.end(); ++sIter){
    UInt sDim = sGrid_->GetEntityDim(*sIter);
    if(gDim != sDim){
      EXCEPTION("Incompatible region dimensions in Volume intersection");
    }
  }
}

template<class INTER>
StdVector<ElemIntersect::VolCenterInfo> VolumeGridIntersection<INTER>::GetVolCenterInfo(){

  PreComputeCandidates();

  //now loop over each candidate and call for the intersection algorithm
  UInt numCandidates = elemCandidates_.size();
  StdVector<ElemIntersect::VolCenterInfo> retInfo;
  UInt numIntersects;
#pragma omp parallel shared(numIntersects)
{

  INTER locAlgo = interAlgo_;
  StdVector<Double> curVols;
  StdVector<Vector<Double> > curCenters;
  StdVector<ElemIntersect::VolCenterInfo> localInfo;
  StdVector<ElemIntersect::VolCenterInfo> elemInfo;
  localInfo.Reserve(std::ceil(numCandidates/numThreads));
#pragma omp for
  for(UInt aCand = 0; aCand < numCandidates; ++aCand){
    //obtain element pointers from grids
    const std::pair<UInt,UInt>& candidate = elemCandidates_[aCand];
    const Elem* elem1 = tGrid_->GetElem(candidate.first);
    const Elem* elem2 = sGrid_->GetElem(candidate.second);
    locAlgo.SetTElem(elem1->elemNum);
    bool isIntersect = locAlgo.Intersect(elem2->elemNum);
    if(isIntersect){
      elemInfo.Clear(true);
      locAlgo.GetVolumeAndCenters(elemInfo);
      for(UInt i=0;i<curVols.GetSize(); ++i){
        numIntersects++;
        localInfo.Push_back(elemInfo[i]);
      }
    }
  }
#pragma omp flush(numIntersects)
#pragma omp single
  {
  retInfo.Reserve(numIntersects);
  }
#pragma omp critical
  {
    for(UInt aRet =0;localInfo.GetSize(); ++aRet){
      retInfo.Push_back(localInfo[aRet]);
    }
  }
}
 return retInfo;

}

template<class INTER>
void VolumeGridIntersection<INTER>::PreComputeCandidates(){
  //try to prepare the intersection
  // i.e. we compute intersection candidates according to bounding boxes
  StdVector<Elem*> targetElems;
  StdVector<Elem*> sourceElems;
  StdVector<Elem*> cElems;

  UInt numSElems = 0;
  UInt numTElems = 0;

  RegionIdType cId;
  std::set<std::string>::iterator tIter = tRegs_.begin();
  std::set<std::string>::iterator sIter = sRegs_.begin();

  for(;tIter != tRegs_.end(); ++tIter){
    cId = tGrid_->GetRegion().Parse(*tIter);
    numTElems += tGrid_->GetNumElems(cId);
  }
  for(;sIter != sRegs_.end(); ++sIter){
    cId = sGrid_->GetRegion().Parse(*sIter);
    numSElems += sGrid_->GetNumElems(cId);
  }

  targetElems.Reserve(numTElems);
  sourceElems.Reserve(numSElems);

  tIter = tRegs_.begin();
  sIter = sRegs_.begin();
  for(;tIter != tRegs_.end(); ++tIter){
    cId = tGrid_->GetRegion().Parse(*tIter);
    tGrid_->GetVolElems(cElems,cId);
    for(UInt aE = 0;aE<cElems.GetSize();++aE){
      targetElems.Push_back(cElems[aE]);
    }
  }
  for(;sIter != sRegs_.end(); ++sIter){
    cId = sGrid_->GetRegion().Parse(*sIter);
    sGrid_->GetVolElems(cElems,cId);
    for(UInt aE = 0;aE<cElems.GetSize();++aE){
      sourceElems.Push_back(cElems[aE]);
    }
  }

  elemCandidates_ = BoundingBoxAssociate::AssociateEntities(targetElems,sourceElems,tGrid_,sGrid_,tol_);

}

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class VolumeGridIntersection<TriaIntersect>;
#endif

}
