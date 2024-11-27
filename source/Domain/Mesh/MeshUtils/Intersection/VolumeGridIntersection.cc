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
#include "General/Environment.hh"
#include <fstream>

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace CoupledField{

template<class INTER>
VolumeGridIntersection<INTER>::VolumeGridIntersection(Grid* target, Grid* source,
                                               std::set<std::string> targetRegions, std::set<std::string> sourceRegions,
                                               Double tolerance)
                       :tGrid_(target),sGrid_(source),
                        tRegs_(targetRegions),sRegs_(sourceRegions), tol_(tolerance){

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
  //std::vector< std::pair<UInt,UInt> > newVec = elemCandidates_;
  StdVector<ElemIntersect::VolCenterInfo> retInfo;
  UInt numIntersects = 0;
  UInt intersectsOutsideElem = 0;
  UInt intersectsSmall = 0;
  UInt numThreads= CFS_NUM_THREADS;

  std::cout << "\t\t\t Processing " << numCandidates << " intersection candidates using " <<  numThreads << " threads." << std::endl;


#pragma omp parallel num_threads(numThreads)
{
  INTER locAlgo(tGrid_,sGrid_);
  StdVector<Double> curVols;
  Vector<Double> locPoint(tGrid_->GetDim());
  StdVector<Vector<Double> > curCenters;
  StdVector<ElemIntersect::VolCenterInfo> localInfo;
  StdVector<ElemIntersect::VolCenterInfo> elemInfo;

  localInfo.Reserve(2*std::ceil(numCandidates/numThreads));
  UInt tENum = 0;
  UInt sENum = 0;

  const Elem* elem1;
  const Elem* elem2;
#pragma omp for reduction(+ : intersectsOutsideElem , numIntersects , intersectsSmall)  schedule(guided,10)
  for(Integer aCand = 0; aCand < (Integer) numCandidates; ++aCand){
    //obtain element pointers from grids
    tENum = elemCandidates_[aCand].first;
    sENum = elemCandidates_[aCand].second;
    elem1 = tGrid_->GetElem(tENum);
    elem2 = sGrid_->GetElem(sENum);
    locAlgo.SetTElem(elem1->elemNum);
    bool isIntersect = locAlgo.Intersect(elem2->elemNum);
    if(isIntersect){
      elemInfo.Clear(true);
      locAlgo.GetVolumeAndCenters(elemInfo);
      for(UInt i=0;i<elemInfo.GetSize(); ++i){
        //reject very small elements
        if(elemInfo[i].volume < 1e-30){
          intersectsSmall++;
          continue;
        }
        // shared_ptr<ElemShapeMap> esm = tGrid_->GetElemShapeMap(elem1,true);
        // tGrid_->UpdateIndividualElemShapeMap(elem1, true); // elem1 = const
        // shared_ptr<ElemShapeMap> esm(elem1->ptrShapeMap);
        shared_ptr<ElemShapeMap> esm = (elem1)->GetElemShapeMap(tGrid_, true);

        esm->Global2Local(locPoint,elemInfo[i].center);
        if(!esm->CoordIsInsideElem(locPoint,1e-6)){
          intersectsOutsideElem++;
          //locAlgo.DumpLastIntersect();
          //std::cout << locPoint << std::endl;
          //exit(0);
          continue;
        }
        numIntersects++;
        elemInfo[i].targetElemNum = elem1->elemNum;
        elemInfo[i].sourceElemNum = elem2->elemNum;
        localInfo.Push_back(elemInfo[i]);
      }
    }
  }
#pragma omp flush(numIntersects)
#pragma omp single
  {
  retInfo.Reserve(numIntersects);
  }
#pragma omp critical (VolumeGridIntersection)
  {
    for(UInt aRet =0;aRet<localInfo.GetSize(); ++aRet){
      retInfo.Push_back(localInfo[aRet]);
    }
  }

}

  //clear the element candidates as it is an expensive structure
  elemCandidates_.clear();

  //give overall volume of intersection
  Double volume = 0;
  for(UInt aRet =0;aRet<retInfo.GetSize(); ++aRet){
   volume += retInfo[aRet].volume;
 }
  //std::ofstream aStream("centers.csv", std::ios::trunc | std::ios::out);
  //for(UInt aRet =0;aRet<retInfo.GetSize(); ++aRet){
  //  aStream << retInfo[aRet].center[0] << "," << retInfo[aRet].center[1] << "," << retInfo[aRet].center[2] << std::endl;
  //}
  //aStream.close();
  std::cout << "\t\t\t Rejected " << intersectsOutsideElem << " intesection elements as the computed center is outside the taget element " << std::endl;
  std::cout << "\t\t\t Rejected " << intersectsSmall << " intesection elements due to very small volume " << std::endl;
  std::cout << "\t\t\t Computed intersection Volume is " << volume << " from " << numIntersects << " intersection elements" << std::endl;

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
template class VolumeGridIntersection<TetraIntersect>;
template class VolumeGridIntersection<ElemIntersect2D>;

}
