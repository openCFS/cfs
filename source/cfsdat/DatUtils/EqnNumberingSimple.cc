// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     EqnNumberingSimple.cc
 *       \brief    <Description>
 *
 *       \date     Aug 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#include "EqnNumberingSimple.hh"

namespace CFSDat{

namespace CF = CoupledField;

EqnMapSimple::EqnMapSimple(CF::ResultInfo::EntityUnknownType type,
                           CF::Grid* myGrid,
                           CF::UInt eqnPerEnt, bool isOneBased)
             : mapType_(type),eqnPerEnt_(eqnPerEnt),numEqns_(0),
               maxNumEqns_(0),zeroOne_(0),isFinalized_(false){

 ptGrid_ = myGrid;
 //obtain maximum number of entities of given type in Grid
 CF::UInt numEntries = 0;
 if(mapType_==CF::ResultInfo::ELEMENT){
   numEntries = ptGrid_->GetNumElems(CF::ALL_REGIONS);
 }else if(mapType_==CF::ResultInfo::NODE){
   numEntries = ptGrid_->GetNumNodes(CF::ALL_REGIONS);
 }else {
   CF::Exception("mapType not supported for simple equation numbering");
 }

 // zero and one based entity numbers
 zeroOne_ = (isOneBased)? 1 : 0;

 maxNumEqns_ = numEntries;
 UInt vecSize = (isOneBased) ? maxNumEqns_ : maxNumEqns_+1;
 entityEquations_.Resize(vecSize);
 entityEquations_.Init(0);
 usedEntity_.clear();
 usedEntity_.resize(vecSize, false);
}

void EqnMapSimple::Finalize(){

  assert(mapType_==CF::ResultInfo::ELEMENT || mapType_==CF::ResultInfo::NODE);

  if(isFinalized_)
    CF::Exception("EqnMapSimple::Finalize() may only be called once");

  CF::UInt numRegions = regions_.GetSize();
  boost::shared_ptr<CF::EntityList> currentList;
  CF::StdVector<CF::UInt> regEntities;
  
  for(CF::UInt aReg = 0; aReg < numRegions; aReg++){

    std::string regName = ptGrid_->GetRegion().ToString(regions_[aReg]);
    if(mapType_==CF::ResultInfo::ELEMENT){
      ptGrid_->GetElemNumsByName(regEntities, regName);
    }else{
      ptGrid_->GetNodesByName(regEntities, regName);
    }

    for(CF::UInt i=0 ; i < regEntities.GetSize() ; i++){
      CF::UInt aENum = regEntities[i];
      if(!usedEntity_[aENum-zeroOne_]){
        usedEntity_[aENum-zeroOne_] = true;
        entityEquations_[aENum-zeroOne_] = numEqns_;
        numEqns_+=eqnPerEnt_;
      }
    }
  }
  isFinalized_ = true;
}


void EqnMapSimple::GetEquation(CF::StdVector<CF::UInt> & eqns,
                               const CF::UInt globalEntNum,
                               CF::ResultInfo::EntityUnknownType type) const {
  // Some mild security checks
  assert(type == mapType_);
  if(globalEntNum > maxNumEqns_){
    CF::Exception("Requested equation number for an invalid entity number (Exceeds maximum)");
  }

  eqns.Resize(eqnPerEnt_);

  CF::UInt start = entityEquations_[globalEntNum-zeroOne_];
  for(CF::UInt aDof = 0;aDof < eqnPerEnt_;aDof++){
    eqns[aDof] = start+aDof;
  }

}

CF::UInt EqnMapSimple::GetEntityIndex(const CF::UInt globalEntNum) const {
  // Some mild security checks
  if(globalEntNum > maxNumEqns_){
    CF::Exception("EqnNumberingSimple::GetEntityIndex: Requested equation number for an invalid entity number (Exceeds maximum)");
  }
  return entityEquations_[globalEntNum-zeroOne_] / eqnPerEnt_;
}

bool EqnMapSimple::IsEntityUsed(const UInt globalEntNum) const {
  return usedEntity_[globalEntNum-zeroOne_];  
}

void EqnMapSimple::GetRegionEquations(CF::StdVector<CF::UInt> & eqns, CF::RegionIdType region) const {
  assert(mapType_==CF::ResultInfo::ELEMENT || mapType_==CF::ResultInfo::NODE);

  CF::UInt aEntNum = 0;
  CF::UInt start = 0;
  CF::StdVector<CF::UInt> regEntities;
  std::string regName = ptGrid_->GetRegion().ToString(region);

  if(mapType_==CF::ResultInfo::ELEMENT){
    ptGrid_->GetElemNumsByName(regEntities, regName);
  }else{
    ptGrid_->GetNodesByName(regEntities, regName);
  }

  eqns.Resize(regEntities.GetSize() * eqnPerEnt_);
  eqns.Init(2000000000);
  for(CF::UInt aEnt = 0 ; aEnt < regEntities.GetSize() ; aEnt++){
    aEntNum = regEntities[aEnt];
    start = entityEquations_[aEntNum-zeroOne_];

    for(CF::UInt aDof = 0;aDof < eqnPerEnt_;aDof++){
      eqns[(aEnt*eqnPerEnt_)+aDof] = start+aDof;
    }
  }
}

void EqnMapSimple::GetSubsetEquations(CF::StdVector<UInt> & eqns, CF::StdVector<UInt> & globalEntNumbers) const{
  UInt aEntNum = 0;
  UInt start = 0;

  eqns.Resize(globalEntNumbers.GetSize() * eqnPerEnt_);
  for(CF::UInt aEnt = 0 ; aEnt < globalEntNumbers.GetSize() ; aEnt++){
    aEntNum = globalEntNumbers[aEnt];
    start = entityEquations_[aEntNum-zeroOne_];

    for(CF::UInt aDof = 0;aDof < eqnPerEnt_;aDof++){
      eqns[aEnt+aDof] = start+aDof;
    }
  }
}

void EqnMapSimple::GetReverseEntityMap(CF::StdVector<UInt>& revMap) const {
  revMap.Resize(GetNumEntities());
  for (UInt i = 1; i <= maxNumEqns_; i++) {
    if (IsEntityUsed(i)) {
      revMap[GetEntityIndex(i)] = i;
    }
  }
}

}

