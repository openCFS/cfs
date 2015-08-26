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

EqnMapSimple::EqnMapSimple(CF::ResultInfo::EntityUnknownType type,
                           boost::shared_ptr<CF::Grid> myGrid,
                           CF::UInt eqnPerEnt, bool isOneBased)
             : mapType_(type),numEqns_(0),isFinalized_(false),
               maxNumEqns_(0),eqnPerEnt_(eqnPerEnt){

 entityEquations_ = NULL;

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
 entityEquations_ = new CF::UInt[numEntries];
 for(CF::UInt i=0;i<numEntries;i++){
   entityEquations_[i] = 0;
 }
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
      entityEquations_[aENum-zeroOne_] = (entityEquations_[aENum-zeroOne_]==0)? numEqns_ : 0;
      numEqns_+=eqnPerEnt_;
    }
  }
  isFinalized_ = true;

}

void EqnMapSimple::GetEquation(CF::StdVector<CF::UInt> & eqns,
                               const CF::UInt globalEntNum,
                               CF::ResultInfo::EntityUnknownType type) const {
  // Some mild security checks
  CF::assert(type == mapType_);
  if(globalEntNum > maxNumEqns_){
    CF::Exception("Requested equation number for an invalid entity number (Exceeds maximum)");
  }

  eqns.Resize(eqnPerEnt_);

  CF::UInt start = entityEquations_[globalEntNum-zeroOne_];
  for(CF::UInt aDof = 0;aDof < eqnPerEnt_;aDof++){
    eqns[aDof] = start+aDof;
  }

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

  for(CF::UInt aEnt = 0 ; aEnt < regEntities.GetSize() ; aEnt++){
    aEntNum = regEntities[aEnt];
    start = entityEquations_[aEntNum-zeroOne_];

    for(CF::UInt aDof = 0;aDof < eqnPerEnt_;aDof++){
      eqns[aEnt+aDof] = start+aDof;
    }
  }
}

}

