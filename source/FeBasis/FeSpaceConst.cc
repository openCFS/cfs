// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//=================================
/*
 * \file   FeSpaceConst.cc
 * \brief  see FeSpaceConst.hh
 *
 * \date   May 6, 2014
 * \author dperchto
 */
//=================================

#include "FeBasis/FeSpaceConst.hh"

namespace CoupledField {

FeSpaceConst::FeSpaceConst(PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid)
  : FeSpace::FeSpace(paramNode, infoNode, ptGrid){

  type_ = CONSTANT;

  allowedEntities_.insert(EntityList::COIL_LIST);

}

FeSpaceConst::~FeSpaceConst(){}

void FeSpaceConst::Init( shared_ptr<SolStrategy> solStrat ){

  solStrat_ = solStrat;

}

BaseFE* FeSpaceConst::GetFe( const EntityIterator ent ){

  // used by GetNumFunctions in FeSpace
  return(NULL);

}

BaseFE* FeSpaceConst::GetFe( const EntityIterator ent,
                             IntScheme::IntegMethod& method,
                             IntegOrder & order ){

  return(NULL);

}

BaseFE* FeSpaceConst::GetFe( UInt elemNum ){

  return(NULL);

}

UInt FeSpaceConst::GetNumFunctions( const EntityIterator ent ){

  // uses GetFe in FeSpace
  return(1);

}

void FeSpaceConst::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){

  this->CheckEntityType(ent);

  // determine equation number by trying to insert the entity to the equation map
  std::pair<std::map<std::string,Integer>::iterator,bool> ret;
  ret = equationMap_.insert( std::pair<std::string,Integer>( ent.GetIdString(),
      (Integer)numEqns_ + 1 ) );
  eqns.Resize(1);
  if( ret.second ){
    // new entry in equation map was created
    numEqns_++;
    numFreeEquations_++;
    eqns[0] = numEqns_;
  } else {
    // entry already existed, equation number is known
    eqns[0] = ret.first->second;
  }

}

void FeSpaceConst::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                            UInt dof ){

  this->GetEqns(eqns, ent);

}

void FeSpaceConst::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                            UInt dof, BaseFE::EntityType ){

  this->GetEqns(eqns, ent);

}

void FeSpaceConst::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                            BaseFE::EntityType ){

  this->GetEqns(eqns, ent);

}

void FeSpaceConst::GetElemEqns(StdVector<Integer>& eqns, const Elem* elem){

  EXCEPTION("This space does not have elements.");

}

void FeSpaceConst::GetElemEqns(StdVector<Integer>& eqns, const Elem* elem, UInt dof){

  EXCEPTION("This space does not have elements.");

}

/*void FeSpaceConst::GetEntityListEqns( StdVector<Integer>& eqns,
                                      shared_ptr<EntityList> ent,
                                      BaseFE::EntityType entType ){

  std::set<Integer> allEqns;
  StdVector<Integer> tmp;
  EntityIterator it = ent->GetIterator();
  for( it.Begin(); ! it.IsEnd(); it++ ) {
    tmp.Init();
    this->GetEqns( tmp, it, entType);
    allEqns.insert(tmp.Begin(), tmp.End());
  }
  eqns.Clear();
  eqns.Resize(allEqns.size());
  eqns.Init();
  std::copy(allEqns.begin(), allEqns.end(), eqns.Begin());

}

void FeSpaceConst::GetEntityListEqns( StdVector<Integer>& eqns,
                                      shared_ptr<EntityList> ent, UInt dof,
                                      BaseFE::EntityType entType ){

  std::set<Integer> allEqns;
  StdVector<Integer> tmp;
  EntityIterator it = ent->GetIterator();
  for( it.Begin(); ! it.IsEnd(); it++ ) {
    tmp.Init();
    GetEqns( tmp, it, dof, entType );
    allEqns.insert(tmp.Begin(), tmp.End());
  }
  eqns.Clear();
  eqns.Resize(allEqns.size());
  std::copy(allEqns.begin(), allEqns.end(), eqns.Begin());

}*/

void FeSpaceConst::Finalize(){

  isFinalized_ = true;

}

void FeSpaceConst::MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support,
                                     shared_ptr<CoefFunction> coefFct,
                                     shared_ptr<BaseFeFunction> feFct,
                                     std::map<Integer, Double>& vals,
                                     bool cache,
                                     const std::set<UInt>& comp ){

  EXCEPTION("Neee.");

}

void FeSpaceConst::MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support,
                                     shared_ptr<CoefFunction> coefFct,
                                     shared_ptr<BaseFeFunction> feFct,
                                     std::map<Integer, Complex>& vals,
                                     bool cache,
                                     const std::set<UInt>& comp ){

  EXCEPTION("Neee.");

}

bool FeSpaceConst::IsSameEntityApproximation( shared_ptr<EntityList> list,
                                              shared_ptr<FeSpace> space ){

  return( space->GetSpaceType() == FeSpace::CONSTANT );

}

void FeSpaceConst::SetRegionElements( RegionIdType region, MappingType mType,
                                      const ApproxOrder& order,
                                      PtrParamNode infoNode ){

  EXCEPTION("This space does not approximate spacially.");

}

void FeSpaceConst::CheckConsistency(){

  // no check so far

}

void FeSpaceConst::SetDefaultIntegration(PtrParamNode infoNode){

  // might be important to guarantee a default case so everything else works

}

void FeSpaceConst::SetDefaultElements(PtrParamNode infoNode){

  // might be important to guarantee a default case so everything else works

}

void FeSpaceConst::MapNodalBCs(){

  // ignore BCs totally

}

void FeSpaceConst::CheckEntityType(const EntityIterator ent) const {

  if( allowedEntities_.find(ent.GetType()) == allowedEntities_.end() ){
    EXCEPTION("Entity type not allowed.");
  }

}

}
