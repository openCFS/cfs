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
#include "PDE/MagEdgePDE.hh"
#include "PDE/MagneticPDE.hh"

namespace CoupledField {

FeSpaceConst::FeSpaceConst(PtrParamNode paramNode, PtrParamNode infoNode, Grid* ptGrid)
  : FeSpace(paramNode, infoNode, ptGrid){

  type_ = CONSTANT;

  allowedEntities_.insert(EntityList::COIL_LIST);

}

FeSpaceConst::~FeSpaceConst(){}

void FeSpaceConst::Init( shared_ptr<SolStrategy> solStrat ){

  solStrat_ = solStrat;

}

BaseFE* FeSpaceConst::GetFe( const EntityIterator ent ){

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

  return(1);

}

void FeSpaceConst::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){

  eqns.Clear();
  eqns.Resize( 1, equationMap_.at(ent.GetIdString()) );

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

  this->GetElemEqns(eqns, elem);

}

void FeSpaceConst::Finalize(){

  // take all entities from the FeFunction and generate equation map based on entity id string
  // every unique id gets an own equation

  shared_ptr<BaseFeFunction> feFct = feFunction_.lock();
  StdVector<shared_ptr<EntityList> > entListVec = feFct->GetEntityList();

  for( UInt k = 0; k < entListVec.GetSize(); k++ ){
    EntityIterator entIt = entListVec[k]->GetIterator();
    this->CheckEntityType(entIt);
    while ( !(entIt.IsEnd()) ){
      std::pair<boost::unordered_map<std::string,Integer>::iterator,bool> ret;
      ret = equationMap_.insert( std::pair<std::string,Integer>(
          entIt.GetIdString(), (Integer)numEqns_ + 1) );
      if( ret.second ){
        // new entry in equation map was created
        numEqns_++;
        numFreeEquations_++;
      }
      entIt++;
    }
  }

  isFinalized_ = true;

}

void FeSpaceConst::MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support,
                                     shared_ptr<CoefFunction> coefFct,
                                     shared_ptr<BaseFeFunction> feFct,
                                     std::map<Integer, Double>& vals,
                                     bool cache,
                                     const std::set<UInt>& comp ){

  EXCEPTION("This FeSpace does not approximate space.");

}

void FeSpaceConst::MapCoefFctToSpace(StdVector<shared_ptr<EntityList> > support,
                                     shared_ptr<CoefFunction> coefFct,
                                     shared_ptr<BaseFeFunction> feFct,
                                     std::map<Integer, Complex>& vals,
                                     bool cache,
                                     const std::set<UInt>& comp ){

  EXCEPTION("This FeSpace does not approximate space.");

}

bool FeSpaceConst::IsSameEntityApproximation( shared_ptr<EntityList> list,
                                              shared_ptr<FeSpace> space ){

  return( space->GetSpaceType() == FeSpace::CONSTANT );

}

void FeSpaceConst::SetRegionElements( RegionIdType region, MappingType mType,
                                      const ApproxOrder& order,
                                      PtrParamNode infoNode ){

  EXCEPTION("This FeSpace does not approximate space.");

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
