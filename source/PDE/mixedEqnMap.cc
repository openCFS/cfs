#include "mixedEqnMap.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh" 
#include "Domain/ansatzFct.hh"
#include "Utils/coordSystem.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "Domain/entityList.hh"

#include <boost/lexical_cast.hpp> 

using std::string; 
using boost::lexical_cast; 

namespace CoupledField {
  // declare class specific logging stream
  DECLARE_LOG(mixedEqnMap)
  DEFINE_LOG(mixedEqnMap, "mixedEqnMap")

  MixedEqnMap::MixedEqnMap ( Grid* ptGrid, PdeIdType pdeId, bool usePenalty):EqnMap(ptGrid,pdeId,usePenalty){
    this->contMap = new EqnMap( ptGrid, pdeId, usePenalty);
    this->disContMap = new DiscontinuousEqnMap( ptGrid, pdeId, usePenalty);
  }

 MixedEqnMap::~MixedEqnMap () {
  }

  void MixedEqnMap::AddResult(ResultInfo& result, shared_ptr<EntityList> list){
    if(result.fctType->IsDiscontinuous()){
      disContMap->AddResult(result,list);
    }else{
      contMap->AddResult(result,list);
    }
  }

  void MixedEqnMap::SetHomoDirichletBCs( HdBcList& hdBcs ){
    for( UInt i = 0; i< hdBcs.GetSize(); i++ ) {
      ResultInfo actResult = *(hdBcs[i]->result);
      if(actResult.fctType->IsDiscontinuous()){
        disContMap->SetHomoDirichletBCs(hdBcs);
      }else{
        contMap->SetHomoDirichletBCs(hdBcs);
      }
    }
  }
      
  void MixedEqnMap::SetInhomDirichletBCs( IdBcList& idBcs ){
    for( UInt i = 0; i< idBcs.GetSize(); i++ ) {
      ResultInfo actResult = *(idBcs[i]->result);
      if(actResult.fctType->IsDiscontinuous()){
        disContMap->SetInhomDirichletBCs(idBcs);
      }else{
        contMap->SetInhomDirichletBCs(idBcs);
      }
    }
  }

  void MixedEqnMap::SetConstraints( ConstraintList& constraints ){
    for( UInt i = 0; i< constraints.GetSize(); i++ ) {
      ResultInfo actResult = *(constraints[i]->result);
      if(actResult.fctType->IsDiscontinuous()){
        disContMap->SetConstraints(constraints);
      }else{
        contMap->SetConstraints(constraints);
      }
    }
  }
  
  void MixedEqnMap::Finalize(){
    contMap->Finalize();
    disContMap->CopyMapInfo(contMap);
    disContMap->Finalize();
  }

  void MixedEqnMap::GetEqns( StdVector<Integer>& eqns, const ResultInfo& result, const EntityIterator& it )const {
    if(result.fctType->IsDiscontinuous()){
      disContMap->GetEqns(eqns,result,it);
    }else{
      contMap->GetEqns(eqns,result,it);
    }
  }

  void MixedEqnMap::GetEqns( StdVector<Integer>& eqns, const ResultInfo& result, const EntityIterator& it, UInt dof )const {
    if(result.fctType->IsDiscontinuous()){
      disContMap->GetEqns(eqns,result,it,dof);
    }else{
      contMap->GetEqns(eqns,result,it,dof);
    }
  }

  Integer MixedEqnMap::GetEqn( const ResultInfo& result, const EntityIterator& it, UInt dof ) const {
    if(result.fctType->IsDiscontinuous()){
     return disContMap->GetEqn(result,it,dof);
    }else{
     return contMap->GetEqn(result,it,dof);
    }
  }

  //! Name mapping function for obtaining a nodal equation
  Integer MixedEqnMap::GetNodeEqn( const ResultInfo& result, UInt nodeNr, UInt dof ){
    if(result.fctType->IsDiscontinuous()){
     return disContMap->GetNodeEqn(result,nodeNr,dof);
    }else{
     return contMap->GetNodeEqn(result,nodeNr,dof);
    }
  }

  void MixedEqnMap::ToInfo(InfoNode* in) const{
    contMap->ToInfo(in);
    disContMap->ToInfo(in);
  }
}
