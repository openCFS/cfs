// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FeSpaceL2.cc
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================


#include "FeSpaceL2.hh"

namespace CoupledField{

  FeSpaceL2::FeSpaceL2(PtrParamNode aNode, PtrParamNode infoNode)
   : FeSpaceH1(aNode, infoNode ) {
    type_ = L2;
    isContinuous_ = false;
  }

  FeSpaceL2::~FeSpaceL2(){
  }

  void FeSpaceL2::AddFeFunction( shared_ptr<BaseFeFunction> fct ){
    feFunction_ = fct;
    return;
  }

  UInt FeSpaceL2::GetNumFunctions( const EntityIterator ent ){
    BaseFE * fe = GetFe(ent);
    return fe->GetNumFncs();
  }

  void FeSpaceL2::GetEqns( StdVector<Integer>& eqns,
                           const EntityIterator ent ){
    ///Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      if(gridToVirtualNodes_.find(node) == gridToVirtualNodes_.end())
        return;

      UInt nDisNodes = gridToVirtualNodes_[node].GetSize();
      eqns.Resize(dofsPerUnknown*nDisNodes);
      eqns.Init();
      for (UInt iNode = 0; iNode < nDisNodes; iNode++ ) {
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[(iNode*dofsPerUnknown) + iDof] =
              nodeMap_[gridToVirtualNodes_[node][iNode]][iDof];
        }
      }
    }else if( ent.GetType() == EntityList::ELEM_LIST ||
             ent.GetType() == EntityList::SURF_ELEM_LIST){
      StdVector<UInt> nodes;
      GetNodesOfElement(nodes,ent.GetElem());
      eqns.Resize( nodes.GetSize() * dofsPerUnknown );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[(iNode*dofsPerUnknown) + iDof] =
              nodeMap_[nodes[iNode]][iDof];
        }
      }
    }else{
      EXCEPTION("In FeSpaceL2::GetEqns():  Supplied an iterator which is not supported by FeSpace");
    }
  }

  void FeSpaceL2::GetEqns( StdVector<Integer>& eqns,
                           const EntityIterator ent,
                           UInt dof ){
    ///Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      if(gridToVirtualNodes_.find(node) == gridToVirtualNodes_.end())
        return;

      UInt nDisNodes = gridToVirtualNodes_[node].GetSize();
      eqns.Resize(nDisNodes);
      eqns.Init();
      for (UInt iNode = 0; iNode < nDisNodes; iNode++ ) {
        eqns[iNode] = nodeMap_[gridToVirtualNodes_[node][iNode]][dof];
      }
    }else if( ent.GetType() == EntityList::ELEM_LIST ||
             ent.GetType() == EntityList::SURF_ELEM_LIST){
      StdVector<UInt> nodes;
      GetNodesOfElement(nodes,ent.GetElem());
      eqns.Resize( nodes.GetSize() );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqns[iNode] = nodeMap_[nodes[iNode]][dof];
      }
    }else{
      EXCEPTION("In FeSpaceL2::GetEqns():  Supplied an iterator which is not supported by FeSpace");
    }
  }

  void FeSpaceL2::GetEqns( StdVector<Integer>& eqns,
                           const EntityIterator ent,
                           UInt dof,
                           BaseFE::EntityType entityType){
    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      if(gridToVirtualNodes_.find(node) == gridToVirtualNodes_.end())
        return;

      UInt nDisNodes = gridToVirtualNodes_[node].GetSize();
      eqns.Resize(nDisNodes);
      eqns.Init();
      for (UInt iNode = 0; iNode < nDisNodes; iNode++ ) {
        eqns[iNode] = nodeMap_[gridToVirtualNodes_[node][iNode]][dof];
      }
    }else if( ent.GetType() == EntityList::ELEM_LIST ||
              ent.GetType() == EntityList::SURF_ELEM_LIST){
      StdVector<UInt> nodes;
      GetNodesOfElement(nodes,ent.GetElem(), entityType);
      eqns.Resize( nodes.GetSize() );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
       eqns[iNode] = nodeMap_[nodes[iNode]][dof];
      }
    } else {
      EXCEPTION("In FeSpaceH1::GetEqns(StdVector,EntityIterator,UInt):  Supplied an iterator which is not supported by FeSpace");
    }
  }

  //! get a nodal equation number
  UInt FeSpaceL2::GetNodeEqn(UInt nodeNr, UInt dof){
    EXCEPTION("NODE EQN NUMBERS NOT AVAILABLE FOR L2Space");
    //return nodeMap_[nodeNr][dof];
    return 0;
  }

  void FeSpaceL2::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem){
    FeSpaceH1::GetElemEqns(eqns,elem);
  }

  //! Get Equation numbers for a specific element
  void FeSpaceL2::GetElemEqns(StdVector<Integer>& eqns,
                              const Elem* elem,
                              UInt dof){
    FeSpaceH1::GetElemEqns(eqns,elem,dof);
  }

  void FeSpaceL2::GetOlasMappings( shared_ptr<SolStrategy> solStrat,
                                   StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                   std::map<UInt,StdVector<std::set<Integer> > >&
                                   minorBlocks ) {
    WARN("GetOlasMappings never testet with L2 space!!!!!!");
    FeSpaceH1::GetOlasMappings(solStrat,sbmBlocks,minorBlocks);
  }
  //! Map Nodal BC Equation NUmbers
  void FeSpaceL2::MapNodalBCs(){
    WARN("MapNodalBCs never testet with L2 space!!!!!!");
    FeSpaceH1::MapNodalBCs();
  }

  //! Map Nodal Equation Numbers
  void FeSpaceL2::MapNodalEqns(UInt phase){
    WARN("MapNodalEqns never testet with L2 space!!!!!!");
    FeSpaceH1::MapNodalEqns(phase);
  }

  void FeSpaceL2::PrintEqnMap(){

  }

}
