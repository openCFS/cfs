// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "FeSpaceNodal.hh"

#include "FeNodal.hh"
#include "H1/H1ElemsLagVar.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

// declare class specific logging stream
DECLARE_LOG(feSpaceNodal)
DEFINE_LOG(feSpaceNodal, "feSpaceNodal")

FeSpaceNodal ::FeSpaceNodal (PtrParamNode paramNode, PtrParamNode infoNode, 
                       Grid* ptGrid )
           : FeSpace( paramNode, infoNode, ptGrid ) {

  isHierarchical_ = false;

}

FeSpaceNodal ::~FeSpaceNodal() {

}

template<typename T>
void FeSpaceNodal::MapCoefFctToSpacePriv(shared_ptr<EntityList> entityList,
                                  shared_ptr<CoefFunction> coefFct,
                                  std::map <Integer, T>& vals,
                                  bool cache,
                                  const std::set<UInt>& comp ) {

  
  // Perform some checks at first:
  // a) if coefficient function is constant -> perform "easy mapping", without
  //    any caching. Here we simply take the constant value and 
  // b) if coefficient function has a general (spatial) dependency,
  //    we create solve FE problem on the subset of the FeSpace, as defined
  //    by the boundary condition, i.e. a mass-bilinearform is created and an
  //    according RHS integrator. The auxilliary data needed can be cached for
  //    repeated access / changing values (e.g. time / frequency dependend boundary+
  //    conditions.
  
  shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer


  if (IS_LOG_ENABLED(feSpaceNodal, trace)) {
    StdVector<UInt> compVec;
    std::set<UInt>::const_iterator it = comp.begin();
    for (; it != comp.end(); ++it ) {
      compVec.Push_back(*it);
    }
    LOG_TRACE(feSpaceNodal) << "Mapping coeffct " << coefFct->ToString() 
                            << " on " << entityList->GetName() << " for dofs "
                            << compVec.ToString() << " for FeFunction "
                            << SolutionTypeEnum.ToString(feFct->GetResultInfo()->resultType);
  }
  
  
  if( coefFct->GetDependency() == CoefFunction::CONSTANT ||
      coefFct->GetDependency() == CoefFunction::TIMEFREQ ) {
    // --------------------------
    //  SIMPLE MAPPING MECHANISM
    // --------------------------
    LocPointMapped lpm;
    StdVector<Integer> eqns, vEqns;
    Vector<T> dummyVec;
    if( coefFct->GetDimType() == CoefFunction::SCALAR ) {
      dummyVec.Resize(1);
      coefFct->GetScalar(dummyVec[0], lpm);
    } else {
      coefFct->GetVector( dummyVec, lpm);
    }
            
    // loop over all dofs
    std::set<UInt>::const_iterator it = comp.begin();
    for( ; it != comp.end(); ++it) {
    
      T val = dummyVec[*it];
      this->GetEntityListEqns( eqns, entityList, *it );
      UInt numEqns = eqns.GetSize();
      
      // --- non-hierachical case ---
      for( UInt i = 0; i < numEqns; ++i ) {
        vals[eqns[i]] = val; 
      }
      
    } // loop: dofs
  } else 
  {
    // ---------------------------
    //  GENERAL MAPPING MECHANISM
    // ---------------------------
    //Get coordinate representation of the entity list
    StdVector< StdVector<UInt> > idxMap;
    StdVector< Vector<Double> > globalCoords;
    StdVector< Vector<T> > valuesAtCoords;
    this->GetCoordinateRepresentation(entityList,idxMap,globalCoords);

    coefFct->GetVectorValuesAtCoords(globalCoords,valuesAtCoords);
    for( UInt aNode = 0; aNode < idxMap.GetSize(); ++aNode ) {
      StdVector<UInt> curEqs = idxMap[aNode];
      // loop over all dofs
      //if we have an empty set we loop over everything
      if(comp.size()==0){
        UInt numDof = valuesAtCoords[aNode].GetSize();
        for(UInt i = 0; i< numDof ; ++i) {
          vals[curEqs[i]] = valuesAtCoords[aNode][i];
        }
      }else{
        std::set<UInt>::const_iterator it = comp.begin();
        for( ; it != comp.end(); ++it) {
          vals[curEqs[*it]] = valuesAtCoords[aNode][*it];
        }
      }
    }
  } // end of general mapping section

}


void FeSpaceNodal::GetCoordinateRepresentation( shared_ptr<EntityList>   eList,
                                                   StdVector< StdVector<UInt> > & idxMap,
                                                   StdVector< Vector<Double> > & coords){

    shared_ptr<BaseFeFunction> feFct = feFunction_.lock();


   StdVector<UInt> nodeVec;
   //StdVector<Vector<Double> > coordinateVec;
   GetNodalCoords(coords,nodeVec,eList);

   idxMap.Resize(nodeVec.GetSize());
   idxMap.Init(StdVector<UInt>(feFct->GetResultInfo()->dofNames.GetSize()));

   StdVector<UInt>::iterator setIter = nodeVec.Begin();
   for(UInt pos = 0;setIter!= nodeVec.End();++setIter,++pos){
     StdVector<Integer> eqns;
     UInt curN = *setIter;
     eqns = GetVNodeEqns(curN);
     for(UInt i=0;i < eqns.GetSize();++i){
       if(abs(eqns[i]) > 0)
         idxMap[pos][i] = std::abs(eqns[i]);
     }
   }
 }

void FeSpaceNodal::GetNodalCoords(StdVector<Vector<Double> > & coords,
                                     StdVector<UInt> & vNodeNums,
                                     const shared_ptr<EntityList> & entities){
   //cover the special case of nodal mapping in which we can directly access the node coordinates
   std::cout << "\t \t Getting Nodal Coords ...." << std::endl;std::cout.flush();
   EntityIterator entIt = entities->GetIterator();
   StdVector<UInt> curNodes;
   StdVector<UInt> curENodes;
   std::map<UInt,bool> vnodeVisited;
   std::cout << "\t \t \t Getting nodes of entities ...." << std::endl;;std::cout.flush();


   if(mapType_ == GRID){
     StdVector<UInt> nodes;
     this->ptGrid_->GetNodesByName(nodes, entities->GetName());
     UInt size = nodes.GetSize();
     UInt dim = this->ptGrid_->GetDim();
     vNodeNums.Resize(size);
     coords.Resize(size,Vector<Double>(dim));
     for(UInt i=0;i<size;i++){
       vNodeNums[i] = this->gridToVirtualNodes_[nodes[i]][0];
       this->ptGrid_->GetNodeCoordinate(coords[i],nodes[i],true);
     }
     std::cout << "USED GIRD REPRESENTATION" << std::endl;
   }else{
     this->GetNodesOfEntities(curNodes,entities,BaseFE::ALL);
     vNodeNums.Clear();
     coords.Clear();
     vNodeNums.Reserve(curNodes.GetSize());
     coords.Reserve(curNodes.GetSize());
     LocPointMapped lpm;
     shared_ptr<ElemShapeMap> esm;
     for(entIt.Begin();!entIt.IsEnd();entIt++){
       //get elem virtual nodes
       this->GetNodesOfElement(curENodes,entIt.GetElem(),BaseFE::ALL);
       //get shape map
       esm = this->ptGrid_->GetElemShapeMap( entIt.GetElem(), true );


       Matrix<Double> nCoord;
       //make here decision...
       if(mapType_ == GRID){
         UInt rows = Elem::shapes[entIt.GetElem()->type].nodeCoords.GetSize();
         UInt cols = Elem::shapes[entIt.GetElem()->type].nodeCoords[0].GetSize();
         nCoord.Resize(rows,cols);
         for(UInt i=0;i<rows;i++){
           for(UInt j=0;j<cols;j++){
             nCoord[i][j] = Elem::shapes[entIt.GetElem()->type].nodeCoords[i][j];
           }
         }
       }else{
         //cast down to Lagrange Var
         FeNodal * myElem = dynamic_cast<FeNodal *>(this->GetFe(entIt));
         myElem->GetLocalDOFCoordinates(nCoord);
       }

       //loop over nodes and assign the LocalPointMapped
       if(nCoord.GetNumRows()!=curENodes.GetSize()){
         EXCEPTION("nodal coordinates do not match equation numbers.May not happen in grid case");
       }
       for(UInt i=0;i<nCoord.GetNumRows();++i){
         if(!vnodeVisited[curENodes[i]]){
           Vector<Double> locCoord(nCoord.GetNumCols());
           Vector<Double> globCoord;
           for(UInt d=0;d<nCoord.GetNumCols();d++){
             locCoord[d] = nCoord[i][d];
           }
           LocPoint lp = locCoord;
           lpm.Set(lp,esm,1.0);
           vNodeNums.Push_back(curENodes[i]);
           lpm.shapeMap->Local2Global(globCoord,locCoord);
           coords.Push_back(globCoord);
           vnodeVisited[curENodes[i]] = true;
         }
       }
       esm.reset();
     }
   }
   std::cout << "\t \t done" << std::endl;;std::cout.flush();
 }



}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template void FeSpaceNodal::MapCoefFctToSpacePriv<Double>( shared_ptr<EntityList> ,shared_ptr<CoefFunction>, std::map <Integer, Double>&,bool,const std::set<UInt>&);
  template void FeSpaceNodal::MapCoefFctToSpacePriv<Complex>(shared_ptr<EntityList> ,shared_ptr<CoefFunction>, std::map <Integer, Complex>&,bool,const std::set<UInt>&);
#endif
