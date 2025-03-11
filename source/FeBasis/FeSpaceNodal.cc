// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include "FeSpaceNodal.hh"



#include "FeNodal.hh"
#include "H1/H1ElemsLagVar.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "PDE/SinglePDE.hh"

namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(feSpaceNodal, "feSpaceNodal")

FeSpaceNodal ::FeSpaceNodal (PtrParamNode paramNode, PtrParamNode infoNode, 
                       Grid* ptGrid )
           : FeSpace( paramNode, infoNode, ptGrid ) {

  isHierarchical_ = false;

}

FeSpaceNodal ::~FeSpaceNodal() {

}

template<typename T>
void FeSpaceNodal::MapCoefFctToSpacePriv(StdVector<shared_ptr<EntityList> > entityLists,
                                         shared_ptr<CoefFunction> coefFct,
                                         shared_ptr<BaseFeFunction> feFct,
                                         std::map <Integer, T>& vals,
                                         bool cache,
                                         const std::set<UInt>& comp,
                                         bool updatedGeo) {

  
  // Perform some checks at first:
  // a) if coefficient function is constant -> perform "easy mapping", without
  //    any caching. Here we simply take the constant value and 
  // b) if coefficient function has a general (spatial) dependency,
  //    we create solve FE problem on the subset of the FeSpace, as defined
  //    by the boundary condition, i.e. a mass-bilinearform is created and an
  //    according RHS integrator. The auxilliary data needed can be cached for
  //    repeated access / changing values (e.g. time / frequency dependend boundary+
  //    conditions.
  
  if (IS_LOG_ENABLED(feSpaceNodal, trace)) {
    StdVector<UInt> compVec;
    std::set<UInt>::const_iterator it = comp.begin();
    for (; it != comp.end(); ++it ) {
      compVec.Push_back(*it);
    }
    std::string entityNames;
    for( UInt i = 0; i < entityLists.GetSize(); ++i ) {
      entityNames += entityLists[i]->GetName() + ", ";
    }
    LOG_TRACE(feSpaceNodal) << "Mapping coeffct " << coefFct->ToString() 
                            << " on " << entityNames << " for dofs "
                            << compVec.ToString() << " for FeFunction "
                            << ( feFct->GetResultInfo() ?
                               SolutionTypeEnum.ToString(feFct->GetResultInfo()->resultType)
                               : "");
  }
  
  if( coefFct->GetDependency() == CoefFunction::CONSTANT ||
      coefFct->GetDependency() == CoefFunction::TIMEFREQ ) {
    // --------------------------
    //  SIMPLE MAPPING MECHANISM
    // --------------------------
    LOG_DBG(feSpaceNodal) << "=> Applying simple mapping algorithm for constant values"; 
    LocPointMapped lpm;
    StdVector<Integer> eqns, vEqns;
    Vector<T> dummyVec;
    std::set<UInt> dofs(comp);
    
    if( coefFct->GetDimType() == CoefFunction::SCALAR ) {
      dummyVec.Resize(1);
      coefFct->GetScalar(dummyVec[0], lpm);
      if ( dofs.size() == 0) {
        dofs.insert(0);
      }
    } else {
      coefFct->GetVector( dummyVec, lpm);
      if ( dofs.size() == 0 ) {
        for ( UInt i=0, numDofs=coefFct->GetVecSize(); i<numDofs; ++i ) {
          dofs.insert(i);
        }
      }
    }

    // loop over all lists 
    for( UInt i = 0; i < entityLists.GetSize(); ++i ) {
      shared_ptr<EntityList> actList = entityLists[i];
      
      // loop over all dofs
      std::set<UInt>::const_iterator it = dofs.begin();
      for( ; it != dofs.end(); ++it) {

        T val = dummyVec[*it];
        this->GetEntityListEqns( eqns, actList, *it );
        UInt numEqns = eqns.GetSize();

        // --- non-hierarchical case ---
        for( UInt i = 0; i < numEqns; ++i ) {
          vals[eqns[i]] = val; 
        }
      } // loop: dofs
    } // loop: lists
  } else  {
    // ---------------------------
    //  GENERAL MAPPING MECHANISM
    // ---------------------------
    LOG_DBG(feSpaceNodal) << "=> Applying general mapping algorithm for general values";
    //Get coordinate representation of the entity list
    StdVector< StdVector<UInt> > idxMap;
    StdVector< Vector<Double> > globalCoords;
    StdVector< Vector<T> > valuesAtCoords;

    bool updateGeo;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock();
    assert(feFct);
    SinglePDE * curPDE = feFct->GetPDE();
    if (curPDE != NULL) {
      // the PDE exists, retrieve the updateGeometry-flag
      updateGeo = feFct->GetPDE()->IsUpdatedGeo();
    } else {
      // the PDE does not exist (can be the case if external data is applied), set the flag to false
      updateGeo = false;
    }
    bool realUpdatedGeo;
    if (updateGeo != updatedGeo) {
      WARN("Geometry update of FeFunction (" << updateGeo << ") and CoefFunction (" << updatedGeo
          << ") not consistent, falling back to non-updated geometry");
      realUpdatedGeo = false; // fall back plan
    } else {
      realUpdatedGeo = updatedGeo;
    }

    // loop over all lists 
    for( UInt i = 0; i < entityLists.GetSize(); ++i ) {
      shared_ptr<EntityList> actList = entityLists[i];
 
      // check if the coefFct lives on the upadted geometry and retrieve the coordinates accordingly
      //coefFct->feFunctions_;
      // get all coordinates of the nodes and their equation numbers
      this->GetCoordinateRepresentation(actList,idxMap,globalCoords,realUpdatedGeo);

      StdVector<shared_ptr<EntityList> > lists(1);
      lists[0] = actList; 

      // Only search for coordinates in elements belonging to the
      // regions where the current FeSpace is defined.
      LOG_DBG2(feSpaceNodal) << "Requesting values at globalCoords\n" << globalCoords.ToString();
      if( coefFct->GetDimType() == CoefFunction::SCALAR ) {
        StdVector<T> tmpVals;
        coefFct->GetScalarValuesAtCoords(globalCoords, tmpVals,
                                         ptGrid_, lists);
        valuesAtCoords.Resize(tmpVals.GetSize(), Vector<T>(1));
        for(UInt i=0;i<tmpVals.GetSize();i++){
          valuesAtCoords[i][0] = tmpVals[i];
        }
      } else {
        coefFct->GetVectorValuesAtCoords(globalCoords,valuesAtCoords,
                                         ptGrid_, lists, realUpdatedGeo);
      }
      //LOG_DBG2(feSpaceNodal) << "GetVectorValuesAtCoords for updatedGeo = " << updateGeo << ":\n" << valuesAtCoords.ToString();

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
        } // if: comp size
      } // loop: nodes
    } // loop: lists
  } // end of general mapping section

}


void FeSpaceNodal::GetCoordinateRepresentation( shared_ptr<EntityList>   eList,
                                                   StdVector< StdVector<UInt> > & idxMap,
                                                   StdVector< Vector<Double> > & coords,
                                                   bool updatedGeo){

    shared_ptr<BaseFeFunction> feFct = feFunction_.lock();


   StdVector<UInt> nodeVec;
   //StdVector<Vector<Double> > coordinateVec;
   GetNodalCoords(coords,nodeVec,eList,updatedGeo);

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
                                  const shared_ptr<EntityList> & entities,
                                  bool updatedGeo){
  LOG_TRACE(feSpaceNodal) << "Getting Nodal Coords ...."; 
  
  //cover the special case of nodal mapping in which we can directly access the node coordinates
  
   EntityIterator entIt = entities->GetIterator();
   StdVector<UInt> curNodes;
   StdVector<UInt> curENodes;
   std::map<UInt,bool> vnodeVisited;
   LOG_DBG(feSpaceNodal) << "Getting nodes of entities";


//   bool updateGeo;
//   shared_ptr<BaseFeFunction> feFct = feFunction_.lock();
//   assert(feFct);
//   updateGeo = feFct->GetPDE()->IsUpdatedGeo();

   if ( updatedGeo ) {
     LOG_DBG(feSpaceNodal) << "Using updated geometry";
   } else {
     LOG_DBG(feSpaceNodal) << "Using initial geometry (not updated)";
   }

   if(mapType_ == GRID){
     StdVector<UInt> nodes;
     this->ptGrid_->GetNodesByName(nodes, entities->GetName());
     UInt size = nodes.GetSize();
     UInt dim = this->ptGrid_->GetDim();
     vNodeNums.Resize(size);
     coords.Resize(size,Vector<Double>(dim));
     for(UInt i=0;i<size;i++){
       vNodeNums[i] = this->gridToVirtualNodes_[nodes[i]][0];
       //this->ptGrid_->GetNodeCoordinate(coords[i],nodes[i],true);
       this->ptGrid_->GetNodeCoordinate(coords[i],nodes[i],updatedGeo);
     }
     LOG_DBG(feSpaceNodal) << "USED GIRD REPRESENTATION";
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
       //esm = this->ptGrid_->GetElemShapeMap( entIt.GetElem(), true );
       esm = this->ptGrid_->GetElemShapeMap( entIt.GetElem(), updatedGeo );


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
   LOG_TRACE(feSpaceNodal) << "Finished Getting Nodal Coords ....";
 }



} // namespace

template void FeSpaceNodal::
MapCoefFctToSpacePriv<Double>( StdVector<shared_ptr<EntityList> > ,
                               shared_ptr<CoefFunction>,
                               shared_ptr<BaseFeFunction> feFct,
                               std::map <Integer, Double>&,
                               bool,const std::set<UInt>&,
                               bool updatedGeo);
template void FeSpaceNodal::
MapCoefFctToSpacePriv<Complex>( StdVector<shared_ptr<EntityList> > ,
                                shared_ptr<CoefFunction>,
                                shared_ptr<BaseFeFunction> feFct,
                                std::map <Integer, Complex>&,
                                bool,const std::set<UInt>&,
                                bool updatedGeo);
