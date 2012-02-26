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

  void FeSpaceL2::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){
    ///Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();
    if ( ent.GetType() == EntityList::NODE_LIST )
    {
      UInt node = ent.GetNode();
      if(gridToVirtualNodes_.find(node) == gridToVirtualNodes_.end())
        return;

      //UInt nDisNodes = gridToVirtualNodes_[node].GetSize();
      //eqns.Resize(dofsPerUnknown*nDisNodes);
      //eqns.Init();
      //TODO: Major Hack!!!!!!! this is just to enable compatibility
      // with the extract result scheme
      eqns.Resize(dofsPerUnknown);
      eqns.Init();
      //for (UInt iNode = 0; iNode < 1; iNode++ ) {
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[iDof] =
              nodeMap_[gridToVirtualNodes_[node][0]][iDof];
        }
      //}
    }else if( ent.GetType() == EntityList::ELEM_LIST ||
             ent.GetType() == EntityList::SURF_ELEM_LIST ||
             ent.GetType() == EntityList::NC_ELEM_LIST){
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

      //UInt nDisNodes = gridToVirtualNodes_[node].GetSize();
      //eqns.Resize(dofsPerUnknown*nDisNodes);
      //eqns.Init();
      //TODO: Major Hack!!!!!!! this is just to enable compatibility
      // with the extract result scheme
      eqns.Resize(1);
      eqns.Init();
      for (UInt iNode = 0; iNode < 1; iNode++ ) {
        eqns[iNode] = nodeMap_[gridToVirtualNodes_[node][iNode]][dof];
      }
    }else if( ent.GetType() == EntityList::ELEM_LIST ||
             ent.GetType() == EntityList::SURF_ELEM_LIST||
             ent.GetType() == EntityList::NC_ELEM_LIST){
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

  void FeSpaceL2::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent, UInt dof, BaseFE::EntityType entityType){
    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();
    if ( ent.GetType() == EntityList::NODE_LIST )
    {
      UInt node = ent.GetNode();
      if(gridToVirtualNodes_.find(node) == gridToVirtualNodes_.end())
        return;

      //UInt nDisNodes = gridToVirtualNodes_[node].GetSize();
      //eqns.Resize(dofsPerUnknown*nDisNodes);
      //eqns.Init();
      //TODO: Major Hack!!!!!!! this is just to enable compatibility
      // with the extract result scheme
      eqns.Resize(1);
      eqns.Init();
      for (UInt iNode = 0; iNode < 1; iNode++ ) {
        eqns[iNode] = nodeMap_[gridToVirtualNodes_[node][iNode]][dof];
      }
    }else if( ent.GetType() == EntityList::ELEM_LIST ||
              ent.GetType() == EntityList::SURF_ELEM_LIST ||
              ent.GetType() == EntityList::NC_ELEM_LIST){
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

  void FeSpaceL2::GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                   std::map<UInt,StdVector<std::set<Integer> > >&
                                   minorBlocks ) {
    FeSpaceH1::GetOlasMappings(sbmBlocks,minorBlocks);
  }
  //! Map Nodal BC Equation NUmbers
  void FeSpaceL2::MapNodalBCs(){
    FeSpaceH1::MapNodalBCs();
  }

  //! Map Nodal Equation Numbers
  void FeSpaceL2::MapNodalEqns(UInt phase){
    FeSpaceH1::MapNodalEqns(phase);
  }

  void FeSpaceL2::PrintEqnMap(){
    // obtain (fctId,eqnNr) -> (sbm,index) mapping from OLAS
        StdVector<UInt> blockNums, indices;
        AlgebraicSys * algSys = feFunction_->GetSystem();
        FeFctIdType fctId = feFunction_->GetFctId();
        algSys->MapCompleteFctIdToIndex(fctId, blockNums, indices);


        // =================================
        // 1) ELEMENT INFORMATION
        // =================================

        // iterate over all elements
        Grid* ptGrid = feFunction_->GetGrid();
        std::map< UInt, ElemVirtualNodes >::iterator elemIt;

        for( elemIt = virtualNodes_.begin();
            elemIt != virtualNodes_.end(); elemIt++ ) {

          // print element information (type, region, connect, edges, faces)
          const Elem * ptElem = ptGrid->GetElem(elemIt->first);

          std::cout << "=============\n"
                    << " Elem #" << elemIt->first << std::endl
                    << "=============\n";
          std::cout << "Type: " << Elem::feType.ToString( ptElem->type ) << std::endl;
          std::cout << "Connect: " << ptElem->connect.ToString( 0 ) << std::endl;

          // Print edge  information
          std::cout << "Edges: ";
          for( UInt i=0, numEdges = ptElem->edges.GetSize(); i < numEdges; ++i ) {
            Integer edgeNum = ptElem->edges[i];
            const Edge& myEdge = ptGrid->GetEdge(std::abs(edgeNum) );
            std::cout << "E #" << edgeNum << " ("
                      << myEdge.nodes[0] << "-> " << myEdge.nodes[1] << "), ";
          }
          std::cout << "\n";

          // Print face  information
          std::cout << "Faces: ";
          for( UInt i=0, numFaces = ptElem->faces.GetSize(); i < numFaces; ++i ) {
            Integer faceNum = ptElem->faces[i];
            const Face& myFace = ptGrid->GetFace(faceNum);
            std::cout << "F #" << faceNum << " (" << myFace.nodes.ToString( 0 ) << "), ";
          }
          std::cout << "\n\n";


          // print header
          std::cout << "\t\t#num\tvNodes\tEqnNrs\tSBM\tindex\n";
          std::cout << "\t\t=================================================================\n";

          std::string prefix = "\t\t\t";

          // print vertex / edge / face / inner information
          StdVector<BaseFE::EntityType> entTypes;
          entTypes = BaseFE::VERTEX, BaseFE::EDGE, BaseFE::FACE, BaseFE::INTERIOR;

          for( UInt iType = 0; iType < entTypes.GetSize(); ++iType ) {

            BaseFE::EntityType type = entTypes[iType];

            // vector containing entity numbers
            StdVector<UInt> entNumbers;
            StdVector<UInt> & vNodes = elemIt->second[type].vNodes;
            StdVector<UInt> & offset = elemIt->second[type].offset;
            StdVector<UInt> rNodes(vNodes.GetSize());
            rNodes.Init(0);

            // fill vector containing entity types
            if( type == BaseFE::VERTEX ) {
              // note: we may only take the number of nodes, which are
              // also used for the calculation element
              entNumbers = ptElem->connect;
              entNumbers.Resize(offset.GetSize());
            } else if( type == BaseFE::EDGE ) {
              entNumbers.Resize(ptElem->edges.GetSize());
              for( UInt i=0; i < ptElem->edges.GetSize(); ++i )
                entNumbers[i] = std::abs(ptElem->edges[i]);
            } else if( type == BaseFE::FACE ) {
              entNumbers.Resize(ptElem->faces.GetSize());
              for( UInt i=0; i < ptElem->faces.GetSize(); ++i )
                entNumbers[i] = ptElem->faces[i];
            } else if( type == BaseFE::INTERIOR ) {
              // only treat "interior", if we have any unknowns at all assigned
              // (= size of vNodes != 0)
              if( vNodes.GetSize() ) {
                entNumbers.Resize(1);
                entNumbers.Init(0);
              }
            }

            // try to find real nodes (currently only for vertices)
            if( type == BaseFE::VERTEX ) {
              for( UInt i = 0; i < vNodes.GetSize(); ++i ) {
                UInt actNode = ptElem->connect[i];
                Integer index = vNodes.Find(gridToVirtualNodes_[actNode][0]);
                if( index > -1 ) {
                  rNodes[index] = actNode;
                }
              }
            }
            // if any nodes are available
            if( vNodes.GetSize() ) {
              std::cout << iType+1 << ") " << BaseFE::entityType.ToString(type) << std::endl;
              std::cout << "========\n";
            } else {
              continue;
            }
            // loop over all entities
            UInt pos = 0;
            for( UInt i = 0; i < entNumbers.GetSize(); ++i ) {
              std::cout << "\t\t#" << std::abs(static_cast<Double>(entNumbers[i])) << "\t";

              // leave, virtual node numbers are assigned
              for( UInt j = 0; j < offset[i]; ++j ) {

                // print virtual node only for first entry
                if( j > 0 ) {
                  std::cout << prefix;
                }
                // print virtual node
                std::cout << vNodes[pos] << "\t";


                // equation numbers (loop)
                StdVector<Integer> & eqns = nodeMap_[vNodes[pos]];
                for( UInt iEqn = 0; iEqn < eqns.GetSize(); ++iEqn ) {
                  Integer & eqn = eqns[iEqn];

                  // indent succeeding equations correctly
                  if( iEqn > 0 ) {
                    std::cout << prefix << "\t";
                  }

                  //equation number
                  std::cout << eqn << "\t";

                  if( eqn > 0 ) {
                    // SBM-Block
                    std::cout << blockNums[eqn-1] << "\t";

                    // index
                    std::cout << indices[eqn-1] << "\n";
                  } else {
                    std::cout << "-\t-\n";
                  }
                }
                pos++;
              } // loop over virtual nodes
            } // loop over entity numbers
            std::cout << "\n";
          } // loop over entity types
          std::cout << "\n\n";
        } // loop over elements


        // =================================
        // 2) Global Equation numbering
        // =================================
        shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();
        // ---------------
        //  NODES
        // ---------------
        std::map< Integer , StdVector<Integer> >::iterator nodeIt = nodeMap_.eqns.begin();
        std::map< Integer , StdVector<BcType> >::iterator nodeBcIt;

        std::cout << "EQUATION MAPPING" << std::endl << std::endl;
        std::cout << "nodeNr \t|"  << " type  | " <<  std::setw (7)
        <<" Comp" << "|\teqnNr  \t| SBM\t|\tindex   |\t BC" << std::endl;
        std::cout << "----------------------------------------------------------------------------"
                  << std::endl;
        while(nodeIt != nodeMap_.eqns.end()){
          nodeBcIt = nodeMap_.BcKeys.find(nodeIt->first);
          for(UInt iDof =0; iDof < nodeIt->second.GetSize(); iDof++){
            // virtual node number (only once for all dofs) and type
            if( iDof == 0) {
              std::cout << nodeIt->first;

              // type of node (print first character of entityType )
              std::cout << "\t|  "
                  << BaseFE::entityType.ToString(nodesType_[nodeIt->first])[0];
            } else {
              std::cout << "        |";
            }

            // component
            std::cout << "\t|" << std::setw (8) << feFctResult->dofNames[iDof];
            // eqn number
            Integer & eqn = nodeIt->second[iDof];
            std::cout << "|\t" << eqn;


            if( eqn == 0) {
              std::cout << "\t|" << std::setw(1) << "-";;

              // index
              std::cout << "\t|" << std::setw(8) << "-";
            } else {
              // sbm-block
              std::cout << "\t|" << std::setw(1) << blockNums[eqn-1];

              // index
              std::cout << "\t|" << std::setw(8) << indices[eqn-1];
            }

            // bc type
            std::cout << "\t|\t";
            if(  nodeBcIt != nodeMap_.BcKeys.end() ) {
              if( nodeBcIt->second[iDof] != NOBC ) {
                std::cout << BcTypeEnum.ToString(nodeBcIt->second[iDof]);
              }
            }
            std::cout << std::endl;
          }
          nodeIt++;
        }
  }

  void FeSpaceL2::GetInteriorSurfaceElems(RegionIdType region,
                                    shared_ptr<EntityList> & surfElems,
                                    shared_ptr<EntityList> & opposingElems)
  {
    // what we do right now is just to create empty shared pointers
    // and fill them later in finalize
    if(interiorElemMap_.find(region) == interiorElemMap_.end()){
      //create a name for this list...
      std::string name1 = this->feFunction_->GetGrid()->GetRegion().ToString(region);
      std::string name2 = name1;
      name1 += "_dgInterior";
      name2 += "_dgInteriorOpposite";
      interiorElemMap_[region] = shared_ptr<NcElemList>(new NcElemList(this->feFunction_->GetGrid(),name1));
      interiorElemMapOpposite_[region] = shared_ptr<NcElemList>(new NcElemList(this->feFunction_->GetGrid(),name2));
    }
    surfElems = interiorElemMap_[region];
    opposingElems = interiorElemMapOpposite_[region];
  }

  void FeSpaceL2::GetInteriorSurfaceElems(RegionIdType region,
                                       shared_ptr<EntityList> & surfElems){
    // what we do right now is just to create empty shared pointers
    // and fill them later in finalize
    //create a name for this list...
    std::string name = this->feFunction_->GetGrid()->GetRegion().ToString(region);
    name += "_dgInterior";
    interiorElemMap_[region] = shared_ptr<NcElemList>(new NcElemList(this->feFunction_->GetGrid(),name));
    surfElems = interiorElemMap_[region];
  }

  //! \copydoc FeSpace::GetExteriorSurfaceElemsOfFeSpace
  void FeSpaceL2::GetExteriorSurfaceElems(RegionIdType region, shared_ptr<EntityList> & surfElems){

    if(exteriorElements_.find(region) == exteriorElements_.end()){
      std::string name = this->feFunction_->GetGrid()->GetRegion().ToString(region);
      name += "_dgExterior";
      exteriorElements_[region] = shared_ptr<NcElemList>(new NcElemList(this->feFunction_->GetGrid(),name));
    }
    surfElems = exteriorElements_[region];
  }

  void FeSpaceL2::CreateSurfaceElems(){
    //TODO: ADD DEBUGGING OUTPUT!!!!!
    //ok lets go. we collect the information from our maps and pass them to the grid
    std::set<RegionIdType> regions;
    std::map<RegionIdType,shared_ptr<NcElemList> >::iterator regIter;
    //its enough to go through the interiorElemMap_
    regIter = interiorElemMap_.begin();
    for(;regIter != interiorElemMap_.end();++regIter){
      regions.insert(regIter->first);
    }
    //create vectors
    StdVector<shared_ptr<NcSurfElem> > interiorElemList;
    StdVector<shared_ptr<NcSurfElem> > exteriorElemList;
    this->feFunction_->GetGrid()->GenerateDGSurfaceElemes(regions,interiorElemList,exteriorElemList);
    //ok first we fill the standard map and, if requested, the opposing map
    for(UInt aElem =0;aElem < interiorElemList.GetSize();++aElem){
      RegionIdType curReg = interiorElemList[aElem]->regionId;
      shared_ptr<NcElemList> curIntList = interiorElemMap_[curReg];
      curIntList->SetElement(interiorElemList[aElem]);
      //check if we need to fill in the opposite list
      if(interiorElemMapOpposite_.find(curReg) != interiorElemMapOpposite_.end()){
        shared_ptr<NcElemList> curOppList = interiorElemMapOpposite_[curReg];
        //we have at least one neighbor...
        curOppList->SetElement(interiorElemList[aElem]->neighbors[0]);
        UInt numNeighbors = interiorElemList[aElem]->neighbors.GetSize();
        for(UInt i = 1; i<numNeighbors;++i){
          //we need to push back in both lists....
          curOppList->SetElement(interiorElemList[aElem]->neighbors[i]);
          curIntList->SetElement(interiorElemList[aElem]);
        }
      }
    }
    //ok now we should have valid interior element lists lets fill the exterior
    for(UInt aElem =0;aElem < exteriorElemList.GetSize();++aElem){
      //obtain region id
      RegionIdType curReg = exteriorElemList[aElem]->regionId;
      shared_ptr<NcElemList> & curIntList = exteriorElements_[curReg];
      curIntList->SetElement(exteriorElemList[aElem]);
    }
    regIter = exteriorElements_.begin();
    for(;regIter != exteriorElements_.end();++regIter){
      this->feFunction_->AddEntityList(regIter->second);
    }

    regIter = interiorElemMap_.begin();
    for(;regIter != interiorElemMap_.end();++regIter){
      this->feFunction_->AddEntityList(regIter->second);
    }

  }
}
