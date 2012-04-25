#include "FeSpaceH1.hh"
#include "H1Elems.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"

namespace CoupledField {

// declare class specific logging stream
 DECLARE_LOG(feSpaceH1)
 DEFINE_LOG(feSpaceH1, "feSpaceH1")

  //! Constructor
  FeSpaceH1::FeSpaceH1(PtrParamNode aNode, PtrParamNode infoNode,
                       Grid* ptGrid) 
  : FeSpace(aNode, infoNode, ptGrid  ) {
    type_ = H1;
  }
  
  FeSpaceH1::~FeSpaceH1(){
  }

  void FeSpaceH1::AddFeFunction( shared_ptr<BaseFeFunction> fct ){
    feFunction_ = fct;
    return;
  }
  
  

  //! Returns the number of (vectorial) unkowns on the element
  UInt FeSpaceH1::GetNumFunctions( const EntityIterator ent ){
    BaseFE * fe = GetFe(ent);
    return fe->GetNumFncs();
  }

  //! Return equation numbers for a all DOFs
  void FeSpaceH1::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){

    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = GetNumDofs();

    //This if clause should be avoided if the functionality for higher order entities
    //is implemented
    //First cover the nodal Case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(dofsPerUnknown);
      eqns.Init();
      if (gridToVirtualNodes_.find(node) != gridToVirtualNodes_.end()){
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[iDof] = nodeMap_[gridToVirtualNodes_[node][0]][iDof];
        }
      } else {
        // In case a node was not found, we reset the eqnarray to
        // size 0
        eqns.Resize(0);
      }
        

    } else if( ent.GetType() == EntityList::ELEM_LIST ||
        ent.GetType() == EntityList::SURF_ELEM_LIST||
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
    } else {
      EXCEPTION("In FeSpaceH1::GetEqns():  Supplied an iterator which is not supported by FeSpace");
    }
  }

  //! Return equation numbers for a specific dof
  void FeSpaceH1::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                        , UInt dof ){ 
    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
//    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();

    //First cover the nodal/grid case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(1);
      eqns.Init();
      //without the if clause we would have a segfault in case of higher order approximaton
      //for quadratic meshes
      if (gridToVirtualNodes_.find(node)!= gridToVirtualNodes_.end()) {
    	  eqns[0] = nodeMap_[gridToVirtualNodes_[node][0]][dof];
      } else {
        // In case a node was not found, we reset the eqnarray to
        // size 0
        eqns.Resize(0);
      }
    } else if( ent.GetType() == EntityList::ELEM_LIST ||
               ent.GetType() == EntityList::SURF_ELEM_LIST||
             ent.GetType() == EntityList::NC_ELEM_LIST){
      StdVector<UInt> nodes;
      GetNodesOfElement(nodes,ent.GetElem());
      eqns.Resize( nodes.GetSize() );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqns[iNode] = nodeMap_[nodes[iNode]][dof];
      }
    } else {
      EXCEPTION("In FeSpaceH1::GetEqns(StdVector,EntityIterator,UInt):  Supplied an iterator which is not supported by FeSpace");
    }
  
  }
  void FeSpaceH1::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                          , UInt dof, BaseFE::EntityType entityType){ 
      //Get result for the feFunction
      shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

      //Get "dimension" of one Unknown
//      UInt dofsPerUnknown = feFctResult->dofNames.GetSize();

      //First cover the nodal/grid case
      if ( ent.GetType() == EntityList::NODE_LIST ) {
        UInt node = ent.GetNode();
        eqns.Resize(1);
        eqns.Init();
        if(gridToVirtualNodes_.find(node)!= gridToVirtualNodes_.end()) {
          eqns[0] = nodeMap_[gridToVirtualNodes_[node][0]][dof];
        } else {
          // In case a node was not found, we reset the eqnarray to
          // size 0
          eqns.Resize(0);
        }
      } else if( ent.GetType() == EntityList::ELEM_LIST ||
                 ent.GetType() == EntityList::SURF_ELEM_LIST||
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
  UInt FeSpaceH1::GetNodeEqn(UInt nodeNr, UInt dof){
    return nodeMap_[nodeNr][dof];
  }

  //! Get Equation numbers for a specific element
  void FeSpaceH1::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem){
    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = GetNumDofs();

    StdVector<UInt> nodes;
    GetNodesOfElement(nodes,elem);
    eqns.Resize( nodes.GetSize() * dofsPerUnknown );
    eqns.Init();
    for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
      for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
        eqns[(iNode*dofsPerUnknown) + iDof] = 
            nodeMap_[nodes[iNode]][iDof];
      }
    }
  }

  //! Get Equation numbers for a specific element
  void FeSpaceH1::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof){
  }

  //! Map Nodal BC Equation Numbers
  void FeSpaceH1::MapNodalBCs(){
  LOG_TRACE(feSpaceH1) << "Mapping Nodal BCs";
    StdVector<UInt> actNodes;

    // check if feFunction was defined
    if( !feFunction_ ) {
      EXCEPTION("No FeFunction set at FeSpace");
    }
    
    // check if feFunction has a result assigned
    if( !feFunction_->GetResultInfo()) {
      EXCEPTION("No resultinfo associated with feFunction");
    }
    
    //Get Grip of HdBC List for the fefunction
    const HdBcList hdbcs = feFunction_->GetHomDirichletBCs();
    HdBcList::const_iterator actHBC;
    UInt dofsPerUnknown = GetNumDofs();

    for(actHBC = hdbcs.Begin(); actHBC != hdbcs.End(); actHBC++) {
      // Get EntityIterator
      GetNodesOfEntities(actNodes,(*actHBC)->entities);
      for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
        UInt vNode = actNodes[iNode];
        
         if( nodeMap_.BcKeys.find(vNode) == nodeMap_.BcKeys.end()){
           //nodeMap_.BcKeys[node] = StdVector<BaseFeFunction::BcType>(dofsPerUnknown,BaseFeFunction::NOBC);
           nodeMap_.BcKeys[vNode] = StdVector<BcType>(dofsPerUnknown);
           nodeMap_.BcKeys[vNode].Init(NOBC);
         }
         LOG_DBG(feSpaceH1) << "\tHDBC for vNode " << vNode << ", dof " <<  (*actHBC)->dof;
         nodeMap_.BcKeys[vNode][(*actHBC)->dof] = HDBC;
         bcCounter_[HDBC]++;
      }
    }

    
    // ===================================================================
    // DISTINCTION REGARDING HIERARCHICAL SPACE
    // ===================================================================
    
    if( !isHierarchical_ ) {

      //Get Grip of IdBC List for the fefunction
      const IdBcList idbcs = feFunction_->GetInHomDirichletBCs();
      IdBcList::const_iterator actIBC;

      for(actIBC = idbcs.Begin(); actIBC != idbcs.End(); actIBC++) {
        // Get all (Virtual) Nodes of the list
        GetNodesOfEntities(actNodes,(*actIBC)->entities);
        for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
          UInt vNode = actNodes[iNode];
          if( nodeMap_.BcKeys.find(vNode) == nodeMap_.BcKeys.end()){
            nodeMap_.BcKeys[vNode] = StdVector<BcType>(dofsPerUnknown);
            nodeMap_.BcKeys[vNode].Init(NOBC);
          }
          // check first, if this node was already processed
          if( nodeMap_.BcKeys[vNode][(*actIBC)->dof] != IDBC) {
            nodeMap_.BcKeys[vNode][(*actIBC)->dof] = IDBC;
            bcCounter_[IDBC]++;
          }
        }
      }

      //Get Grip of constraint List for the fefunction
      const ConstraintList constraints = feFunction_->GetConstraints();
      ConstraintList::const_iterator actConstr;
      for(actConstr = constraints.Begin(); actConstr != constraints.End(); actConstr++) {
        StdVector<UInt> slaveNodes;
        GetNodesOfEntities(slaveNodes,(*actConstr)->slaveEntities);
        UInt masterDof = (*actConstr)->masterDof;
        UInt slaveDof = (*actConstr)->slaveDof;
        UInt mNode = slaveNodes[0];

        for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
          if( nodeMap_.BcKeys.find(slaveNodes[iNode]) == nodeMap_.BcKeys.end()){
            nodeMap_.BcKeys[slaveNodes[iNode]] = StdVector<BcType>(dofsPerUnknown);
            nodeMap_.BcKeys[slaveNodes[iNode]].Init(NOBC);
          }
          nodeMap_.BcKeys[slaveNodes[iNode]][slaveDof] = CS;
          nodeMap_.constraintNodes[std::pair<Integer,Integer>(slaveNodes[iNode],slaveDof)] = 
              std::pair<Integer,Integer>(mNode,masterDof);
          bcCounter_[CS]++;
        }
      }
    } // is !hierarchical 
    
    if( isHierarchical_) {
      
      /* 
       * In a hierarchical space, only the equations corresponding to a vertex
       * get fixed. The edge / face / inner equations are set to zero.
       */
      const IdBcList idbcs = feFunction_->GetInHomDirichletBCs();
      IdBcList::const_iterator actIBC;
      for(actIBC = idbcs.Begin(); actIBC != idbcs.End(); actIBC++) {
        // Get all (Virtual) Nodes of the list
        GetNodesOfEntities(actNodes,(*actIBC)->entities);
        for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){

          UInt vNode = actNodes[iNode];
          //TODO find the source
          //make it zero based
          if( nodeMap_.BcKeys.find(vNode) == nodeMap_.BcKeys.end()){
            nodeMap_.BcKeys[vNode] = StdVector<BcType>(dofsPerUnknown);
            nodeMap_.BcKeys[vNode].Init(NOBC);
          }

          // check, if this node belongs to a vertex
          // yes set "true" IDBC
          //  no: set HDBC
          if( nodesType_[vNode] == BaseFE::VERTEX ) {

            // check first, if this node was already processed
            if( nodeMap_.BcKeys[vNode][(*actIBC)->dof] != IDBC) {
              nodeMap_.BcKeys[vNode][(*actIBC)->dof] = IDBC;
              bcCounter_[IDBC]++;
            }
          } else {
            nodeMap_.BcKeys[vNode][(*actIBC)->dof] = HDBC ;
          }
        } // loop over virtual nodes
      } // loop over idbcs

      //Get Grip of constraint List for the fefunction
      const ConstraintList constraints = feFunction_->GetConstraints();
      ConstraintList::const_iterator actConstr;
      for(actConstr = constraints.Begin(); actConstr != constraints.End(); actConstr++) {
        StdVector<UInt> slaveNodes;
        StdVector<UInt> slaveVertexNodes;
        GetNodesOfEntities(slaveNodes,(*actConstr)->slaveEntities);
        GetNodesOfEntities(slaveVertexNodes,(*actConstr)->slaveEntities, 
                           BaseFE::VERTEX);
        UInt masterDof = (*actConstr)->masterDof;
        UInt slaveDof = (*actConstr)->slaveDof;
        UInt mNode = slaveNodes[0];

        for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
          UInt slaveNode = slaveNodes[iNode];
          if( nodeMap_.BcKeys.find(slaveNode) == nodeMap_.BcKeys.end()){
            nodeMap_.BcKeys[slaveNode] = StdVector<BcType>(dofsPerUnknown);
            nodeMap_.BcKeys[slaveNode].Init(NOBC);
          }

          // check, if this node belongs to a vertex
          // yes set "true" CS
          //  no: set HDBC
          if( slaveVertexNodes.Find( slaveNode ) != -1 ) {
            nodeMap_.BcKeys[slaveNode][slaveDof] = CS;
            nodeMap_.constraintNodes[std::pair<Integer,Integer>(slaveNode,slaveDof)] = 
                std::pair<Integer,Integer>(mNode,masterDof);
            bcCounter_[CS]++; 
          } else {
            nodeMap_.BcKeys[slaveNode][slaveDof] = HDBC;
          }
        } // loop over nodes
      } // loop over constraints
    }

    //DEBUG output reaenable along with logging
    //std::map< Integer,StdVector<FeSpace::BcType> >::iterator iter = nodeMap_.BcKeys.begin();
    //for(iter; iter!= nodeMap_.BcKeys.end();iter++){
    //  std::cout << "The node #" << iter->first << " has the following flags: " << std::endl;
    //  for(UInt i = 0; i < iter->second.GetSize() ; i++){
    //    std::cout <<  iter->second[i] << std::endl;
    //  }
    //  std::cout << std::endl;
    //}
  }
  //=======================================================
  //Perform Nodal equation Numbering
  //=======================================================
  void FeSpaceH1::MapNodalEqns(UInt phase){
    LOG_TRACE(feSpaceH1) << "Mapping nodal Eqns, Phase " << phase;
    UInt dofsPerUnknown = 0;
    shared_ptr<ResultInfo> feFctResult;
    StdVector< shared_ptr<EntityList> > fctEntList;
    StdVector< shared_ptr<EntityList> >::iterator entIt;
    boost::unordered_map< Integer , StdVector<BcType> >::iterator bcIt;
    UInt actNode = 0;
    NodeTypeMap::const_iterator it;
    switch(phase){
      case 1:
        // number the nodal equations if they are not contained in the BcKeys
        // due to our virtual node array we have a valid structure

        //Get result for the feFunction
        feFctResult = feFunction_->GetResultInfo();
        
        // Get number of dofs
        dofsPerUnknown = GetNumDofs();
        
        //now loop over all nodes and assign an equation number
         it= nodesType_.begin();
        for( ; it != nodesType_.end(); ++it ) {
          actNode = it->first;
          if(nodeMap_.eqns.find(actNode) == nodeMap_.eqns.end()){
            nodeMap_[actNode] = StdVector<Integer>(dofsPerUnknown);
            nodeMap_[actNode].Init(-1);
          }

          for(UInt iDof = 0; iDof < dofsPerUnknown; iDof ++){
            if(nodeMap_.BcKeys.find(actNode) != nodeMap_.BcKeys.end()){
              if(nodeMap_.BcKeys[actNode][iDof] == NOBC){
                nodeMap_[actNode][iDof] = ++numEqns_;
                numFreeEquations_++;
              }
            }else if(nodeMap_[actNode][iDof] == -1){
              nodeMap_[actNode][iDof] = ++numEqns_;
              numFreeEquations_++;
            }
          }
        }
        
        break;
        
      case 2:
        bcIt = nodeMap_.BcKeys.begin();
        while(bcIt != nodeMap_.BcKeys.end()){
          //TODO: take the kind of BC into account
          //OPTIONAL: Delete the BcKeys Map afterwards?!
          UInt vNode = bcIt->first;
          StdVector<BcType> & dofs = bcIt->second;
          
          LOG_DBG3(feSpaceH1) << "\tvNode: " << vNode;
          LOG_DBG3(feSpaceH1) << "\teqns: " << nodeMap_[vNode].ToString();
          for(UInt iDof =0; iDof < dofs.GetSize(); iDof++){
            if(dofs[iDof]!=NOBC){
               //nodeMap_[bcIt->first] = ++numEqns_;
              if(dofs[iDof] == IDBC){
                nodeMap_[vNode][iDof] = ++numEqns_;
              } else if(dofs[iDof] == HDBC){
                nodeMap_[vNode][iDof] = 0 ;
              } 
              //else if(bcIt->second[iDof] == CONSTRAINT){
              //  Integer masterNode = nodeMap_.slaveMasterNodes[bcIt->first];
              //  nodeMap_[bcIt->first] = -1*masterNode;
              //}
            }
          }
          bcIt++;
        }
        {
          std::map<std::pair<Integer,Integer>, std::pair<Integer,Integer>  >::iterator  conNodes = nodeMap_.constraintNodes.begin();
          for(; conNodes !=nodeMap_.constraintNodes.end() ; conNodes++){

            std::pair<Integer,Integer> slaveP = conNodes->first;
            std::pair<Integer,Integer> masterP = conNodes->second;
            
            // attention: we can not rely on the fact, that the nodes mentioned in the
            // constraint nodes list are really contained in our region list
            if( nodeMap_.eqns.find(masterP.first) == nodeMap_.eqns.end() )  {
              EXCEPTION("Master node of constraints could not be found");
            }
            
            if( nodeMap_.eqns.find(slaveP.first) == nodeMap_.eqns.end() )  {
              EXCEPTION("Slave node of constraints could not be found");
            }

            nodeMap_[slaveP.first][slaveP.second] = -1*nodeMap_[masterP.first][masterP.second];
          }
        }
        break;
      default:
        EXCEPTION("In FeSpaceH1::MapNodalEqns(): Supplied wrong argument for the numbering phase");
        break;
    }
  }

  
  void FeSpaceH1::GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                   std::map<UInt,StdVector<std::set<Integer> > >&
                                   minorBlocks ) {
    
    // currently we only support "standard" solution strategy
    if( solStrat_->GetType() != SolStrategy::STD_STRATEGY ) {
      EXCEPTION( "Currently we just support the standard solution strategy for H1.");
    }
    
    FeFctIdType fctId = feFunction_->GetFctId();
    
    // Check, if static condensation is to be performed
    bool statCond = solStrat_->UseStaticCondensation();
    if( !statCond ) {
      // -------------------------
      //  No Static Condensation
      // -------------------------
      // we have just one big matrix block, where all equations are put in
      sbmBlocks.Resize(1);
      std::set<Integer> & allEqns = sbmBlocks[0][fctId];
      for( UInt i = 0; i < numEqns_; ++i ) {
        allEqns.insert(i+1);
      }
    } else {

      // ---------------------
      //  Static Condensation
      // ---------------------
      // Create 2 SBM block sets:
      // 1st block: Contains all nodal, edge and face contributions
      // 2nd block: contains all element interior equations (in 2D: faces)
      
      // push back empty set -> sbm block 0 has no minor blocks
      minorBlocks[1]=(StdVector<std::set<Integer> >());
      
      // Preliminary first set: all equations
      // If we have lateron no interior equations, we can directly return
      // this set
      std::set<Integer>  nonIntEqns;
      for( UInt i = 0; i < numEqns_; ++i ) {
        nonIntEqns.insert(i+1);
      }
      
      // Determine type of entity type to condensate
      // 2D: faces, 3D: interior
      BaseFE::EntityType intType;
      if( ptGrid_->GetDim() == 3) {
        intType = BaseFE::INTERIOR;
      } else {
        intType = BaseFE::FACE;
      }

      // list with all equation numbers for interior block
      std::set<Integer> intBlock;

      // Loop over all elements
      boost::unordered_map< UInt, ElemVirtualNodes >::iterator elemIt = virtualNodes_.begin();
      for( ; elemIt != virtualNodes_.end(); ++elemIt ) {

        const UInt elemNum = elemIt->first;
        const Elem * elem = ptGrid_->GetElem(elemNum);
        UInt dim = Elem::shapes[elem->type].dim;;
        LOG_DBG(feSpaceH1) << "\nDim of elem #" 
            << elemNum << ": " << dim << std::endl;

        
        ElemVirtualNodes& vn = elemIt->second;
        // collect all inner nodes 
        StdVector<UInt> & innerNodes = vn[intType].vNodes;
        LOG_DBG(feSpaceH1) << "innerNode has size " << innerNodes.GetSize() 
                                                             << std::endl;
        std::set<Integer> innerEqns;
        for(UInt i = 0; i < innerNodes.GetSize(); ++i ) {
          // collect entries for SBM block
          intBlock.insert( nodeMap_[innerNodes[i]].Begin(),
                           nodeMap_[innerNodes[i]].End() );
          // collect entries for minor block
          innerEqns.insert(nodeMap_[innerNodes[i]].Begin(),
                           nodeMap_[innerNodes[i]].End());
        } // loop over interior nodes
        
        if( innerEqns.size() ) {
          minorBlocks[1].Push_back(innerEqns);
        }
      } // loop over elements
      
      
      
      // take set-difference:
      // so far we have collected all inner equations, now subtract them
      // from the set of "all" equations, to get non-inner equations
      std::set<Integer> diff;
      std::insert_iterator<std::set<Integer> > insert_it (diff, diff.begin());
      std::set_difference( nonIntEqns.begin(), nonIntEqns.end(),
                           intBlock.begin(), intBlock.end(),
                           insert_it );
      sbmBlocks.Resize(2);
      sbmBlocks[0][fctId] = diff;
      sbmBlocks[1][fctId] = intBlock;
//      }
      
      
    } // if static condensation
    
    
    
  }

  void FeSpaceH1::PrintEqnMap(){

    // obtain (fctId,eqnNr) -> (sbm,index) mapping from OLAS
    StdVector<UInt> blockNums, indices;
    AlgebraicSys * algSys = feFunction_->GetSystem();
    FeFctIdType fctId = feFunction_->GetFctId();
    algSys->MapCompleteFctIdToIndex(fctId, blockNums, indices);
    
    std::string resultName = 
        SolutionTypeEnum.ToString(feFunction_->GetResultInfo()->resultType);
    
    std::cout << " *****************************************\n";
    std::cout << "  F E - S P A C E - I N F O R M A T I O N \n";
    std::cout << " *****************************************\n";
    std::cout << " Physical Quantity: " << resultName << std::endl;
    std::cout << " FeFunction Id: " << fctId << std::endl << std::endl;
    
    
    // =================================
    // 1) ELEMENT INFORMATION
    // =================================
    
    // iterate over all elements
    Grid* ptGrid = feFunction_->GetGrid();
    boost::unordered_map< UInt, ElemVirtualNodes >::iterator elemIt;
    
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

          // check, if any virtual nodes are assigned at all
          if( offset[i] == 0 ) {
            std::cout << "-\t-\t-\t-\n";
          }
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
            // immediately begin new line, if entity has no equations
            if( eqns.GetSize() == 0 ){
              std::cout << "\n";
            }
            
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
    boost::unordered_map< Integer , StdVector<Integer> >::const_iterator nodeIt = nodeMap_.eqns.begin();
    boost::unordered_map< Integer , StdVector<BcType> >::iterator nodeBcIt;

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
        const Integer & eqn = nodeIt->second[iDof];
        std::cout << "|\t" << eqn;


        if( eqn == 0) {
          std::cout << "\t|" << std::setw(1) << "-";; 

          // index
          std::cout << "\t|" << std::setw(8) << "-";
        } else {
          // sbm-block  
          std::cout << "\t|" << std::setw(1) << blockNums[std::abs(eqn)-1]; 

          // index
          std::cout << "\t|" << std::setw(8) << indices[std::abs(eqn)-1];
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
  
  UInt FeSpaceH1::GetNumDofs() const {
    // Just return number of components
    assert(feFunction_ );
    return feFunction_->GetResultInfo()->dofNames.GetSize();
  }

 

}
