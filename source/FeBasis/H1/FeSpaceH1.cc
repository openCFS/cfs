#include "FeSpaceH1.hh"
#include "H1Elems.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/algsys/AlgebraicSys.hh"

namespace CoupledField {

// declare class specific logging stream
 DECLARE_LOG(feSpaceH1)
 DEFINE_LOG(feSpaceH1, "feSpaceH1")

  //! Constructor
  FeSpaceH1::FeSpaceH1(PtrParamNode aNode, PtrParamNode infoNode) 
  : FeSpace(aNode, infoNode ) {
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
    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();

    //This if clause should be avoided if the functionality for higher order entities
    //is implemented
    //First cover the nodal Case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(dofsPerUnknown);
      eqns.Init();
      for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
        eqns[iDof] = nodeMap_[gridToVirtualNodes_[node]][iDof];
      }
    } else if( ent.GetType() == EntityList::ELEM_LIST ||
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
      eqns[0] = nodeMap_[gridToVirtualNodes_[node]][dof];

    } else if( ent.GetType() == EntityList::ELEM_LIST ||
               ent.GetType() == EntityList::SURF_ELEM_LIST){
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
        eqns[0] = nodeMap_[gridToVirtualNodes_[node]][dof];

      } else if( ent.GetType() == EntityList::ELEM_LIST ||
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
  UInt FeSpaceH1::GetNodeEqn(UInt nodeNr, UInt dof){
    return nodeMap_[nodeNr][dof];
  }

  //! Get Equation numbers for a specific element
  void FeSpaceH1::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem){
    //Get result for the feFunction
    shared_ptr<ResultInfo> feFctResult = feFunction_->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = feFctResult->dofNames.GetSize();

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

  //! Map Nodal BC Equation NUmbers
  void FeSpaceH1::MapNodalBCs(){
    StdVector<UInt> actNodes;
    
    //Get Grip of HdBC List for the fefunction
    const HdBcList hdbcs = feFunction_->GetHomDirichletBCs();
    HdBcList::const_iterator actHBC;
    UInt dofsPerUnknown = feFunction_->GetResultInfo()->dofNames.GetSize();

    for(actHBC = hdbcs.Begin(); actHBC != hdbcs.End(); actHBC++) {
      // Get EntityIterator
      GetNodesOfEntities(actNodes,(*actHBC)->entities);
      for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
         if( nodeMap_.BcKeys.find(actNodes[iNode]) == nodeMap_.BcKeys.end()){
           //nodeMap_.BcKeys[node] = StdVector<BaseFeFunction::BcType>(dofsPerUnknown,BaseFeFunction::NOBC);
           nodeMap_.BcKeys[actNodes[iNode]] = StdVector<BcType>(dofsPerUnknown);
           nodeMap_.BcKeys[actNodes[iNode]].Init(NOBC);
         }
         nodeMap_.BcKeys[actNodes[iNode]][(*actHBC)->dof] = HDBC;
         bcCounter_[HDBC]++;
      }
    }

    //Get Grip of IdBC List for the fefunction
    const IdBcList idbcs = feFunction_->GetInHomDirichletBCs();
    IdBcList::const_iterator actIBC;

    for(actIBC = idbcs.Begin(); actIBC != idbcs.End(); actIBC++) {
      // Get all (Virtual) Nodes of the list
      GetNodesOfEntities(actNodes,(*actIBC)->entities);
      for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
         //TODO find the source
         //make it zero based
         if( nodeMap_.BcKeys.find(actNodes[iNode]) == nodeMap_.BcKeys.end()){
           nodeMap_.BcKeys[actNodes[iNode]] = StdVector<BcType>(dofsPerUnknown);
           nodeMap_.BcKeys[actNodes[iNode]].Init(NOBC);
         }
         // check first, if this node was already processed
         if( nodeMap_.BcKeys[actNodes[iNode]][(*actIBC)->dof] != IDBC) {
           nodeMap_.BcKeys[actNodes[iNode]][(*actIBC)->dof] = IDBC;
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
        nodeMap_.constraintNodes[std::pair<Integer,Integer>(slaveNodes[iNode],slaveDof)] = std::pair<Integer,Integer>(mNode,masterDof);
        bcCounter_[CS]++;
      }
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
    UInt dofsPerUnknown = 0;
    shared_ptr<ResultInfo> feFctResult;
    StdVector< shared_ptr<EntityList> > fctEntList;
    StdVector< shared_ptr<EntityList> >::iterator entIt;
    std::map< Integer , StdVector<BcType> >::iterator bcIt;
    UInt actNode = 0;

    switch(phase){
      case 1:
        // number the nodal equations if they are not contained in the BcKeys
        // due to our virtual node array we have a valid structure

        //Get result for the feFunction
        feFctResult = feFunction_->GetResultInfo();
        // Get number of dofs
        dofsPerUnknown = feFctResult->dofNames.GetSize();
        
        //now loop over all nodes and assign an equation number
        for(UInt curNode = 0; curNode < nodes_.GetSize(); curNode++){
          actNode = nodes_[curNode];

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
          for(UInt iDof =0; iDof < bcIt->second.GetSize(); iDof++){
            if(bcIt->second[iDof]!=NOBC){
               //nodeMap_[bcIt->first] = ++numEqns_;
              if(bcIt->second[iDof] == IDBC){
                nodeMap_[bcIt->first] = ++numEqns_;
              } else if(bcIt->second[iDof] == HDBC){
                nodeMap_[bcIt->first] = 0 ;
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
            nodeMap_[slaveP.first][slaveP.second] = -1*nodeMap_[masterP.first][masterP.second];
          }
        }
        break;
      default:
        EXCEPTION("In FeSpaceH1::MapNodalEqns(): Supplied wrong argument for the numbering phase");
        break;
    }
  }


  void FeSpaceH1::PrintEqnMap(){
    

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
      std::cout << "\t\t#num\tgridNodes\tvNodes\tEqnNrs\tSBM-Block\tindex\n";
      std::cout << "\t\t=================================================================\n";

      std::string prefix = "\t\t\t\t\t ";
      
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
          entNumbers = ptElem->connect;
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
            Integer index = vNodes.Find(gridToVirtualNodes_[actNode]); 
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
          std::cout << "\t\t#" << std::abs(static_cast<Double>(entNumbers[i])) << "\t\t";

          // leave, virtual node numbers are assigned
          for( UInt j = 0; j < offset[i]; ++j ) {

            // print virtual node only for first entry            
            if( j > 0 ) {
              std::cout << prefix;
            } else {
              // physicalnode
              if( rNodes[pos] > 0 ) {
              std::cout << rNodes[pos] << "\t";
              } else {
                std::cout <<  "-\t ";
              }
            }
            // print virtual node
            std::cout << vNodes[pos] << "\t";
            

            // equation numbers (loop)
            StdVector<Integer> & eqns = nodeMap_[vNodes[pos]];
            for( UInt iEqn = 0; iEqn < eqns.GetSize(); ++iEqn ) {
              Integer & eqn = eqns[iEqn];

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

    //      //  VIRTUAL NODES 
    //      // ---------------
    //      std::cout << "a) Nodal mapping\n" 
    //      << "----------------\n";
    //      StdVector<UInt> & vVertexNodes = elemIt->second[BaseFE::VERTEX].vNodes;
    //      StdVector<UInt> & vEdgeNodes = elemIt->second[BaseFE::EDGE].vNodes;
    //      StdVector<UInt> & vFaceNodes = elemIt->second[BaseFE::FACE].vNodes;
    //      StdVector<UInt> & vInnerNodes = elemIt->second[BaseFE::INTERIOR].vNodes;
    //      
    //      StdVector<UInt> rVertexNodes(vVertexNodes.GetSize());
    //      StdVector<UInt> rEdgeNodes(vEdgeNodes.GetSize());
    //      StdVector<UInt> rFaceNodes(vFaceNodes.GetSize());
    //      StdVector<UInt> rInnerNodes(vInnerNodes.GetSize());
    //      rVertexNodes.Init(0);
    //      rEdgeNodes.Init(0);
    //      rFaceNodes.Init(0);
    //      rInnerNodes.Init(0);
    //      
    //
    //      // loop over nodal connectivity of element and try to find corresponding virtual node
    //      for( UInt i = 0; i < ptElem->connect.GetSize(); ++i ) {
    //        UInt actNode = ptElem->connect[i];
    //        
    //        // vertex nodes
    //        Integer index = vVertexNodes.Find(gridToVirtualNodes_[actNode]); 
    //        if( index > -1 ) {
    //          rVertexNodes[index] = actNode;
    //        }
    //        // edge nodes
    //        index = vEdgeNodes.Find(gridToVirtualNodes_[actNode]); 
    //        if( index > -1 ) {
    //          rEdgeNodes[index] = vEdgeNodes[index];
    //        }
    //        // face nodes
    //        index = vFaceNodes.Find(gridToVirtualNodes_[actNode]); 
    //        if( index > -1 ) {
    //          rFaceNodes[index] = vFaceNodes[index];
    //        }
    //        // inner nodes
    //        index = vInnerNodes.Find(gridToVirtualNodes_[actNode]); 
    //        if( index > -1 ) {
    //          rInnerNodes[index] = vInnerNodes[index];
    //        }
    //      }
    //      std::cout << "vVertexNodes: " << vVertexNodes.ToString(0) << std::endl;
    //      std::cout << "real nodes:   " << rVertexNodes.ToString(0) << std::endl << std::endl;
    //      
    //      std::cout << "vEdgeNodes:   " << vEdgeNodes.ToString(0) << std::endl;
    //      std::cout << "real nodes:   " << rEdgeNodes.ToString(0) << std::endl << std::endl;
    //      
    //      std::cout << "vFaceNodes:   " << vFaceNodes.ToString(0) << std::endl;
    //      std::cout << "real nodes:   " << rFaceNodes.ToString(0) << std::endl << std::endl;
    //      
    //      std::cout << "vInnerNodes:  " << vInnerNodes.ToString(0) << std::endl;
    //      std::cout << "real nodes:   " << rInnerNodes.ToString(0) << std::endl << std::endl;
    //    
    //    std::cout << std::endl;
    //    }  


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
          // type of node
          std::cout << "\t|  " 
              << BaseFE::entityType.ToString(nodesType_[nodeIt->first])[0];
        } else {
          std::cout << "       |\t";
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

 

}
