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

  FeSpaceL2::FeSpaceL2(PtrParamNode aNode, PtrParamNode infoNode, Grid* ptGrid)
   : FeSpaceNodal(aNode, infoNode, ptGrid ) {
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
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();

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
    FeSpace::GetElemEqns(eqns,elem);
  }

  //! Get Equation numbers for a specific element
  void FeSpaceL2::GetElemEqns(StdVector<Integer>& eqns,
                              const Elem* elem,
                              UInt dof){
    FeSpace::GetElemEqns(eqns,elem,dof);
  }

  void FeSpaceL2::GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                   std::map<UInt,StdVector<std::set<Integer> > >&
                                   minorBlocks ) {
    FeSpace::GetOlasMappings(sbmBlocks,minorBlocks);
  }
  //! Map Nodal BC Equation NUmbers
  void FeSpaceL2::MapNodalBCs(){
    // Note: this a complete copy of thee FeSpaceH1Nodal::MapNodalBCs()

    StdVector<UInt> actNodes;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    // check if feFunction was defined
    if( !feFct ) {
      EXCEPTION("No FeFunction set at FeSpace");
    }

    // check if feFunction has a result assigned
    if( !feFct->GetResultInfo()) {
      EXCEPTION("No resultinfo associated with feFunction");
    }

    //Get Grip of HdBC List for the fefunction
    const HdBcList hdbcs = feFct->GetHomDirichletBCs();
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
        // loop over all dofs
        std::set<UInt>::const_iterator dofIt = (*actHBC)->dofs.begin();
        for( ; dofIt != (*actHBC)->dofs.end(); ++dofIt) { 
          nodeMap_.BcKeys[vNode][*dofIt] = HDBC;
          bcCounter_[HDBC]++;
        } // dofs
      }
    }

    //Get Grip of IdBC List for the fefunction
    const IdBcList idbcs = feFct->GetInHomDirichletBCs();
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
        // loop over all dofs
        std::set<UInt>::const_iterator dofIt = (*actIBC)->dofs.begin();
        for( ; dofIt != (*actIBC)->dofs.end(); ++dofIt) { 
          // check first, if this node was already processed
          if( nodeMap_.BcKeys[vNode][*dofIt] != IDBC) {
            nodeMap_.BcKeys[vNode][*dofIt] = IDBC;
            bcCounter_[IDBC]++;
          }
        } // dofs
      }
    }

    //Get Grip of constraint List for the fefunction
    const ConstraintList constraints = feFct->GetConstraints();
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

  //! Map Nodal Equation Numbers
  void FeSpaceL2::MapNodalEqns(UInt phase){
    FeSpace::MapNodalEqns(phase);
  }
  
   void FeSpaceL2::GetNodesOfElement( StdVector<UInt>& nodes,
                                    const Elem* ptElem,
                                    BaseFE::EntityType entType){
     UInt elemNum = ptElem->elemNum;
     if(virtualNodes_.find(elemNum) ==virtualNodes_.end()){

       EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element #"
           << ptElem->elemNum << " of region " 
           <<  ptGrid_->GetRegion().ToString(ptElem->regionId));
     }
     if(entType == BaseFE::ALL){
       nodes.Resize(virtualNodes_[elemNum][BaseFE::VERTEX].vNodes.GetSize()+
                    virtualNodes_[elemNum][BaseFE::EDGE].vNodes.GetSize()+
                    virtualNodes_[elemNum][BaseFE::FACE].vNodes.GetSize()+
                    virtualNodes_[elemNum][BaseFE::INTERIOR].vNodes.GetSize());
       ElemVirtualNodes::iterator nodeIt = virtualNodes_[elemNum].begin();
       UInt c = 0;
       while(nodeIt !=virtualNodes_[elemNum].end()){

         StdVector<UInt> & entNodes =  nodeIt->second.vNodes;
         for (UInt i = 0; i < entNodes.GetSize(); i++ ){
           nodes[c++] =  entNodes[i];
         }
         nodeIt++;
       }
     }else{
       nodes.Resize(virtualNodes_[elemNum][entType].vNodes.GetSize());
       const StdVector<UInt>& entNodes =  virtualNodes_[elemNum][entType].vNodes;
       for (UInt i = 0; i < entNodes.GetSize(); i++ ){
         nodes[i] =  entNodes[i];
       }
     }
   }
   
   void FeSpaceL2::CreatePolynomialNodes(){
     //follow the following algorithm
     // - loop over every element
     //  - get vertices and if not already mapped assign new virtual nodes
     //  - get edges and if not already mapped assign new virtual nodes
     //  - get faces and if not already mapped assign new virtual nodes
     //  - assign interior node to the element
     //  - now fill the Virtual node map with the information according to element orientation
     // - finally delete all intermediate arrays

     StdVector< shared_ptr<EntityList> > fctEntList;

     //ok these data structures are a bit messy but this is the most obvious typ
     // and they are temporary anyway. furthermore we run here once in a lifetime
     // 1st key: v/e/f number
     // 2nd key: element number
     // 3rc key: virtual nodes
     std::map< UInt, std::map<UInt, StdVector<Integer> > > vertexNodes;
     std::map< UInt, std::map<UInt, StdVector<Integer> > > edgenodes;
     std::map< UInt, std::map<UInt, StdVector<Integer> > > facenodes;
     std::map<UInt, StdVector<Integer> > interiornodes;
     //Different lists have to be treated differently
     EntityList::ListType actListType = EntityList::NO_LIST;


     shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
     assert(feFct);
     fctEntList = feFct->GetEntityList();

     //get the highest possible node number
     //UInt offset = feFunction_->GetGrid()->GetNumNodes();
     UInt offset =0;

     //Stores the current order
     // MappingType curMap = GRID;

     //changed algorithm first we add the volume elements and later on the surfaces
     // so we loop twice over the entities

     // loop over all entitylists (i.e. regions)
     for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){


       actListType = fctEntList[actList]->GetType();
       //determine mapping type and order of current entity list
       //if we do not find the name of the reigon in our map, we fall back to the default
       //BE CARFUL, if the user has specified something for the volume region but not for surface regions,
       // it could be that the fallback to default leads to errors!

       //check if we got what we expected

       if ( ! (actListType == EntityList::ELEM_LIST) &&
           ! (actListType == EntityList::SURF_ELEM_LIST) &&
           ! (actListType == EntityList::NC_ELEM_LIST))  {
         continue;
       }
       //cast down to element list
       EntityList* actElemList = fctEntList[actList].get();
       //      RegionIdType curReg = actElemList->GetRegion();

       // curMap = POLYNOMIAL;

       // Get iterator of current element list
       EntityIterator entIt = actElemList->GetIterator();

       // loop over all elements
       for(entIt.Begin(); !entIt.IsEnd();entIt++){

         // Fetch current finite element. This is performed by the specialized 
         // version, so this element "knows" already about its order, unknowns etc.
         BaseFE* ptFe = GetFe( entIt );
         const Elem* actEl = entIt.GetElem();  

         StdVector<UInt> permutations; // initially size 0
         UInt elemNum = actEl->elemNum;

         //===========================================================
         //Assign the BaseFE::VERTEX node numbers
         //===========================================================
         {
           UInt numVertexNodes = 0;
           UInt eNum = elemNum;
           UInt numVert = Elem::shapes[actEl->type].numVertices;
           StdVector<UInt> elemNodes = actEl->connect;

           EntityTypeNodes & vtn =  virtualNodes_[actEl->elemNum][BaseFE::VERTEX];

           // check, if the vertices of this element were already numbered
           if( vtn.vNodes.GetSize() == 0 ) {
             for ( UInt iVert= 0; iVert< numVert; iVert++ ) {
               UInt vertexNum = elemNodes[iVert];
               ptFe->GetNodalPermutation(permutations,actEl,BaseFE::VERTEX,iVert);
               numVertexNodes = permutations.GetSize();

               // Check if the vertex is already numbered.
               if( vertexNodes[vertexNum][eNum].GetSize() == 0 ) {

                 vertexNodes[vertexNum][eNum].Resize(numVertexNodes);
                 vertexNodes[vertexNum][eNum].Init();
                 for( UInt vertNode = 0; vertNode < numVertexNodes; ++vertNode ) {
                   vertexNodes[vertexNum][eNum][vertNode] = ++offset;
                   nodesType_[offset] = BaseFE::VERTEX;
                 }
               }


               for( UInt i = 0; i < numVertexNodes; ++i ) {
                 vtn.vNodes.Push_back(vertexNodes[vertexNum][eNum][permutations[i] ]);
               }
               vtn.offset.Push_back( permutations.GetSize() );

               gridToVirtualNodes_[vertexNum].Push_back(offset);
             } // loop over vertices
           } // mapping of vertex nodes
         }
         feFct->GetGrid()->MapEdges();
         feFct->GetGrid()->MapFaces();

         ElemShape actShape = Elem::shapes[actEl->type];

         //===========================================================
         //Assign the Edge node numbers
         //===========================================================
         {
           UInt numEdgeNodes = 0;
           UInt eNum = elemNum;
           EntityTypeNodes & etn =  virtualNodes_[actEl->elemNum][BaseFE::EDGE];

           // check if edges of this element were already numbered
           if( etn.vNodes.GetSize() == 0 ) {
             for ( UInt iEdge=0; iEdge < actShape.numEdges; iEdge++) {
               UInt edgeNum = std::abs(actEl->edges[iEdge]);
               //get the permutation Vector
               ptFe->GetNodalPermutation(permutations,actEl,BaseFE::EDGE,iEdge);
               numEdgeNodes = permutations.GetSize();
          
               // Check if the edge is already numbered.
               // Additionally, if we have the case of discontinuous approximation,
               // we number the nodes separately for every element anyway.
               if(edgenodes[edgeNum][eNum].GetSize() == 0 ) {
                 //here we assume spectral element approximation and we have
                 //order-1 nodes on the edge
                 edgenodes[edgeNum][eNum].Resize(numEdgeNodes);
                 edgenodes[edgeNum][eNum].Init();
                 for ( UInt edgeNode = 0;edgeNode < numEdgeNodes ;edgeNode++ ) {
                   edgenodes[edgeNum][eNum][edgeNode] = ++offset;
                   nodesType_[offset] = BaseFE::EDGE;
                 }
               }

               //fill the virtual Nodes in the correct ordering

               for ( UInt i = 0; i < numEdgeNodes ; i++ ) {
                 etn.vNodes.Push_back(edgenodes[edgeNum][eNum][ permutations[i] ]);
               }
               etn.offset.Push_back( permutations.GetSize() );
             } // loop over edges
           }
         }
         //===========================================================
         //Assign the Face node numbers
         //===========================================================
         {
           UInt numFaceNodes = 0;
           UInt eNum = elemNum;
           EntityTypeNodes & ftn =  virtualNodes_[actEl->elemNum][BaseFE::FACE];

           // check if faces of this element ware already numbered
           if( ftn.vNodes.GetSize() == 0 ) {
             for ( UInt iFace=0; iFace < actShape.numFaces; iFace++) {
               UInt faceNum = actEl->faces[iFace];
               //get the permutation Vector
               ptFe->GetNodalPermutation(permutations,actEl,BaseFE::FACE,iFace);
               numFaceNodes = permutations.GetSize();

               // Check if the face is already numbered.
               // Additionally, if we have the case of discontinuous approximation,
               // we number the nodes separately for every element separately anyway.
               if(facenodes[faceNum][eNum].GetSize() == 0 ){
                 facenodes[faceNum][eNum].Resize(numFaceNodes);
                 for ( UInt faceNode = 0;faceNode < numFaceNodes ;faceNode++ ) {
                   facenodes[faceNum][eNum][faceNode] = ++offset;
                   nodesType_[offset] = BaseFE::FACE;
                 }
               } else {
                 // additional check, this face got already mapped with different size
                 if( facenodes[faceNum][eNum].GetSize() != numFaceNodes ) {
                   WARN("Face #" << faceNum << " for element #" << eNum 
                        << " got " << facenodes[faceNum][eNum].GetSize()
                        << " faceNodes, whereas we want to set "
                        << numFaceNodes << " entries for element #" <<
                        elemNum );
                 }
               }
               //fill the virtual Nodes in the correct ordering

               for ( UInt i = 0; i < numFaceNodes ; i++ ) {
                 ftn.vNodes.Push_back(facenodes[faceNum][eNum][ permutations[i] ]);
               }
               ftn.offset.Push_back( permutations.GetSize() );
             }
           }
         }
         //===========================================================
         //Assign the Interior node numbers
         //===========================================================
         {
           //get the permutation Vector just for the number of nodes
           ptFe->GetNodalPermutation(permutations,actEl,BaseFE::INTERIOR,0);
           UInt numIntNodes = permutations.GetSize();
           EntityTypeNodes & itn =  virtualNodes_[actEl->elemNum][BaseFE::INTERIOR];

           //Check if the current element got already numbered
           if(interiornodes[actEl->elemNum].GetSize() == 0){
             interiornodes[actEl->elemNum].Resize(numIntNodes);
             for ( UInt intNode = 0;intNode < numIntNodes ;intNode++ ) {
               interiornodes[actEl->elemNum][intNode] = ++offset;
               nodesType_[offset] = BaseFE::INTERIOR;
             }
           }
           //fill the virtual Nodes in the correct ordering

           for ( UInt i = 0; i  < numIntNodes ; i++ ) {
             itn.vNodes.Push_back(interiornodes[actEl->elemNum][ permutations[i] ]);
           }
           itn.offset.Push_back( permutations.GetSize());
         } // loop elements
       }
     } // loop entity lists


     // trim all vectors to get unused memory back
     boost::unordered_map< UInt, ElemVirtualNodes >::iterator it = 
         virtualNodes_.begin();
     for( ; it != virtualNodes_.end(); ++it ) {
       ElemVirtualNodes & evn = it->second;
       ElemVirtualNodes::iterator evnIt = evn.begin();
       for( ; evnIt != evn.end(); evnIt++ ) {
         //        evnIt->second.vNodes.Trim();
         //        evnIt->second.offset.Trim();
       }

     }


   }

  void FeSpaceL2::PrintEqnMap(){
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    // obtain (fctId,eqnNr) -> feFct mapping from OLAS
    StdVector<UInt> blockNums, indices;
    AlgebraicSys * algSys = feFct->GetSystem();
    FeFctIdType fctId = feFct->GetFctId();
    algSys->MapCompleteFctIdToIndex(fctId, blockNums, indices);


    // =================================
    // 1) ELEMENT INFORMATION
    // =================================

    // iterate over all elements
    Grid* ptGrid = feFct->GetGrid();
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
        StdVector<UInt> edgeNodes;
        Integer edgeNum = ptElem->edges[i];
        ptElem->GetEdgeNodes( std::abs(edgeNum) , edgeNodes );
        std::cout << "E #" << edgeNum << " (" 
            << edgeNodes[0] << "-> " << edgeNodes[1] << "), ";
      }
      std::cout << "\n";

      // Print face  information
      std::cout << "Faces: ";
      for( UInt i=0, numFaces = ptElem->faces.GetSize(); i < numFaces; ++i ) {
        StdVector<UInt> faceNodes;
        UInt faceNum = ptElem->faces[i];
        ptElem->GetFaceNodes( faceNum, faceNodes );
        std::cout << "F #" << faceNum << " (" << faceNodes.ToString( 0 ) << "), ";
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
    shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();
    // ---------------
    //  NODES
    // ---------------
    boost::unordered_map< Integer , StdVector<Integer> >::iterator nodeIt = nodeMap_.eqns.begin();
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
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
        assert(feFct);
    if(interiorElemMap_.find(region) == interiorElemMap_.end()){
      //create a name for this list...
      std::string name1 = feFct->GetGrid()->GetRegion().ToString(region);
      std::string name2 = name1;
      name1 += "_dgInterior";
      name2 += "_dgInteriorOpposite";
      interiorElemMap_[region] = shared_ptr<NcSurfElemList>(new NcSurfElemList(feFct->GetGrid(),name1));
      interiorElemMapOpposite_[region] = shared_ptr<NcSurfElemList>(new NcSurfElemList(feFct->GetGrid(),name2));
    }
    surfElems = interiorElemMap_[region];
    opposingElems = interiorElemMapOpposite_[region];
  }

  void FeSpaceL2::GetInteriorSurfaceElems(RegionIdType region,
                                       shared_ptr<EntityList> & surfElems){
    // what we do right now is just to create empty shared pointers
    // and fill them later in finalize
    //create a name for this list...
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
        assert(feFct);
    std::string name = feFct->GetGrid()->GetRegion().ToString(region);
    name += "_dgInterior";
    interiorElemMap_[region] = shared_ptr<NcSurfElemList>(new NcSurfElemList(feFct->GetGrid(),name));
    surfElems = interiorElemMap_[region];
  }

  //! \copydoc FeSpace::GetExteriorSurfaceElemsOfFeSpace
  void FeSpaceL2::GetExteriorSurfaceElems(RegionIdType region, shared_ptr<EntityList> & surfElems){
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
        assert(feFct);
    if(exteriorElements_.find(region) == exteriorElements_.end()){
      std::string name = feFct->GetGrid()->GetRegion().ToString(region);
      name += "_dgExterior";
      exteriorElements_[region] = shared_ptr<NcSurfElemList>(new NcSurfElemList(feFct->GetGrid(),name));
    }
    surfElems = exteriorElements_[region];
  }

  void FeSpaceL2::CreateSurfaceElems(){
    //TODO: ADD DEBUGGING OUTPUT!!!!!
    //ok lets go. we collect the information from our maps and pass them to the grid
    std::set<RegionIdType> regions;
    std::map<RegionIdType,shared_ptr<NcSurfElemList> >::iterator regIter;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    //its enough to go through the interiorElemMap_
    regIter = interiorElemMap_.begin();
    for(;regIter != interiorElemMap_.end();++regIter){
      regions.insert(regIter->first);
    }
    //create vectors
    StdVector<shared_ptr<NcSurfElem> > interiorElemList;
    StdVector<shared_ptr<NcSurfElem> > exteriorElemList;
    feFct->GetGrid()->GenerateDGSurfaceElemes(regions,interiorElemList,exteriorElemList);
    //ok first we fill the standard map and, if requested, the opposing map
    for(UInt aElem =0;aElem < interiorElemList.GetSize();++aElem){
      RegionIdType curReg = interiorElemList[aElem]->regionId;
      shared_ptr<NcSurfElemList> curIntList = interiorElemMap_[curReg];
      curIntList->AddElement(interiorElemList[aElem]);
      //check if we need to fill in the opposite list
      if(interiorElemMapOpposite_.find(curReg) != interiorElemMapOpposite_.end()){
        shared_ptr<NcSurfElemList> curOppList = interiorElemMapOpposite_[curReg];
        //we have at least one neighbor...
        curOppList->AddElement(interiorElemList[aElem]->neighbors[0]);
        UInt numNeighbors = interiorElemList[aElem]->neighbors.GetSize();
        for(UInt i = 1; i<numNeighbors;++i){
          //we need to push back in both lists....
          curOppList->AddElement(interiorElemList[aElem]->neighbors[i]);
          curIntList->AddElement(interiorElemList[aElem]);
        }
      }
    }
    //ok now we should have valid interior element lists lets fill the exterior
    for(UInt aElem =0;aElem < exteriorElemList.GetSize();++aElem){
      //obtain region id
      RegionIdType curReg = exteriorElemList[aElem]->regionId;
      shared_ptr<NcSurfElemList> & curIntList = exteriorElements_[curReg];
      curIntList->AddElement(exteriorElemList[aElem]);
    }
    regIter = exteriorElements_.begin();
    for(;regIter != exteriorElements_.end();++regIter){
      feFct->AddEntityList(regIter->second);
    }

    regIter = interiorElemMap_.begin();
    for(;regIter != interiorElemMap_.end();++regIter){
      feFct->AddEntityList(regIter->second);
    }

  }
}
