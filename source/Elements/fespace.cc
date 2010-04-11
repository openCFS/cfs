// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "fespace.hh"
#include "DataInOut/Logging/cfslog.hh"


namespace CoupledField {


  // declare class specific logging stream
  DECLARE_LOG(feSpace)
  DEFINE_LOG(feSpace, "feSpace");
  
  
  //! Constructor
  FeSpace::FeSpace(){
    
    isFinalized_ = false;
    isContinuous_ = true;
    isoOrder_ = 0;
    numEqns_ = 0;
    numUnknowns_ = 0;
    numFreeEquations_ = 0;
    mapType_ = GRID;

    bcCounter_[NOBC] = 0;
    bcCounter_[HDBC] = 0;
    bcCounter_[IDBC] = 0;
    bcCounter_[CS] = 0;

    //now call the function defining its integration scheme

  }

  FeSpace::~FeSpace(){
  }

  void FeSpace::SetIsoOrder( UInt order ) {
    // ToDo: Add some consistency checks
    isoOrder_ = order;
  }

  void FeSpace::PreCalcShapeFncs(){
     //now precalculate all available integration points
     //stupid but simple
     //get grip of the integrationScheme object
     shared_ptr<IntegrationScheme> integScheme = feFunction_->GetGrid()->GetIntegrationScheme();

     StdVector<LocPoint> integPoints;
     std::map<Elem::FEType, BaseFE* >::iterator elemIt = refElems_.begin();
     while(elemIt != refElems_.end()){
       integScheme->GetAllIntegrationPoints(integPoints,elemIt->first);
       elemIt->second->SetFunctionsAtIp(integPoints);
       elemIt++;
     }
   }

  void FeSpace::GetNodesOfEntities( StdVector<UInt>& nodes,
                                    shared_ptr<EntityList> ent,
                                    BaseFE::EntityType entType ) {

    // Case 1: No "virtual" nodes are present. Therefore we can directly use
    //          the grid-nodes of the geometric element
    if (mapType_ == GRID){ // get name of entitylist
      std::string name= ent->GetName();
      feFunction_->GetGrid()->GetNodesByName( nodes, name );
      
      // Case 2: We rely on the the information of the virtualNodes_ array.
      //         which gets filled depending on the number of "unknown" nodes 
      //          initially.
    } else if (mapType_ == POLYNOMIAL) {
      EntityList::ListType defType = ent->GetType();
      EntityIterator entIt = ent->GetIterator();
      nodes.Resize(0);
      UInt nCount = 0; 
      switch (defType){
        case EntityList::REGION_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<Elem*> eList;
            feFunction_->GetGrid()->GetElems(eList,entIt.GetRegion());
            StdVector<UInt> eNodes;
            for (UInt i = 0; i < eList.GetSize(); i ++ ){

              GetNodesOfElement(eNodes,eList[i],entType);
              for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
                nodes.Push_back(eNodes[aNode]);
              }
            }
          }
        break;
        case EntityList::ELEM_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            GetNodesOfElement(eNodes,entIt.GetElem(),entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::SURF_ELEM_LIST:
          Warning(" FeSpaceH1::GetNodesOfEntities(): Going to treat a SURF_ELEM_LIST as a ELEM_LIST...");
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            GetNodesOfElement(eNodes,entIt.GetElem(),entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::NODE_LIST:
          nCount = 0;
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            if(gridToVirtualNodes_.find(entIt.GetNode()) != gridToVirtualNodes_.end()){
              nodes.Push_back( gridToVirtualNodes_[entIt.GetNode()]);
            }
          }
        break;
        default:
          EXCEPTION("FeSpace::GetNodesOfEntities: Called with a list which is not supported ");
          break;
      }
    }
  }

  void FeSpace::GetNodesOfElement( StdVector<UInt>& nodes,
                          const Elem* ptElem,
                          BaseFE::EntityType entType){
    UInt elemNum = ptElem->elemNum;
    if(virtualNodes_.find(elemNum) ==virtualNodes_.end()){
      EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element Number");
    }
    if(entType == BaseFE::ALL){
      nodes.Resize(virtualNodes_[elemNum][BaseFE::VERTEX].GetSize()+
                   virtualNodes_[elemNum][BaseFE::EDGE].GetSize()+
                   virtualNodes_[elemNum][BaseFE::FACE].GetSize()+
                   virtualNodes_[elemNum][BaseFE::INTERIOR].GetSize());
      ElemVirtualNodes::iterator nodeIt = virtualNodes_[elemNum].begin();
      UInt c = 0;
      while(nodeIt !=virtualNodes_[elemNum].end()){

        StdVector<UInt> entNodes =  nodeIt->second;
        for (UInt i = 0; i < entNodes.GetSize(); i++ ){
          nodes[c++] =  entNodes[i];
        }
        nodeIt++;
      }
    }else{
      nodes.Resize(virtualNodes_[elemNum][entType].GetSize());
      StdVector<UInt> entNodes =  virtualNodes_[elemNum][entType];
      for (UInt i = 0; i < entNodes.GetSize(); i++ ){
        nodes[i] =  entNodes[i];
      }
    }
  }
  
  void FeSpace::CreateVirtualNodes(){
    //follow the following algorithm
    // - loop over every element
    //  - get edges and if not already mapped assign new virtual nodes
    //  - get faces and if not already mapped assign new virtual nodes
    //  - assign interior node to the element
    //  - now fill the Virtual node map with the information according to element orientation
    // - finally delete all intermediate arrays
    
    StdVector< shared_ptr<EntityList> > fctEntList;
    std::map<UInt, StdVector<Integer> > edgenodes;
    std::map<UInt, StdVector<Integer> > facenodes;
    std::map<UInt, StdVector<Integer> > interiornodes;

    LOG_TRACE(feSpace) << "starting to create virtual nodes";
    
    fctEntList = feFunction_->GetEntityList();

    //get the highest possible node number
    UInt offset = feFunction_->GetGrid()->GetNumNodes();

    // loop over all entitylists (i.e. regions)
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){

      LOG_DBG(feSpace) << "treating entityList '" << fctEntList[actList]->GetName() << "'";
      // Get iterator of current element list
      EntityIterator entIt = fctEntList[actList]->GetIterator();


      //a little helper to create the nodes_ array
      // ToDo: Check, if this is really necessary
      if( mapType_ == GRID ) {
        StdVector<UInt> curNodes;
        GetNodesOfEntities( curNodes, fctEntList[actList] );
        for ( UInt aNode= 0; aNode < curNodes.GetSize(); aNode++ ) {
          if(gridToVirtualNodes_.find(curNodes[aNode]) == gridToVirtualNodes_.end()){
            gridToVirtualNodes_[curNodes[aNode]] =  curNodes[aNode];
          }
        }
      }

      
      // loop over all elements
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        //check if we got what we expected
        if ( ! (entIt.GetType() == EntityList::ELEM_LIST ||
            entIt.GetType() == EntityList::SURF_ELEM_LIST) )  {
          break;
          // This is not in general a problem. We can have mixed boundary conditions
          //EXCEPTION("FeSpaceH1Lagrange::CreateVirtualNodes(): Calling the method with an unsupported EntityListType");
        }

        // Fetch current finite element. This is performed by the specialized 
        // version, so this element "knows" already about its order, unknowns etc.
        const Elem* actEl = entIt.GetElem();
        BaseFE* ptFe = GetFe( entIt ); 


        // Check, if the element has nodal unknowns at all
        UInt nodalUnknowns = ptFe->GetNumFncs( BaseFE::VERTEX );

        // Continue, if no nodal unknowns are present
        if( nodalUnknowns == 0 ) continue;

        //distinguish between Grid or polinomial based mapping
        if( mapType_ == GRID ) {
          const StdVector<UInt> & elemNodes = actEl->connect;
          virtualNodes_[actEl->elemNum][BaseFE::VERTEX] = elemNodes;
        } else if (mapType_ == POLYNOMIAL) {
          //===========================================================
          //Assign the BaseFE::VERTEX node numbers
          //===========================================================
          
          // In the case of disconinuous elements, I do not think that
          // the following code works.
          if ( !isContinuous_ ) {
            Warning("The nodal mapping is not yet adapted to discontinuous mapping.");
          }
          
          StdVector<UInt> elemNodes = actEl->connect;
          UInt numVert = Elem::shapes[actEl->type].numVertices;
          virtualNodes_[actEl->elemNum][BaseFE::VERTEX].Resize(numVert);
          for ( UInt aNode= 0; aNode < numVert; aNode++ ) {
            virtualNodes_[actEl->elemNum][BaseFE::VERTEX][aNode] = elemNodes[aNode];
            if(gridToVirtualNodes_.find(elemNodes[aNode]) == gridToVirtualNodes_.end()){
              gridToVirtualNodes_[elemNodes[aNode]] =  elemNodes[aNode];
              nodes_.Push_back(elemNodes[aNode]);
            }
          }
          //Create the permutation array
          StdVector<UInt> permutations; // initially size 0
          //if(isoOrder_ > 1){
          feFunction_->GetGrid()->MapEdges();
          feFunction_->GetGrid()->MapFaces();
          //}
          //===========================================================
          //Assign the Edge node numbers
          //===========================================================
          UInt numEdgeNodes = 0;
          ElemShape actShape = Elem::shapes[actEl->type];
          for ( UInt iEdge=0; iEdge < actShape.numEdges; iEdge++) {
            UInt edgeNum = std::abs(actEl->edges[iEdge]);
            //get the permutation Vector
            ptFe->GetNodalPermutation(permutations,actEl,BaseFE::EDGE,iEdge); 
            numEdgeNodes = permutations.GetSize(); 

            // Check if the edge is already numbered. 
            // Additionally, if we have the case of discontinuous approximation,
            // we number the nodes separately for every element separately anyway.
            if(edgenodes[edgeNum].GetSize() == 0 &&  isContinuous_) {
              //here we assume spectral element approximation and we have 
              //order-1 nodes on the edge
              edgenodes[edgeNum].Resize(numEdgeNodes);
              for ( UInt edgeNode = 0;edgeNode < numEdgeNodes ;edgeNode++ ) {
                edgenodes[edgeNum][edgeNode] = ++offset;
                nodes_.Push_back(offset);
              }
            }

            //fill the virtual Nodes in the correct ordering
            for ( UInt i = 0; i < numEdgeNodes ; i++ ) {
              virtualNodes_[actEl->elemNum][BaseFE::EDGE].Push_back(edgenodes[edgeNum][ permutations[i] ]);
            }
          }

          //===========================================================
          //Assign the Face node numbers
          //===========================================================
          UInt numFaceNodes = 0; 
          for ( UInt iFace=0; iFace < actShape.numFaces; iFace++) {
            UInt faceNum = actEl->faces[iFace];
            //get the permutation Vector
            ptFe->GetNodalPermutation(permutations,actEl,BaseFE::FACE,iFace); 
            numFaceNodes = permutations.GetSize(); 

            // Check if the face is already numbered. 
            // Additionally, if we have the case of discontinuous approximation,
            // we number the nodes separately for every element separately anyway.
            if(facenodes[faceNum].GetSize() == 0 && isContinuous_ ){
              //here we assume spectral element approximation and we have 
              //order-1 nodes on the edge
              facenodes[faceNum].Resize(numFaceNodes);
              for ( UInt faceNode = 0;faceNode < numFaceNodes ;faceNode++ ) {
                facenodes[faceNum][faceNode] = ++offset;
                nodes_.Push_back(offset);
              }
            }
            //fill the virtual Nodes in the correct ordering
            for ( UInt i = 0; i < numFaceNodes ; i++ ) {
              virtualNodes_[actEl->elemNum][BaseFE::FACE].Push_back(facenodes[faceNum][ permutations[i] ]);
            }
          }

          //===========================================================
          //Assign the Interior node numbers
          //===========================================================
          //get the permutation Vector just for the number of nodes
          ptFe->GetNodalPermutation(permutations,actEl,BaseFE::INTERIOR,0); 
          UInt numIntNodes = permutations.GetSize(); 

          //Check if the current element got already numbered
          if(interiornodes[actEl->elemNum].GetSize() == 0){
            //here we assume spectral element approximation and we have 
            //order-1 nodes on the edge
            interiornodes[actEl->elemNum].Resize(numIntNodes);
            for ( UInt intNode = 0;intNode < numIntNodes ;intNode++ ) {
              facenodes[actEl->elemNum][intNode] = ++offset;
              nodes_.Push_back(offset);
            }
          }
          //fill the virtual Nodes in the correct ordering
          for ( UInt i = 0; i  < numIntNodes ; i++ ) {
            virtualNodes_[actEl->elemNum][BaseFE::INTERIOR].Push_back(interiornodes[actEl->elemNum][ permutations[i] ]);
          }
        } // if POLYNOMIAL

      } // loop elements 
    } // loop entity lists

    //another little helper to create the nodes_ array
    if( mapType_ == GRID ) {
      nodes_.Resize(gridToVirtualNodes_.size());
      UInt counter = 0;
      for(std::map<UInt,UInt>::const_iterator it = gridToVirtualNodes_.begin(); it != gridToVirtualNodes_.end(); ++it) {
        nodes_[counter++] = it->second;
      }
    }
    LOG_TRACE(feSpace) << "finished creation of virtual nodes";
  }
  
  void FeSpace::CreateVirtualEdges(){
    // Perform following algorithm:
    // - Map edges of grid
    // - Loop over all elements
    // - Check, if element has edge degrees
    // - If edge was not mapped or discontinuous => give local edge number
    
    LOG_TRACE(feSpace) << "starting to create virtual edges"; 
    
    // fetch regions
    StdVector< shared_ptr<EntityList> > fctEntList;
    fctEntList = feFunction_->GetEntityList();
    
    // Make sure, that edges get mapped
    feFunction_->GetGrid()->MapEdges();
    
    std::set<UInt> mappedEdges; // remember all mapped edges
    UInt actEdge = 0;
    
    // loop over all entitylists (i.e. regions)
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){

      LOG_DBG(feSpace) << "treating entityList '" << fctEntList[actList]->GetName() << "'";

      // loop over all elements
      EntityIterator entIt = fctEntList[actList]->GetIterator();
      for(entIt.Begin(); !entIt.IsEnd();entIt++){
        
        // Fetch current finite element. This is performed by the specialized 
        // version, so this element "knows" already about its order, unknowns etc.
        const Elem* actEl = entIt.GetElem();
        BaseFE* ptFe = GetFe( entIt );
        LOG_DBG3(feSpace) << "treating element #" << actEl->elemNum;
        
        // Check, if the element has edge unknowns at all
        UInt edgeUnknowns = ptFe->GetNumFncs( BaseFE::EDGE );
        if( edgeUnknowns == 0 ) continue;
        
        UInt numEdges = Elem::shapes[actEl->type].numEdges;
        
        // Safety check, if element got already mapped
        if( virtualEdges_[actEl->elemNum].GetSize() != 0) 
          continue;
        
        virtualEdges_[actEl->elemNum].Resize(numEdges);

        // Loop over edges and perform mapping
        for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ) {
          UInt edgeNum = std::abs(actEl->edges[iEdge]);
          LOG_DBG3(feSpace) << "treating element edge" << edgeNum;
          // check, if global edge is already numbered or if we are discontinuous
          UInt virtEdgeNum = 0;
          std::map<UInt,UInt>::iterator it = gridToVirtualEdges_.find(edgeNum);
          if( it != gridToVirtualEdges_.end() && isContinuous_) {
            virtEdgeNum = it->second;
          } else {
            // in case of continuous mapping, we simply take the global edge number
            if( isContinuous_) {
              virtEdgeNum = edgeNum;
            } else {
              virtEdgeNum = ++actEdge;
            }
            gridToVirtualEdges_[edgeNum] = virtEdgeNum;
            LOG_DBG3(feSpace) << "Mapping global edge E" << edgeNum << " to " << actEdge;
          } // if
          virtualEdges_[actEl->elemNum][iEdge] = virtEdgeNum;
        } // loop over edges
      } // loop over elements
    } // loop over entitylists
    
    LOG_TRACE(feSpace) << "finished creation of virtual edges";
  } 

  void FeSpace::CreateVirtualFaces(){
    // Perform following algorithm:
    // - Map faces of grid
    // - Loop over all elements
    // - Check, if element hat face degrees of freedoms
    // - If face was not mapped or discontinuous => give local face number
    LOG_TRACE(feSpace) << "starting to create virtual faces"; 

    // fetch regions
    StdVector< shared_ptr<EntityList> > fctEntList;
    fctEntList = feFunction_->GetEntityList();

    // Make sure, that faces get mapped
    feFunction_->GetGrid()->MapFaces();

    std::set<UInt> mappedFaces; // remember all mapped faces
    UInt actFace = 0;

    // loop over all entitylists (i.e. regions)
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){

      LOG_DBG(feSpace) << "treating entityList '" << fctEntList[actList]->GetName() << "'";

      // loop over all elements
      EntityIterator entIt = fctEntList[actList]->GetIterator();
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        // Fetch current finite element. This is performed by the specialized 
        // version, so this element "knows" already about its order, unknowns etc.
        const Elem* actEl = entIt.GetElem();
        BaseFE* ptFe = GetFe( entIt );

        // Check, if the element has edge unknowns at all
        UInt faceUnknowns = ptFe->GetNumFncs( BaseFE::FACE );
        if( faceUnknowns == 0 ) continue;
        
        
        UInt numFaces = Elem::shapes[actEl->type].numFaces;

        // Safety check, if element got already mapped
        if( virtualFaces_[actEl->elemNum].GetSize() != 0) 
          continue;

        virtualFaces_[actEl->elemNum].Resize(numFaces);

        // Loop over faces and perform mapping
        for( UInt iFace = 0; iFace < numFaces; ++iFace ) {
          UInt faceNum = actEl->faces[iFace];

          // check, if global face is already numbered or if we are discontinuous
          UInt virtFaceNum = 0;
          std::map<UInt,UInt>::iterator it = gridToVirtualFaces_.find(faceNum);
          if( it != gridToVirtualFaces_.end() && isContinuous_) {
            virtFaceNum = it->second;
          } else {
            virtFaceNum = ++actFace;
            gridToVirtualFaces_[faceNum] = actFace;
          } // if
          virtualFaces_[actEl->elemNum][iFace] = virtFaceNum;
        } // loop over faces
      } // loop over elements
    } // loop over entitylists

    LOG_TRACE(feSpace) << "finished creation of virtual faces";
  }
  
  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************

  // Definition of finite element space types
  static EnumTuple spaceTypeTuples[] = {
    EnumTuple(FeSpace::CONST,    "CONST"), 
    EnumTuple(FeSpace::H1_LO,    "H1_LO"),
    EnumTuple(FeSpace::H1_HI,    "H1_HI"),
    EnumTuple(FeSpace::HCURL_LO, "HCURL_LO"),
    EnumTuple(FeSpace::HCURL_HI, "HCURL_HI"),
    EnumTuple(FeSpace::HDIV_LO,  "HDIV_LO"),
    EnumTuple(FeSpace::HDIV_HI,  "HDIV_HI"),
    EnumTuple(FeSpace::L2_LO,    "L2_LO"),
    EnumTuple(FeSpace::L2_HI,    "L2_HI")
    
  };
  Enum<FeSpace::SpaceType> FeSpace::spaceType = \
     Enum<FeSpace::SpaceType>("Types of FE Spaces",
         sizeof(spaceTypeTuples) / sizeof(EnumTuple),
         spaceTypeTuples);
  
  // Definition of types of boundary conditions
  static EnumTuple bcTypeTuples[] = {
    EnumTuple(FeSpace::NOBC,    "NOBC"), 
    EnumTuple(FeSpace::HDBC,    "HDBC"),
    EnumTuple(FeSpace::IDBC,    "IDBC"),
    EnumTuple(FeSpace::CS,    "CS")
  };
  
  Enum<FeSpace::BcType> FeSpace::bcType = \
     Enum<FeSpace::BcType>("Types of boundary conditions",
         sizeof(bcTypeTuples) / sizeof(EnumTuple),
         bcTypeTuples);
} // end of namespace 
