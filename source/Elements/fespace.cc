// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "fespace.hh"
#include "DataInOut/Logging/cfslog.hh"

// include headers of sub-classes for factory method
#include "fespaceH1Lagrange.hh"
#include "fespaceH1Hi.hh"
#include "fespaceHCurlHi.hh"

namespace CoupledField {


  // declare class specific logging stream
  DECLARE_LOG(feSpace)
  DEFINE_LOG(feSpace, "feSpace");
  
  
  //! Constructor
  FeSpace::FeSpace(ParamNode * paramNode){
    
    myParam_ = paramNode;
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

  }

  FeSpace::~FeSpace(){
  }

  shared_ptr<FeSpace> FeSpace::CreateInstance(ParamNode* aNode  ) {
    
    // we obtain the <FeFunction> element
    std::string spaceStr, polyStr, orderStr;
    aNode->Get("spaceType",spaceStr);
    aNode->Get("polyType",polyStr);
    aNode->Get("order",orderStr);
    SpaceType spaceType = SpaceTypeEnum.Parse(spaceStr);
    PolyType polyType = PolyTypeEnum.Parse(polyStr);
    
    // switch depending on space type
    shared_ptr<FeSpace> ret;
    switch( spaceType ) {
      case H1:
        if( polyType == LAGRANGE )  {
          ret.reset(new FeSpaceH1Lagrange(aNode));
        } else {  
          ret.reset(new FeSpaceH1Hi(aNode));
        }
        break;
      case HCURL:
        if( polyType == LAGRANGE )  {
          // create explicit lower order HCurl space
        } else {
          ret.reset(new FeSpaceHCurlHi(aNode));
        }
        break;
      case CONST:
      case HDIV:
      case L2:
        EXCEPTION("Function space of type '"
            << spaceStr << "' not yet implemented!");
        break;
      case UNDEF_SPACE:
        EXCEPTION("Can not create unknown function space");
        break;
    }
    return ret;    
  }
  
  void FeSpace::SetIsoOrder( UInt order ) {
    // ToDo: Add some consistency checks
    isoOrder_ = order;
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

      EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element #"
          << ptElem->elemNum << " of region " 
          <<      domain->GetGrid()->RegionIdToName(ptElem->regionId));
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
    //  - get vertices and if not already mapped assign new virtual nodes
    //  - get edges and if not already mapped assign new virtual nodes
    //  - get faces and if not already mapped assign new virtual nodes
    //  - assign interior node to the element
    //  - now fill the Virtual node map with the information according to element orientation
    // - finally delete all intermediate arrays
    
    StdVector< shared_ptr<EntityList> > fctEntList;
    std::map<UInt, StdVector<Integer> > vertexNodes;
    std::map<UInt, StdVector<Integer> > edgenodes;
    std::map<UInt, StdVector<Integer> > facenodes;
    std::map<UInt, StdVector<Integer> > interiornodes;

    LOG_TRACE(feSpace) << "starting to create virtual nodes";
    
    fctEntList = feFunction_->GetEntityList();

    //get the highest possible node number
    //UInt offset = feFunction_->GetGrid()->GetNumNodes();
    UInt offset =0;

    // loop over all entitylists (i.e. regions)
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){

      LOG_DBG(feSpace) << "treating entityList '" << fctEntList[actList]->GetName() << "'";
      // Get iterator of current element list
      EntityIterator entIt = fctEntList[actList]->GetIterator();


      //a little helper to create the nodes_ array
      // ToDo: Check, if this is really necessary
      // Note: This is not generally true: If we have HCurl elements, we do
      // not have nodal unknowns at all, i.e. the mapping type "GRID" doest
      // not really mean anything except the order of the elements to be used
      // should be the same as the geometric order of the grid.
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
        //UInt nodalUnknowns = ptFe->GetNumFncs( BaseFE::VERTEX );

        // Continue, if no nodal unknowns are present
        //if( nodalUnknowns == 0 ) continue;

        //distinguish between Grid or polinomial based mapping
        if( mapType_ == GRID ) {
          const StdVector<UInt> & elemNodes = actEl->connect;
          virtualNodes_[actEl->elemNum][BaseFE::VERTEX] = elemNodes;
        } else if (mapType_ == POLYNOMIAL) {
          //===========================================================
          //Assign the BaseFE::VERTEX node numbers
          //===========================================================
          LOG_DBG2(feSpace) << "mapping vertex nodes";
          StdVector<UInt> permutations; // initially size 0
          UInt numVertexNodes = 0;
          
          UInt numVert = Elem::shapes[actEl->type].numVertices;
          StdVector<UInt> elemNodes = actEl->connect;
          for ( UInt iVert= 0; iVert< numVert; iVert++ ) {
            UInt vertexNum = elemNodes[iVert];
            ptFe->GetNodalPermutation(permutations,actEl,BaseFE::VERTEX,iVert);            
            numVertexNodes = permutations.GetSize();

            // Check if the vertex is already numbered.
            if( vertexNodes[vertexNum].GetSize() == 0 && isContinuous_ ) {
              vertexNodes[vertexNum].Resize(numVertexNodes);

              // Note: if we have just one "vertexnode" (i.e. every vertex is
              // just related to one unknown, we use to global node number
//              if( numVertexNodes == 1) {
//                vertexNodes[vertexNum][0] = vertexNum;
//                nodes_.Push_back(vertexNum);
//                // in any case, we have to increase the offset
//                offset = feFunction_->GetGrid()->GetNumNodes();
//              } else {
                for( UInt vertNode = 0; vertNode < numVertexNodes; ++vertNode ) {
                  vertexNodes[vertexNum][vertNode] = ++offset;
                  LOG_DBG3(feSpace) << "adding " << offset << " to node_";
                  nodes_.Push_back(offset);
//                }
              }
            } // is continuous
            
            for( UInt i = 0; i < numVertexNodes; ++i ) {
              virtualNodes_[actEl->elemNum][BaseFE::VERTEX].Push_back(vertexNodes[vertexNum][permutations[i] ]);
              LOG_DBG3(feSpace) << "adding " << vertexNodes[vertexNum][permutations[i] ]
                                << " as virtual vertex node to element " << actEl->elemNum;
            }
            if(gridToVirtualNodes_.find(vertexNum) == gridToVirtualNodes_.end()){
              LOG_DBG3(feSpace) << "gridToVirtualNodes[" << vertexNum << "] = " << offset;
              gridToVirtualNodes_[vertexNum] =  offset;
            }
          } // loop over vertices
          //Create the permutation array
          
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
              interiornodes[actEl->elemNum][intNode] = ++offset;
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
  
  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************

  // Definition of finite element space types
   static EnumTuple spaceTypeTuples[] = {
     EnumTuple(FeSpace::UNDEF_SPACE, "Undef"), 
     EnumTuple(FeSpace::CONST,       "Const"), 
     EnumTuple(FeSpace::H1,          "H1"),
     EnumTuple(FeSpace::HCURL,       "HCurl"),
     EnumTuple(FeSpace::HDIV,        "HDiv"),
     EnumTuple(FeSpace::L2,          "L2"),
   };
   Enum<FeSpace::SpaceType> FeSpace::SpaceTypeEnum = \
      Enum<FeSpace::SpaceType>("Types of FE Spaces",
          sizeof(spaceTypeTuples) / sizeof(EnumTuple),
          spaceTypeTuples);
   
   // Definition of polynomial types
    static EnumTuple polyTypeTuples[] = {
      EnumTuple(FeSpace::UNDEF_POLY, "Undef"), 
      EnumTuple(FeSpace::LAGRANGE,   "Lagrange"), 
      EnumTuple(FeSpace::LEGENDRE,   "Legendre"),
      EnumTuple(FeSpace::JACOBI,     "Jacobi"),
    };
    Enum<FeSpace::PolyType> FeSpace::PolyTypeEnum = \
       Enum<FeSpace::PolyType>("Types of FE Spaces",
           sizeof(polyTypeTuples) / sizeof(EnumTuple),
           polyTypeTuples);
   
  
  // Definition of types of boundary conditions
  static EnumTuple bcTypeTuples[] = {
    EnumTuple(FeSpace::NOBC,    "NOBC"), 
    EnumTuple(FeSpace::HDBC,    "HDBC"),
    EnumTuple(FeSpace::IDBC,    "IDBC"),
    EnumTuple(FeSpace::CS,    "CS")
  };
  
  Enum<FeSpace::BcType> FeSpace::BcTypeEnum = \
     Enum<FeSpace::BcType>("Types of boundary conditions",
         sizeof(bcTypeTuples) / sizeof(EnumTuple),
         bcTypeTuples);
} // end of namespace 
