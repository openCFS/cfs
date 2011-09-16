// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include "Elements/fespaceHCurlHi.hh"
#include "Elements/HCurlElemsHi.hh"
#include "DataInOut/Logging/cfslog.hh"

/*
 The FeSpace always knows just vertex/nodal, edge, face
 and inner degrees of freedoms. Depending on the grid and
 element dimension, the mapping is accordingly.
 This is why we should introduce some kind of mapping
   
   (GridDim,ElemDim)->UnknownType
 
 The question to be solved is "when" we have to know this
 mapping. Technically it is okay, that we do not have inner degrees
 at all e.g. in 2D where we have geometrically just vertices
 and faces. We only need the information in case we eliminate inner
 degrees of freedoms. 
 Maybe the best solution is to introduce a mapping "innerDofs" in the
 form
 
          ElemDim
 GridDim   1D  2D  3D
 ---------------------
   1D      E   -   -
   2D      -   F   -
   3D      -   -   I
   
 So technically, inner degrees just occur if the element
 AND the grid dimension are the same!! For all other cases, 
 we use the normal geomtric entity. This could reduce the mapping to
 just a vector / map:
 
 Dim  InnerUnknowns
 ------------------
 1D        E
 2D        F
 3D        I
 
 
 
 1D-Grid
 =========
 Eqns     Element-Dimension
          1D  2D  3D
 -----------------------------
 Vertex    V   -   -
 Edge      -   -   -
 Face      -   -   -
 Inner     E   -   -    I

 2D-Grid
 =========
 Eqns     Element-Dimension
          1D  2D  3D
 -----------------------------
 Vertex    V   V   -
 Edge      E   E   -
 Face      -   -   -
 Inner     -   F   -    I

 3D-Grid
 =========
 Eqns     Element-Dimension
           1D  2D  3D
 -----------------------------
 Vertex    V   V   V
 Edge      E   E   E
 Face      -   F   F
 Inner     -   -   I    I
 
 
*/

// declare class specific logging stream
DECLARE_LOG(feSpaceHCurlHi)
DEFINE_LOG(feSpaceHCurlHi, "feSpaceHCurlHi")
namespace CoupledField{

  //! Constructor
  FeSpaceHCurlHi::FeSpaceHCurlHi( PtrParamNode aNode, 
                                  PtrParamNode infoNode)
  : FeSpaceH1(aNode, infoNode) {
    mapType_ = POLYNOMIAL;
    type_ = HCURL;
    isHierarchical_ = true;
    polyType_ = LEGENDRE;
    
    infoNode_ = infoNode->Get("hCurlHierarchical");
  }

  //! Destructor
  FeSpaceHCurlHi::~FeSpaceHCurlHi(){
  }

  void FeSpaceHCurlHi::Init() {
    // read order of function space
    // read, if map type should be isotropic

    /*ParamNode * orderNode = NULL;
    orderNode = myParam_->Get("order");
    if( orderNode ) {
     if( orderNode->Has("grid") ) {
       isoOrder_ = 0; // has no real meaning here
       SetMapType( GRID );
     }
     if( orderNode->Has("uniform")) {
       isoOrder_ = orderNode->Get("uniform")->AsUInt();
       SetMapType(POLYNOMIAL);
     }
    }*/
    //read in polyLists and integLists for easier access later
    ReadIntegList();
    ReadPolyList();
  }
  
  void FeSpaceHCurlHi::SetStrategy( SolStrategyType strategy, UInt step ) {

    solStrategy_ = strategy;
    solStep_ = step;

    // Check: If we have TWO_LEVEL strategy and we are in the 2nd step we 
    // have to do the following:
    if(solStrategy_ == STRAT_TWO_LEVEL && solStep_ == 1 ) {
      
      // in the 1st step we use hard-codedo order 0
      //SetMapType(POLYNOMIAL);
      //isoOrder_ = 0;
      mapType_ = POLYNOMIAL;
    }
    
    if( solStrategy_ == STRAT_TWO_LEVEL && solStep_ == 2 ) {
      isFinalized_ = false;
      numEqns_ = 0;
      numFreeEquations_ = 0;
      numUnknowns_ = 0;
      bcCounter_.clear();
      gridToVirtualNodes_.clear();
      nodes_.Clear();
      nodesType_.clear();
      virtualNodes_.clear();
      nodeMap_.eqns.clear();
      nodeMap_.BcKeys.clear();
      nodeMap_.constraintNodes.clear();
    
    // Re-initialize structure
      this->Init();
      this->Finalize();
    }
  }


  void FeSpaceHCurlHi::SetRegionElements(RegionIdType region, 
                                         MappingType mType,
                                         const Matrix<Integer>& order,
                                         PtrParamNode infoNode ){
    
    //This method may not be called after the space is finalized!
    if(isFinalized_){
      Exception("FeSpace::SetRegionMapping is called after finalization");
    }

    //ToDo: save the information...
    // QUERY FOR USER PARAMS IS STILL TO COME
    refElems_[region][Elem::ET_QUAD4]  = new FeHCurlHiQuad();
    refElems_[region][Elem::ET_HEXA8]  = new FeHCurlHiHex();

    //now set the order
    if(order.GetNumCols() != 1 || order.GetNumRows() != 1){
      Exception("FeSpaceHCurlHi::SetRegionElements : Only Iso-Order is supported right now");
    }
    std::map<Elem::FEType, FeHCurlHi* >::iterator i = refElems_[region].begin();
    for( ; i != refElems_[region].end(); ++i ) {
      i->second->SetIsoOrder(order[0][0]+orderOffset_);
    }
    
    infoNode->Get("order")->SetValue(order[0][0]);

  }


  void FeSpaceHCurlHi::CheckConsistency(){

  }

  //! sets the default integration scheme and order
  void FeSpaceHCurlHi::SetDefaultIntegration(PtrParamNode infoNode ){
    regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].order = Matrix<Integer>(1,1);
    regionIntegration_[ALL_REGIONS].order[0][0] = 0;
    regionIntegration_[ALL_REGIONS].mode = RELATIVE;
  }

  UInt FeSpaceHCurlHi::GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                      UInt entityNum, UInt comp ) {
    if( type == BaseFE::NODE) {
      return 1;
    } else if( type == BaseFE::EDGE ) {
      // Currently we support just isotropic order
      //StdVector<UInt> & orders = edgeOrder_[elemNum];
      EXCEPTION("THis needs to be implemented");
      return 0;
    } else if( type == BaseFE::FACE ) {
      // IMPLEMENT ME
    }
    return 0;
  }

  UInt FeSpaceHCurlHi::GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                         UInt entityNum ) {
    if( type == BaseFE::NODE) {
      return 1;
    } else if( type == BaseFE::EDGE ) {
      // Currently we support just isotropic order
      //StdVector<UInt> & orders = edgeOrder_[elemNum];
      //return *(std::max_element( orders.Begin(), orders.End() ));
      EXCEPTION("THis needs to be implemented");
      return 0;
    } else if( type == BaseFE::FACE ) {
      // IMPLEMENT ME
    }
    return 0;
  }

  //! Map equations i.e. intialize object
  void FeSpaceHCurlHi::Finalize(){
    /* Basic idea:
     * 1. Create the VirtualNode Array
     * 2. Make the polynomial order of edges / faces consistent
     * 3. Map boundary conditions
     * 4. Map equations only based on the virtualNodeArray
     */

    // For the HCurl elements, we always need edges / faces
    feFunction_->GetGrid()->MapEdges();
    feFunction_->GetGrid()->MapFaces();

    CreateVirtualNodes();

    // Apply minimum rule to edges / faces to adjust polynomial order
    AdjustEntityOrder();
    
    //Determine boundary Unknowns
    MapNodalBCs();
    MapNodalEqns(1);
    MapNodalEqns(2);

    // Just for debugging purpose
    //PrintEqnMap();

    
    // Print information to .info.xml
    
    isFinalized_ = true;
  }
  BaseFE* FeSpaceHCurlHi::GetFe( const EntityIterator ent, 
                                 shared_ptr<IntScheme>& intScheme ) {
    BaseFE * ret = GetFe(ent);

    // Set correct integration order
    RegionIdType eRegion;// =  ent.GetElem()->regionId;
    if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
      eRegion = ent.GetSurfElem()->ptVolElem1->regionId;
    } else {
      eRegion = ent.GetElem()->regionId;
    }
       
    intScheme = intScheme_;
    IntScheme::IntegMethod  method;
    Matrix<Integer> order;
    this->GetIntegration(ret, eRegion, method, order);
    // Note: The order is currently more or less hard-coded for isotropic order
    intScheme->SetOrder( method, order[0][0] );

    return ret;

  }
  
  BaseFE* FeSpaceHCurlHi::GetFe( const EntityIterator ent ){
    RegionIdType eRegion;// =  ent.GetElem()->regionId;
    if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
        eRegion = ent.GetSurfElem()->ptVolElem1->regionId;
      } else {
        eRegion = ent.GetElem()->regionId;
      }
    
    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ent.GetElem()->type) == refElems_[eRegion].end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }

    BaseFE * myFe = refElems_[eRegion][ent.GetElem()->type];

    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entity orders.

    return myFe;
  }
  
  BaseFE* FeSpaceHCurlHi::GetFe( UInt elemNum ){
    const Elem * ptElem = feFunction_->GetGrid()->GetElem(elemNum); 
    RegionIdType eRegion = ptElem->regionId;
    

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ptElem->type) == refElems_[eRegion].end()){
      EXCEPTION("FeSpaceHCurlHi::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    BaseFE * myFe = refElems_[eRegion][ptElem->type];

    //NOTE: THE ORDER IS ALREADY SET IN THIS CASE
    //      BUT for more felexibility we should set it here
    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entitiy orders.
    //
    //myFe->SetIsoOrder( (UInt)regionMappings_[eRegion].second[0][0]);

    return myFe;
  }

  void FeSpaceHCurlHi::AdjustEntityOrder() {

    // Has to be re-implemented
//    UInt numComp = feFunction_->GetResultInfo()->dofNames.GetSize();
//
//    // loop over all elements
//    std::map< UInt, StdVector<UInt> >::iterator it = virtualEdges_.begin();
//    for( ; it != virtualEdges_.end(); it++ ) {
//
//      UInt elemNum = it->first;
//
//      // loop over all edges
//      StdVector<UInt > & edges = it->second;
//      UInt numEdges = it->second.GetSize();
//      for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ) {
//        UInt actEdge = edges[iEdge];
//
//        // check, if edge got already mapped
//        if( edgeOrder_.find(actEdge) == edgeOrder_.end() ) {
//          edgeOrder_[actEdge].Resize(numComp);
//          edgeOrder_[actEdge].Init(0);          
//        }
//
//        // Loop over all components
//        for( UInt iComp = 0; iComp < numComp; ++iComp ) {
//
//          // Check if the local order for the edge is larger than the 
//          // one already assigned
//          UInt locOrder = GetEntityOrder(elemNum, BaseFE::EDGE, iEdge, iComp+1 );
//          if(  locOrder > edgeOrder_[actEdge][iComp] ) {
//            edgeOrder_[actEdge][iComp] = locOrder;
//          } // if
//        } // loop over components
//      } // loop over edges
//    } // loop over elements
  }
  
  void FeSpaceHCurlHi::CreateVirtualNodes() {

    //follow the following algorithm
        
    // 1st loop over every element
    //  - get edges 
    //  - assign virtual nodes to lowest order unknowns (i.e. we only label
    //    the Nedelec unknowns)
    // 2nd loop over every element (can be omitted, if gradients get skipped)
    //  - get edges
    //  - assign virtual nodes to higher order edge unknowns (gradients)
    // 3rd loop over every element
    //  - label all face unknowns
    // 4th loop over every element
    //  - label all interior unknowns
    
    StdVector< shared_ptr<EntityList> > fctEntList;
    std::map<UInt, StdVector<Integer> > vertexNodes;
    std::map<UInt, StdVector<Integer> > edgenodes;
    std::map<UInt, StdVector<Integer> > facenodes;
    std::map<UInt, StdVector<Integer> > interiornodes;

    
    fctEntList = feFunction_->GetEntityList();
    StdVector<UInt> permutations; 
    
    //get the highest possible node number
    //UInt offset = feFunction_->GetGrid()->GetNumNodes();
    UInt offset =0;
    //UInt curRegIsoOrder = 0;
    MappingType curMap = GRID;
    // ------------------------------------------
    //  1st loop: lowest order Nedelec unknowns
    // ------------------------------------------
    // loop over all entitylists (i.e. regions)
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){
      LOG_DBG(feSpaceHCurlHi) << "treating entityList '" 
          << fctEntList[actList]->GetName() << "'";
      EntityIterator entIt = fctEntList[actList]->GetIterator();

      // loop over all elements
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        //check if we got what we expected
        if ( ! (entIt.GetType() == EntityList::ELEM_LIST ||
            entIt.GetType() == EntityList::SURF_ELEM_LIST) )  {
          break;
          EXCEPTION("FeSpaceHCurlHi::CreateVirtualNodes(): Calling the method with an unsupported EntityListType");
        }
        //cast down to element list
        //ElemList* actElemList = dynamic_cast<ElemList*>(fctEntList[actList].get());
        //RegionIdType curReg = actElemList->GetRegion();

        curMap = POLYNOMIAL;

        const Elem* actEl = entIt.GetElem();
        BaseFE* ptFe = GetFe( entIt ); 
        //distinguish between Grid or polynomial based mapping
        if( curMap == GRID ) {
          EXCEPTION("This makes no sense");
        } else if (curMap == POLYNOMIAL) {


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
                nodesType_[offset] = BaseFE::EDGE;
              }
            }

            //fill the virtual Nodes in the correct ordering
            for ( UInt i = 0; i < numEdgeNodes ; i++ ) {
              virtualNodes_[actEl->elemNum][BaseFE::EDGE].Push_back(edgenodes[edgeNum][ permutations[i] ]);
            }
          } // loop over edges
        } // if POLYNOMIAL
      } // loop over elements
    } // loop over entitylists
    
    // remember node range for lowest order edge functions
    nedelecNodeRange_[0] = 0;
    nedelecNodeRange_[1] = offset-1;
    
    
    // ------------------------------------------
    //  2nd loop: higher order edge unknowns
    // ------------------------------------------
    // .. we skip this here, as we skip the higher order gradients of the
    // edges anyway
    
    // ------------------------------------------
    //  3rd loop: higher order face unkowns
    // ------------------------------------------
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){
      LOG_DBG(feSpaceHCurlHi) << "treating entityList '" 
          << fctEntList[actList]->GetName() << "'";
      EntityIterator entIt = fctEntList[actList]->GetIterator();



      // loop over all elements
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        //check if we got what we expected
        if ( ! (entIt.GetType() == EntityList::ELEM_LIST ||
            entIt.GetType() == EntityList::SURF_ELEM_LIST) )  {
          break;
          EXCEPTION("FeSpaceHCurlHi::CreateVirtualNodes(): Calling the method with an unsupported EntityListType");
        }
        const Elem* actEl = entIt.GetElem();
        BaseFE* ptFe = GetFe( entIt ); 

        //cast down to element list
        //ElemList* actElemList = dynamic_cast<ElemList*>(fctEntList[actList].get());
        //RegionIdType curReg = actElemList->GetRegion();

        curMap = POLYNOMIAL;

        //distinguish between Grid or polynomial based mapping
        if( curMap == GRID ) {
          EXCEPTION("This makes no sense");
        } else if (curMap == POLYNOMIAL) {


          //===========================================================
          //Assign the Face node numbers
          //===========================================================
          UInt numFaceNodes = 0; 
          ElemShape actShape = Elem::shapes[actEl->type];
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
                nodesType_[offset] = BaseFE::FACE;
              }
            }
            //fill the virtual Nodes in the correct ordering
            for ( UInt i = 0; i < numFaceNodes ; i++ ) {
              virtualNodes_[actEl->elemNum][BaseFE::FACE].Push_back(facenodes[faceNum][ permutations[i] ]);
            }
          }
        } // if POLYNOMIAL
      } // loop over elements
    } // loop over entitylists
    
    // remember node range for lowest order edge functions
    faceNodeRange_[0] = nedelecNodeRange_[1]+1;
    faceNodeRange_[1] = offset-1;
    
    // ------------------------------------------
    //  4th loop: higher order interior unknowns
    // ------------------------------------------
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){
      LOG_DBG(feSpaceHCurlHi) << "treating entityList '" << fctEntList[actList]->GetName() << "'";
      EntityIterator entIt = fctEntList[actList]->GetIterator();

      // loop over all elements
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        //check if we got what we expected
        if ( ! (entIt.GetType() == EntityList::ELEM_LIST ||
            entIt.GetType() == EntityList::SURF_ELEM_LIST) )  {
          break;
          EXCEPTION("FeSpaceHCurlHi::CreateVirtualNodes(): Calling the method with an unsupported EntityListType");
        }
        const Elem* actEl = entIt.GetElem();
        BaseFE* ptFe = GetFe( entIt ); 
        //cast down to element list
        //ElemList* actElemList = dynamic_cast<ElemList*>(fctEntList[actList].get());
        //RegionIdType curReg = actElemList->GetRegion();

        curMap = POLYNOMIAL;
        //distinguish between Grid or polynomial based mapping
        if( curMap == GRID ) {
          EXCEPTION("This makes no sense");
        } else if (curMap == POLYNOMIAL) {

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
              nodesType_[offset] = BaseFE::INTERIOR;
            }
          }
          //fill the virtual Nodes in the correct ordering
          for ( UInt i = 0; i  < numIntNodes ; i++ ) {
            virtualNodes_[actEl->elemNum][BaseFE::INTERIOR].Push_back(interiornodes[actEl->elemNum][ permutations[i] ]);
          }
        } // if POLYNOMIAL
      } // loop over elements
    } // loop over entitylists

    LOG_TRACE(feSpaceHCurlHi) << "finished creation of virtual nodes";
   
    // remember node range for lowest order edge functions
    faceNodeRange_[0] = nedelecNodeRange_[1]+1;
    faceNodeRange_[1] = offset-1;
    
  }

  void FeSpaceHCurlHi::SetDefaultElements(PtrParamNode infoNode ){
    //but it could be, that the PDE requires a minimum order of elements...
    Matrix<Integer> order(1,1);
    if(orderOffset_>0){
      order[0][0] = orderOffset_;
    }else{
      order[0][0] = 1;
    }
    SetRegionElements(ALL_REGIONS,POLYNOMIAL,order, infoNode );
  }
} // end of namespace
