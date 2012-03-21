#include "FeSpaceH1Hi.hh"
#include "H1ElemsHi.hh"

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
#include "DataInOut/Logging/LogConfigurator.hh"

DECLARE_LOG(feSpaceH1Hi)
DEFINE_LOG(feSpaceH1Hi, "feSpaceH1Hi")
namespace CoupledField{

  //! Constructor
  FeSpaceH1Hi::FeSpaceH1Hi(PtrParamNode aNode, PtrParamNode infoNode,
                           Grid* ptGrid )
  : FeSpaceH1(aNode, infoNode, ptGrid ) {

    type_ = H1;
    isHierarchical_ = true;
    polyType_ = LEGENDRE;
    mapType_ = POLYNOMIAL; 
    
    infoNode_ = infoNode->Get("h1Hierarchical");
    
    // important: trigger mapping of edges / faces
    ptGrid_->MapEdges();
    ptGrid_->MapFaces();
  }

  //! Destructor
  FeSpaceH1Hi::~FeSpaceH1Hi(){
  }
  
  //! Initialize class
  void FeSpaceH1Hi::Init( shared_ptr<SolStrategy> solStrat ) {
    
    solStrat_ = solStrat;
    // read order of function space
    // read, if map type should be isotropic

    /*  ParamNode * orderNode = NULL;
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

  void FeSpaceH1Hi::Finalize(){
    /* Basic idea:
     * 1. Adjust entity order of all faces / edges
     * 2. Create the VirtualNode Array
     * 3. Make the polynomial order of edges / faces consistent
     * 4. Map boundary conditions
     * 5. Map equations only based on the virtualNodeArray
     */
    AdjustEntityOrder();
    CreateVirtualNodes();

    //Determine boundary Unknowns
    MapNodalBCs();
    MapNodalEqns(1);
    MapNodalEqns(2);

    isFinalized_ = true;
  }
  
  BaseFE* FeSpaceH1Hi::GetFe( const EntityIterator ent ,
                              IntScheme::IntegMethod& method,
                              IntegOrder & order  ) {
    BaseFE * ret = GetFe(ent);

    // Set correct integration order
    RegionIdType eRegion;// =  ent.GetElem()->regionId;
    if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
      eRegion = GetVolElem(ent.GetSurfElem()) ->regionId;
    } else {
      eRegion = ent.GetElem()->regionId;
    }

    this->GetIntegration(ret, eRegion, method, order);

    return ret;
  }

  BaseFE* FeSpaceH1Hi::GetFe( const EntityIterator ent ){

    // Note: if the element is a surface element, we must omit the regionId
    // and look for the neighbor. Which one to take? Well, we had the 
    // discussion already ....
    RegionIdType eRegion = NO_REGION_ID;
    if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
      eRegion = GetVolElem(ent.GetSurfElem()) ->regionId;
    } else {
      eRegion = ent.GetElem()->regionId;
    }
    LOG_DBG3(feSpaceH1Hi) << "Returning FE #" << ent.GetElem()->elemNum
        << " of region "  << ptGrid_->GetRegion().ToString(eRegion);
    
    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      LOG_DBG3(feSpaceH1Hi) << "\t-> No reference element found, use default region";
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ent.GetElem()->type) == refElems_[eRegion].end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }

    // Fetch reference element and set correct order
    FeH1Hi * myFe = refElems_[eRegion][ent.GetElem()->type];
    SetElemOrder( ent.GetElem(), myFe, regionOrder_[eRegion], true );

    LOG_DBG(feSpaceH1Hi) << "Returning FE #" << ent.GetElem()->elemNum
        << " with " << myFe->BaseFE::GetNumFncs() << " functions";
    
    return myFe;
  }
  
  BaseFE* FeSpaceH1Hi::GetFe( UInt elemNum ){
    const Elem * ptElem = feFunction_->GetGrid()->GetElem(elemNum); 
    RegionIdType eRegion = ptElem->regionId;

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ptElem->type) == refElems_[eRegion].end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    // Fetch reference element and set correct order
    FeH1Hi * myFe = refElems_[eRegion][ptElem->type];
    SetElemOrder( ptElem, myFe, regionOrder_[eRegion], true );
    
    return myFe;
  }


  void FeSpaceH1Hi::SetRegionElements(RegionIdType region, 
                                      MappingType mType,
                                      const Matrix<Integer>& order,
                                      PtrParamNode infoNode ){
    
    LOG_DBG3(feSpaceH1Hi) << "FeSpaceH1HI: SetRegionElements for Region " <<
        ptGrid_->GetRegion().ToString(region) << std::endl;
    
    //This method may not be called after the space is finalized!
    if(isFinalized_){
      Exception("FeSpaceH1Hi::SetRegionMapping is called after finalization");
    }

    //TODO: Save the information somehow
    // QUERY FOR USER PARAMS IS STILL TO COME
    // Generate reference elements for first order geoemtric element types
    refElems_[region][Elem::ET_LINE2]  = new FeH1HiLine();
    refElems_[region][Elem::ET_QUAD4]  = new FeH1HiQuad();
    refElems_[region][Elem::ET_TRIA3]  = new FeH1HiTria();
    refElems_[region][Elem::ET_HEXA8]  = new FeH1HiHex();
    refElems_[region][Elem::ET_WEDGE6] = new FeH1HiWedge();
    
    // Generate reference elements for second order geometric element types
    refElems_[region][Elem::ET_LINE3]  = new FeH1HiLine();
    refElems_[region][Elem::ET_QUAD8]  = new FeH1HiQuad();
    refElems_[region][Elem::ET_TRIA6]  = new FeH1HiTria();
    refElems_[region][Elem::ET_HEXA20]  = new FeH1HiHex();
    refElems_[region][Elem::ET_WEDGE15] = new FeH1HiWedge();


    // store order for this region
    regionOrder_[region] = order;
    
    
    // loop over all elements of this region
    LOG_DBG(feSpaceH1Hi) << "Assigning edge / face order";
    StdVector<Elem*> elems;
    ptGrid_->GetElems( elems, region );
    UInt numElems = elems.GetSize();
    for( UInt i = 0; i < numElems; ++i ) {
      
      const Elem& el = *(elems[i]);
      LOG_DBG3(feSpaceH1Hi) << "Treating element #" << el.elemNum;
      
      // Fetch reference element and set correct order
      FeH1Hi * myFe = refElems_[region][el.type];
      
      // important: do not adjust the entity order here, as the
      // edge / face order is not yet initialized
      SetElemOrder( elems[i], myFe, regionOrder_[region], false );
      
      // a) loop over all edges
      // -----------------------
      UInt numEdges = el.edges.GetSize();
      LOG_DBG3(feSpaceH1Hi) << "Checking " << numEdges << " edges ";
      const StdVector<UInt>& edgeOrders = myFe->GetEdgeOrder();
      for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ){
        UInt edgeNum = std::abs(el.edges[iEdge]);
        UInt edgeOrder = edgeOrders[iEdge];
        
        // check if edge got already mapped
        if( orderEdges_.find(edgeNum) == orderEdges_.end() ) {
          orderEdges_[edgeNum] = edgeOrder;
          LOG_DBG3(feSpaceH1Hi) << "\tedge #" << edgeNum << ", order: " 
                                << edgeOrder;
        } else {
          UInt oldOrder = orderEdges_[edgeNum];
          LOG_DBG3(feSpaceH1Hi) << "\tedge #" << edgeNum 
                                << " -> Already set to " << oldOrder; 
          // check if previous set order is different from current
          if( edgeOrder != oldOrder) {
            orderEdges_[edgeNum] = std::max( edgeOrder, oldOrder );
            LOG_DBG3(feSpaceH1Hi) << "\t\t=> Re-Set to " << orderEdges_[edgeNum];
            adjustedEdges_.insert( edgeNum );
          }
        }
      } // loop over edges
      
      
      // b) loop over all faces (only in 3D case )
      // -----------------------------------------
      if( ptGrid_->GetDim() == 3 ) {
        UInt numFaces = el.faces.GetSize();
        LOG_DBG3(feSpaceH1Hi) << "Checking " << numFaces << " faces ";
        const StdVector<boost::array<UInt,2> >& faceOrders = 
            myFe->GetFaceOrder();
        for( UInt iFace = 0; iFace < numFaces; ++iFace ){
          UInt faceNum = el.faces[iFace];
          boost::array<UInt,2> faceOrder = faceOrders[iFace];

          // check if face got already mapped
          if( orderFaces_.find(faceNum) == orderFaces_.end() ) {
            orderFaces_[faceNum] = faceOrder;
            LOG_DBG3(feSpaceH1Hi) << "\tface #" << faceNum << ", order: " 
                                  << faceOrder[0] << ", " << faceOrder[1];
          } else {
            boost::array<UInt,2> oldOrder = orderFaces_[faceNum];
            LOG_DBG3(feSpaceH1Hi) << "\tface #" << faceNum 
                << " -> Already set to " << oldOrder[0] << ", " 
                << oldOrder[1];
            // check if previous set order is different from current
            if( faceOrder != oldOrder) {
              orderFaces_[faceNum][0] = std::max( faceOrder[0], oldOrder[0] );
              orderFaces_[faceNum][1] = std::max( faceOrder[1], oldOrder[1] );
              LOG_DBG3(feSpaceH1Hi) << "\t\t=> Re-Set to " 
                                    << orderFaces_[faceNum][0] << ", "
                                    << orderFaces_[faceNum][1];
              adjustedFaces_.insert( faceNum );
            }
          }
        }  // loop over faces
      } 

    } // loop over elements
    infoNode->Get("order")->SetValue(order[0][0]);
  }

  void FeSpaceH1Hi::CheckConsistency(){
    // nothing to do here

  }

  void FeSpaceH1Hi::SetDefaultElements( PtrParamNode infoNode ){
    //but it could be, that the PDE requires a minimum order of elements...
    Matrix<Integer> order(1,1);
    if(orderOffset_>0){
      order[0][0] = orderOffset_;
    }else{
      order[0][0] = 1;
    }
    SetRegionElements( ALL_REGIONS, POLYNOMIAL, order, infoNode );
  }

  //! sets the default integration scheme and order
  void FeSpaceH1Hi::SetDefaultIntegration( PtrParamNode infoNode ){
    regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].order.SetIsoOrder( 0 );
    regionIntegration_[ALL_REGIONS].mode = RELATIVE;
  }
  
  
  void FeSpaceH1Hi::
  SetElemOrder( const Elem* ptEl, FeH1Hi* ptFe, 
                const Matrix<Integer>& order,
                bool applyMaxRule ) {
    LOG_DBG(feSpaceH1Hi) << "In SetElemOrder for elem " << ptEl->elemNum;
    // Stage 1: Set order as given from the region template.
    // This order template specified either an isotropic polynomial degree
    // or a direction dependent order (for each unknown).
    
    // check for isotropy -> This has to be changed to a more
    // general representation soon
    if( order.GetNumRows() == 1 ) {
      ptFe->SetIsoOrder( order[0][0] );
    } else {
      EXCEPTION( "Anisotropic order not implemented yet");
    }
    
    // Stage 2: Adjust edge / face order according to minimum /
    // maximum rule, in case the element has neighboring elements with
    // different order. We wil use the maximum rule
    if( applyMaxRule && orderEdges_.size() ) {
      // loop over all edges
      UInt numEdges = ptEl->edges.GetSize();
      for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ) {
        UInt edgeNum = std::abs( ptEl->edges[iEdge] );
          // check if edge got adjusted
        if( orderEdges_.find(edgeNum) != orderEdges_.end() ) {
          LOG_DBG3(feSpaceH1Hi) << "Setting edge " << edgeNum
              << " to order " << orderEdges_[edgeNum];
          ptFe->SetEdgeOrder( iEdge, orderEdges_[edgeNum] );
        }
      }
      
    }
    
    if( applyMaxRule && orderFaces_.size() 
        && ptGrid_->GetDim() == 3  ) {
      // loop over all faces
      UInt numFaces = ptEl->faces.GetSize();
      for( UInt iFace = 0; iFace < numFaces; ++iFace ) {
        UInt faceNum = ptEl->faces[iFace];
        // check if face got adjusted
        if( orderFaces_.find(faceNum) != orderFaces_.end() ) {
          LOG_DBG3(feSpaceH1Hi) << "Setting face " << faceNum
              << " to order " << orderFaces_[faceNum][0]
              << ", " << orderFaces_[faceNum][1];
          ptFe->SetFaceOrder( iFace, orderFaces_[faceNum] );
        }
      }

    }
  }
  
  void FeSpaceH1Hi::AdjustEntityOrder() {
    LOG_TRACE(feSpaceH1Hi) << "Adjusting entity order";
    
    // This method is called in the ned
    
    
    // Create temporary map for edges
    boost::unordered_map<UInt, UInt> edgeOrders;
    
    // Loop over all adjusted edges
    LOG_DBG(feSpaceH1Hi) << "Adjusting edge order";
    boost::unordered_set<UInt>::const_iterator edgeIt = adjustedEdges_.begin();
    for( ; edgeIt != adjustedEdges_.end(); ++edgeIt ) {
      LOG_DBG3(feSpaceH1Hi) << "\tOrder of edge #" << *edgeIt 
                            << " is " <<  orderEdges_[*edgeIt];
      edgeOrders[*edgeIt] = orderEdges_[*edgeIt];
    } // loop over edges

    // store only adjusted edges back
    orderEdges_ = edgeOrders;
    
    
    if( ptGrid_->GetDim() < 3) return;
    
    // Create temporary map for faces
    boost::unordered_map<UInt, boost::array<UInt,2> > faceOrders;
    
    // Loop over all adjusted faces
    LOG_DBG(feSpaceH1Hi) << "Adjusting face order";
    boost::unordered_set<UInt>::const_iterator faceIt = adjustedFaces_.begin();
    for( ; edgeIt != adjustedFaces_.end(); ++faceIt ) {
      LOG_DBG3(feSpaceH1Hi) << "\tOrder of face #" << *faceIt 
                            << " is " <<  orderFaces_[*faceIt][0]
                            << ", " << orderFaces_[*faceIt][1];
      faceOrders[*faceIt] = orderFaces_[*faceIt];
    } // loop over faces

    // store only adjusted faces back
    orderFaces_ = faceOrders;
    
    
    
    
  }
  

} // end of namespace
