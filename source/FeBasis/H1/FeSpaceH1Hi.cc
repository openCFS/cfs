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
  FeSpaceH1Hi::FeSpaceH1Hi(PtrParamNode aNode, PtrParamNode infoNode )
  : FeSpaceH1(aNode, infoNode ) {

    type_ = H1;
    isHierarchical_ = true;
    polyType_ = LEGENDRE;
    
    infoNode_ = infoNode->Get("h1Hierarchical");
  }

  //! Destructor
  FeSpaceH1Hi::~FeSpaceH1Hi(){
  }
  
  //! Initialize class
  void FeSpaceH1Hi::Init() {
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

  UInt FeSpaceH1Hi::GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
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

  UInt FeSpaceH1Hi::GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
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
  void FeSpaceH1Hi::Finalize(){
    /* Basic idea:
     * 1. Create the VirtualNode Array
     * 2. Make the polynomial order of edges / faces consistent
     * 3. Map boundary conditions
     * 4. Map equations only based on the virtualNodeArray
     */

//    UInt numEqns_ = 0;
//    UInt numFreeEquations_ = 0;
    CreateVirtualNodes();

    // Apply minimum rule to edges / faces to adjust polynomial order
    AdjustEntityOrder();
        
    //Determine boundary Unknowns
    MapNodalBCs();
    MapNodalEqns(1);
    MapNodalEqns(2);

    isFinalized_ = true;
  }
  
  BaseFE* FeSpaceH1Hi::GetFe( const EntityIterator ent, 
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
    Matrix<Integer>  order;
    this->GetIntegration(ret, eRegion, method, order);
    // Note: The order is currently more or less hard-coded for isotropic order
    intScheme->SetOrder( method, order[0][0] );

    return ret;
  }

  BaseFE* FeSpaceH1Hi::GetFe( const EntityIterator ent ){

    // Note: if the element is a surface element, we must omit the regionId
    // and look for the neighbor. Which one to take? Well, we had the 
    // discussion already ....
    RegionIdType eRegion = NO_REGION_ID;
    if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
      eRegion = ent.GetSurfElem()->ptVolElem1->regionId;
    } else {
      eRegion = ent.GetElem()->regionId;
    }
    LOG_DBG3(feSpaceH1Hi) << "Returning FE #" << ent.GetElem()->elemNum
        << " of region " 
        << domain->GetGrid()->GetRegion().ToString(eRegion);
    
    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      LOG_DBG3(feSpaceH1Hi) << "\t-> No reference element found, use default region";
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ent.GetElem()->type) == refElems_[eRegion].end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }

    BaseFE * myFe = refElems_[eRegion][ent.GetElem()->type];

    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entitiy orders.
//    myFe->SetIsoOrder( isoOrder_);

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
    BaseFE * myFe = refElems_[eRegion][ptElem->type];

    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entitiy orders.
//    std::cerr << "isoOrder is " << isoOrder_ << std::endl;
////    myFe->SetIsoOrder( isoOrder_);

    return myFe;
  }

  void FeSpaceH1Hi::AdjustEntityOrder() {

    
    // Has to be reimplemented
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

  void FeSpaceH1Hi::SetRegionElements(RegionIdType region, 
                                      MappingType mType,
                                      const Matrix<Integer>& order,
                                      PtrParamNode infoNode ){
    
    LOG_DBG3(feSpaceH1Hi) << "FeSpaceH1HI: SetRegionElements for Region " <<
        domain->GetGrid()->GetRegion().ToString(region) << std::endl;
    
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
    
    // Generate reference elements for second order geoemtric element types
    refElems_[region][Elem::ET_LINE3]  = new FeH1HiLine();
    refElems_[region][Elem::ET_QUAD8]  = new FeH1HiQuad();
    refElems_[region][Elem::ET_TRIA6]  = new FeH1HiTria();
    refElems_[region][Elem::ET_HEXA20]  = new FeH1HiHex();
    refElems_[region][Elem::ET_WEDGE15] = new FeH1HiWedge();

    //now set the order
    if(order.GetNumCols() != 1 || order.GetNumRows() != 1){
      Exception("FeSpaceHCurlHi::SetRegionElements : Only Iso-Order is supported right now");
    }
    std::map<Elem::FEType, FeH1Hi* >::iterator i = refElems_[region].begin();
    for( ; i != refElems_[region].end(); ++i ) {
      i->second->SetIsoOrder(order[0][0]+orderOffset_);
    }
    
    // Important: indicate, that we do not re-use the originial grid based
    // nodes
    mapType_ = POLYNOMIAL;
    
    infoNode->Get("order")->SetValue(order[0][0]);
  }

  void FeSpaceH1Hi::CheckConsistency(){

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
    regionIntegration_[ALL_REGIONS].order = Matrix<Integer>(1,1);
    regionIntegration_[ALL_REGIONS].order[0][0] = 0;
    regionIntegration_[ALL_REGIONS].mode = RELATIVE;
  }
  
  // Conduct test for continuuity of basis functions
  void FeSpaceH1Hi::TestContinuity(){
    
    // Specify entity type (V/E/F)
    BaseFE::EntityType et = BaseFE::VERTEX;
    
    // Loop over all entities of that type
    UInt numEntities = 0;
    Grid * grid = domain->GetGrid();
    switch(et) {
      case BaseFE::VERTEX:
        numEntities = grid->GetNumNodes();
        break;
      case BaseFE::EDGE:
        numEntities = grid->GetNumEdges();
        break;
      case BaseFE::FACE:
        numEntities = grid->GetNumFaces();
        break;
      default:
        EXCEPTION("Not implemented");
    }
    
    UInt numElems = grid->GetNumElems();
    // Loop over all entities
    for( UInt iEnt = 0; iEnt < numEntities; ++iEnt ) {

      // Loop over all elements

      for(UInt iElem = 0; iElem < numElems; ++iElem ) {

        // if element contains that entity:
//        switch(et) {
//          
//        }

        //    -> Loop over all points of that element (distinguish shape type)

        //    -> Calculate basis functions

        //    -> extract the one of the given entity

        //    -> print result into file
      } // loop over elements
      // 
    } // loop over entities
    
  }
} // end of namespace
