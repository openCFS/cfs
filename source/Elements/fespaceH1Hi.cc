// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#include "Elements/fespaceH1Hi.hh"
#include "Elements/H1ElemsHi.hh"

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
namespace CoupledField{

  //! Constructor
  FeSpaceH1Hi::FeSpaceH1Hi(PtrParamNode aNode)
  : FeSpaceH1(aNode) {

    type_ = H1;
    isHierarchical_ = true;
    polyType_ = LEGENDRE;
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
    CreateRefElems();
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
//    CreateVirtualEdges();
//    CreateVirtualFaces();

    // Apply minimum rule to edges / faces to adjust polynomial order
    AdjustEntityOrder();
    
    
        
    //Determine boundary Unknowns
    MapNodalBCs();
    MapNodalEqns(1);
    MapNodalEqns(2);

    // Just for debugging purpose
    //PrintEqnMap();

    isFinalized_ = true;
  }

  BaseFE* FeSpaceH1Hi::GetFe( const EntityIterator ent ){
    RegionIdType eRegion =  ent.GetElem()->regionId;

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ent.GetElem()->type) == refElems_[eRegion].end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }

    BaseFE * myFe = refElems_[eRegion][ent.GetElem()->type];

    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entitiy orders.
    myFe->SetIsoOrder( isoOrder_);

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
    myFe->SetIsoOrder( isoOrder_);

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

  void FeSpaceH1Hi::SetRegionElements(RegionIdType region, MappingType mType,Matrix<Integer> order){
    //This method may not be called after the space is finalized!
    if(isFinalized_){
      Exception("FeSpaceH1Hi::SetRegionMapping is called after finalization");
    }

    //TODO: Save the information somehow
    // QUERY FOR USER PARAMS IS STILL TO COME
    refElems_[region][Elem::ET_LINE2]  = new FeH1HiLine();
    refElems_[region][Elem::ET_QUAD4]  = new FeH1HiQuad();
    refElems_[region][Elem::ET_HEXA8]  = new FeH1HiHex();

    //now set the order
    if(order.GetNumCols() != 1 || order.GetNumRows() != 1){
      Exception("FeSpaceHCurlHi::SetRegionElements : Only Iso-Order is supported right now");
    }
    std::map<Elem::FEType, BaseFE* >::iterator i = refElems_[region].begin();
    for( ; i != refElems_[region].end(); ++i ) {
      i->second->SetIsoOrder(order[0][0]+orderOffset_);
    }
  }


  void FeSpaceH1Hi::SetRegionIntegration(RegionIdType region, IntScheme::IntegMethod method, Matrix<Integer> order){
    //TODO:Implementation of defaults (ALL_REGIONS) and XML
    regionIntegration_[region].first = method;
    regionIntegration_[region].second = order;
  }

  void FeSpaceH1Hi::CheckConsistency(){

  }
  void FeSpaceH1Hi::ProcessPolyRegionNode(PtrParamNode node, RegionIdType region){
    Matrix<Integer> order(1,1);
    order[0][0] = -1;
    PtrParamNode isoOrderNode = node->Get("isoOrder", ParamNode::PASS );
    PtrParamNode anIsoOrderNode = node->Get("anIsoOrder", ParamNode::PASS );

    if(isoOrderNode){
      Integer isoOrder = isoOrderNode->As<Integer>();
      order[0][0] = isoOrder;
    }else if(anIsoOrderNode){
      //TO BE DONE
      EXCEPTION("Anisotropic element orders are not supported");
    }else{
      WARN("Did not find a order node. setting it to 1");
      order[0][0] = 1;
    }
    SetRegionElements(region,POLYNOMIAL,order);
  }

  void FeSpaceH1Hi::CreateDefaultElements(){
    //but it could be, that the PDE requires a minimum order of elements...
    Matrix<Integer> order(1,1);
    if(orderOffset_>0){
      order[0][0] = orderOffset_;
    }else{
      order[0][0] = 1;
    }
    SetRegionElements(ALL_REGIONS,POLYNOMIAL,order);
  }

  //! sets the default integration scheme and order
  void FeSpaceH1Hi::SetDefaultIntegration(){
    regionIntegration_[ALL_REGIONS].first = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].second = Matrix<Integer>(1,1);
    regionIntegration_[ALL_REGIONS].second[0][0] = 2;
  }
} // end of namespace
