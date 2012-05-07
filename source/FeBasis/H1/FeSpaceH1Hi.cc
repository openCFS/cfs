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
  : FeSpaceHi(aNode, infoNode, ptGrid ) {

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
    std::map< RegionIdType, std::map<Elem::FEType, FeH1Hi* > >::iterator regionIt;
    regionIt = refElems_.begin();
    for( ; regionIt != refElems_.end(); ++regionIt ) {
      std::map<Elem::FEType, FeH1Hi* > & elems = regionIt->second;
      std::map<Elem::FEType, FeH1Hi* >::iterator elemIt = elems.begin();
      for( ; elemIt != elems.end(); ++elemIt ) {
        delete elemIt->second;
      }
    }
  }
  
  //! Initialize class
  void FeSpaceH1Hi::Init( shared_ptr<SolStrategy> solStrat ) {
    
    solStrat_ = solStrat;
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
    
    // In case we have vectorial unknowns and the polynomial degree is
    // anisotropic, we have to fix certain dofs, as the shape functions
    // comming from the element are always scalar.
    FixHigherOrderAnisoDofs();
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
    std::map<RegionIdType,ApproxOrder>::iterator it = regionOrder_.find(eRegion);
    assert( it != regionOrder_.end() );
    SetElemOrder( ent.GetElem(), myFe, it->second, true );

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
    std::map<RegionIdType,ApproxOrder>::iterator it = regionOrder_.find(eRegion);
    SetElemOrder( ptElem, myFe, it->second, true );
    
    return myFe;
  }


  void FeSpaceH1Hi::SetRegionElements(RegionIdType region, 
                                      MappingType mType,
                                      const ApproxOrder& order,
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
    
    SetRegionOrder( region, order );
    
    infoNode->Get("order")->SetValue(order.ToString());
    
   
  }

  void FeSpaceH1Hi::CheckConsistency(){
    // nothing to do here

  }

  void FeSpaceH1Hi::SetDefaultElements( PtrParamNode infoNode ){
    //but it could be, that the PDE requires a minimum order of elements...
    ApproxOrder order(ptGrid_->GetDim());
    order.SetIsoOrder(1);
    if(orderOffset_>0){
      order.SetIsoOrder(orderOffset_);
    }
    SetRegionElements( ALL_REGIONS, POLYNOMIAL, order, infoNode );
  }

  //! sets the default integration scheme and order
  void FeSpaceH1Hi::SetDefaultIntegration( PtrParamNode infoNode ){
    regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].order.SetIsoOrder( 0 );
    regionIntegration_[ALL_REGIONS].mode = RELATIVE;
  }

  FeHi* FeSpaceH1Hi::GetFeHi( RegionIdType region, Elem::FEType type ) {
    FeHi * ret = NULL;
    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(region) == refElems_.end()){
      region = ALL_REGIONS;
    }

    if(refElems_[region].find(type) == refElems_[region].end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    ret = refElems_[region][type]; 
    return ret;
  }



} // end of namespace
