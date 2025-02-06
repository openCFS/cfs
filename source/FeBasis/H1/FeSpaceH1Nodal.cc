// =====================================================================================
// 
//       Filename:  fespaceH1Nodal.cc
// 
//    Description:  This implementes the space of lagrange elements the most
//                  important method here is the finalize method which fills the
//                  virtual node map.
// 
//        Version:  0.1
//        Created:  01/21/2010 10:42:37 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "FeSpaceH1Nodal.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

DEFINE_LOG(feSpaceH1Nodal, "feSpaceH1Nodal")
namespace CoupledField{

  //! Constructor
  FeSpaceH1Nodal::FeSpaceH1Nodal(PtrParamNode aNode,
                                 PtrParamNode infoNode,
                                 Grid* ptGrid )
  : FeSpaceNodal(aNode, infoNode, ptGrid ) {
    type_ = H1;
    isHierarchical_ = false;
    mapType_ = GRID;
    polyType_ = LAGRANGE;
    order_ = 0;

    infoNode_ = infoNode->Get("h1Nodal");
  }

  //! Destructor
  FeSpaceH1Nodal::~FeSpaceH1Nodal(){
    std::map< RegionIdType, std::map<Elem::FEType, FeH1* > >::iterator regionIt;
    regionIt = refElems_.begin();
    for( ; regionIt != refElems_.end(); ++regionIt ) {
      std::map<Elem::FEType, FeH1* > & elems = regionIt->second;
      std::map<Elem::FEType, FeH1* >::iterator elemIt = elems.begin();
      for( ; elemIt != elems.end(); ++elemIt ) {
        delete elemIt->second;
      }
    }
#ifdef USE_OPENMP
    std::map< RegionIdType, TLMap<Elem::FEType, FeH1* > >::iterator regIt = TL_RefElems_.begin();
    while(regIt != TL_RefElems_.end()){
      regIt->second.Clear();
      ++regIt;
    }
#endif
  }

  void FeSpaceH1Nodal::Init( shared_ptr<SolStrategy> solStrat ) {
    
    solStrat_ = solStrat;
    
    //read in polyLists and integLists for easier access later
    ReadIntegList();
    ReadPolyList();
  }

  BaseFE* FeSpaceH1Nodal::GetFe( const EntityIterator ent ,
                                 IntScheme::IntegMethod& method,
                                 IntegOrder & order  ) {
    BaseFE * ret = GetFe(ent);

    // Set correct integration order
    RegionIdType eRegion = GetVolElem(ent.GetElem())->regionId;

    this->GetIntegration(ret, eRegion, method, order);
    return ret;

  }

  BaseFE* FeSpaceH1Nodal::GetFe( const EntityIterator ent ){

    if(ent.GetType() != EntityList::ELEM_LIST &&
        ent.GetType() != EntityList::SURF_ELEM_LIST &&
        ent.GetType() != EntityList::NC_ELEM_LIST){
      EXCEPTION("This version of GetFe expects a element iterator")
    }


    // Note: if the element is a surface element, we must omit the regionId
    // and look for the neigbor. Which one to take? Well, we had the
    // discussion already ....
    RegionIdType eRegion = GetVolElem(ent.GetElem())->regionId;

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    Elem::FEType eType = ent.GetElem()->type;

    if(refElems_[eRegion].find(eType) == refElems_[eRegion].end()){
#ifdef _WIN32
      EXCEPTION(__FUNCTION__
                << ": requested fetype which is not supported by space");
#else
      EXCEPTION(__PRETTY_FUNCTION__
                << ": requested fetype which is not supported by space");
#endif
    }
#ifdef USE_OPENMP
    BaseFE * myFe;
    if(isFinalized_ && omp_get_num_threads()>1)
      myFe = TL_RefElems_[eRegion][ent.GetElem()->type];
    else
      myFe = refElems_[eRegion][ent.GetElem()->type];
#else
    BaseFE * myFe = refElems_[eRegion][ent.GetElem()->type];
#endif

    /* No need to set the order here, as this is already done once and for all
     * in the SetMapType() method. For higher order spaces with non-uniform
     * polynomial order, this necessary.
     */
    // myFe->SetIsoOrder( isoOrder_);

    return myFe;
  }

  BaseFE* FeSpaceH1Nodal::GetFe( UInt elemNum ){
    // request a strong pointer
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock();
    assert(feFct);
    const Elem * ptElem = feFct->GetGrid()->GetElem(elemNum);
    RegionIdType eRegion = GetVolElem(ptElem)->regionId;

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ptElem->type) == refElems_[eRegion].end()){
#ifdef _WIN32
      EXCEPTION(__FUNCTION__
                << ": requested FEType ("<<ptElem->type<<") is not supported by space");
#else
      EXCEPTION(__PRETTY_FUNCTION__
                << ": requested FEType ("<<ptElem->type<<") is not supported by space");
#endif
    }
#ifdef USE_OPENMP
    BaseFE * myFe;
    if(isFinalized_ && omp_get_num_threads()>1)
      myFe = TL_RefElems_[eRegion][ptElem->type];
    else
      myFe = refElems_[eRegion][ptElem->type];
#else
    BaseFE * myFe = refElems_[eRegion][ptElem->type];
#endif
    return myFe;

  }



  void FeSpaceH1Nodal::PreCalcShapeFncs(){
    //now pre-calculate shape functions for the
    // desired element orders
    //get grip of the integrationScheme object

    // leave, if element space is hierarchical
    if (isHierarchical_)return;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    shared_ptr<IntScheme> integScheme = feFct->GetGrid()->GetIntegrationScheme();

    std::map<Integer,LocPoint> integPoints;
    std::map<RegionIdType, std::map<Elem::FEType, FeH1* > >::iterator regIt = refElems_.begin();

    Elem::ShapeType shape;
    while(regIt != refElems_.end()){
      std::map<Elem::FEType, FeH1* >::iterator elemIt = regIt->second.begin();
      while(elemIt != regIt->second.end()){
        shape = Elem::GetShapeType(elemIt->first);
        //now obtain iterator to the integ Regions associated with this polyRegion
        std::set<RegionIdType>::iterator integIter = polyToIntegMap[regIt->first].begin();
        while(integIter != polyToIntegMap[regIt->first].end()){
          IntegOrder curOrder;
          IntScheme::IntegMethod curMethod;
          GetIntegration(elemIt->second,*integIter,curMethod,curOrder);
          integScheme->GetIntegrationPoints(integPoints,shape,curMethod,curOrder);
          dynamic_cast<FeNodal*>(elemIt->second)->SetFunctionsAtIp(integPoints);
          integIter++;
        }
        elemIt++;
      }
      regIt++;
    }
  }


  //! Map equations i.e. intialize object
  void FeSpaceH1Nodal::Finalize(){
    /* Basic idea:
     * 1. Set correct element order at elements
     * 2. Create the VirtualNode Array
     * 3. Map boundary conditions
     * 4. Map equations only based on the virtualNodeArray
     */

    //      UInt numEqns_ = 0;
    //      UInt numFreeEquations_ = 0;
    CreateVirtualNodes();
    //Determine boundary Unknowns
    MapNodalBCs();
    MapNodalEqns(1);
    MapNodalEqns(2);
    
    // TEMPORARY: print information
    PrintEqnMap();
    CheckConsistency();
#ifdef USE_OPENMP
    std::map< RegionIdType, std::map<Elem::FEType, FeH1* > >::iterator regIt = refElems_.begin();
    while(regIt != refElems_.end()){
      TL_RefElems_[regIt->first] = regIt->second;
      ++regIt;
    }
#endif
    isFinalized_ = true;
  }

  void FeSpaceH1Nodal::SetRegionElements(RegionIdType region,
                                         MappingType mType,
                                         const ApproxOrder& order,
                                         PtrParamNode infoNode ){

    LOG_DBG(feSpaceH1Nodal) << "Setting region elements";
    LOG_DBG2(feSpaceH1Nodal) << "\tregion: "
        << ptGrid_->GetRegion().ToString(region);
    LOG_DBG2(feSpaceH1Nodal) << "\tmappingType: " << mType;
    LOG_DBG2(feSpaceH1Nodal) << "\torder: " << order.ToString();

    //This method may not be called after the space is finalized!
    if(isFinalized_){
      Exception("FeSpace::SetRegionMapping is called after finalization");
    }
    
    // Ensure, that we have an isotropic approximation order
    if( !order.IsIsotropic() ) {
      EXCEPTION( "Lagrangian polynomials can be only used with "
          << "isotropic approximation order");
    }
    
    if(mType == GRID){
      
      UInt gridOrder = ptGrid_->IsQuadratic() ? 2 : 1;
      infoNode->Get("order")->SetValue(gridOrder);
      order_ = gridOrder;
      // Check if a polynomial order was already set and if they match
      if( order_ != 0 && order_ != gridOrder ) {
        EXCEPTION( "A Lagrangian-based space can only have one specific "
            << "polynomial order" );
      }
      
      if( refElems_[region].empty() ) 
      {
        refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLine1();
        refElems_[region][Elem::ET_TRIA3]  = new FeH1LagrangeTria1();
        refElems_[region][Elem::ET_TRIA6]  = new FeH1LagrangeTria2();
        refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuad1();
        refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHex1();
        refElems_[region][Elem::ET_WEDGE6] = new FeH1LagrangeWedge1();
        refElems_[region][Elem::ET_PYRA5]  = new FeH1LagrangePyra1();
        refElems_[region][Elem::ET_TET4]  = new FeH1LagrangeTet1();
        refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLine2();
        refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
        refElems_[region][Elem::ET_QUAD9]  = new FeH1LagrangeQuad9();
        refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex2();
        refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHex27();
        refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
        refElems_[region][Elem::ET_WEDGE18] = new FeH1LagrangeWedge18();
        refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra2();
        refElems_[region][Elem::ET_PYRA14] = new FeH1LagrangePyra14();
        refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet2();
      }
      else 
      {
        //WARN( "Reference elements for region " << region << " already defined!" );
      }

    } else if (mType == POLYNOMIAL) {
      
      UInt isoOrder = order.GetIsoOrder();
      
      // Check if a polynomial order was already set and if they match
      //if( order_ != 0 && order_ != isoOrder ) {
      //  EXCEPTION( "A Lagrangian-based space can only have one specific "
      //      << "polynomial order" );
      //}

      mapType_ = POLYNOMIAL;
      infoNode->Get("order")->SetValue(isoOrder);
      if(order.GetPolyMapType() == ApproxOrder::SERENDIPITY_TYPE) {
        infoNode->Get("polyMapType")->SetValue("serendipity");
      } else {
        infoNode->Get("polyMapType")->SetValue("tensor");
      }
      
      
      order_ = isoOrder;

      // In case we have polynomial order 1, we can simply use the explicit
      // type elements. In all other cases, we have to check, if TENSOR_TYPE
      // or SERENDIPITY_TYPE elements are desired and create the related ones.
      if( isoOrder == 1 ) {
        refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLine1();
        refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLine1();
        
        refElems_[region][Elem::ET_TRIA3]  = new FeH1LagrangeTria1();
        refElems_[region][Elem::ET_TRIA6]  = new FeH1LagrangeTria1();
        
        refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuad1();
        refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuad1();
        refElems_[region][Elem::ET_QUAD9]  = new FeH1LagrangeQuad1();

        refElems_[region][Elem::ET_TET4]   = new FeH1LagrangeTet1();
        refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet1();
        
        refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHex1();
        refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex1();
        refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHex1();
        
        refElems_[region][Elem::ET_WEDGE6]  = new FeH1LagrangeWedge1();
        refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge1();
        refElems_[region][Elem::ET_WEDGE18] = new FeH1LagrangeWedge1();
        
        refElems_[region][Elem::ET_PYRA5]  = new FeH1LagrangePyra1();
        refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra1();
        refElems_[region][Elem::ET_PYRA14] = new FeH1LagrangePyra1();
      } 
      
      if(isoOrder == 2) {
        if(order.GetPolyMapType() == ApproxOrder::SERENDIPITY_TYPE) {
        // In this case we can utilize the explicit second order elements 
        // of serendipity type
          refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLine2();
          refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLine2();

          refElems_[region][Elem::ET_TRIA3]  = new FeH1LagrangeTria2();
          refElems_[region][Elem::ET_TRIA6]  = new FeH1LagrangeTria2();

          refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuad2();
          refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
          refElems_[region][Elem::ET_QUAD9]  = new FeH1LagrangeQuad2();

          refElems_[region][Elem::ET_TET4]   = new FeH1LagrangeTet2();
          refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet2();

          refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHex2();
          refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex2();
          refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHex2();

          refElems_[region][Elem::ET_WEDGE6]  = new FeH1LagrangeWedge2();
          refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
          refElems_[region][Elem::ET_WEDGE18] = new FeH1LagrangeWedge2();

          refElems_[region][Elem::ET_PYRA5]  = new FeH1LagrangePyra2();
          refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra2();
          refElems_[region][Elem::ET_PYRA14] = new FeH1LagrangePyra2();
          
        } else if(order.GetPolyMapType() == ApproxOrder::TENSOR_TYPE) {
          // In this case we can utilize the explicit second order elements 
          // of serendipity type
          refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
          refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLineVar();

          refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
          refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
          refElems_[region][Elem::ET_QUAD9]  = new FeH1LagrangeQuadVar();

          refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
          refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHexVar();
          refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHexVar();

          std::map<Elem::FEType, FeH1* >::iterator i = refElems_[region].begin();
          for( ; i != refElems_[region].end(); ++i ) {
            FeH1LagrangeVar * ptFe = dynamic_cast<FeH1LagrangeVar*>(i->second);
            ptFe->SetIsoOrder(isoOrder+orderOffset_);
          }
          // triangular and tetrahedral elements have no mid-face nodes,
          // so here the serendipity type and tensor-type elements
          // are the same
          refElems_[region][Elem::ET_TRIA3]  = new FeH1LagrangeTria2();
          refElems_[region][Elem::ET_TRIA6]  = new FeH1LagrangeTria2();

          refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet2();
          refElems_[region][Elem::ET_TET4]   = new FeH1LagrangeTet2();

          // Use explicit tensor type hexas
//          refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHex27();
//          refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex27();
//          refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHex27();
          
          // Use explicit tensor type wedges
          refElems_[region][Elem::ET_WEDGE6]  = new FeH1LagrangeWedge18();
          refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge18();
          refElems_[region][Elem::ET_WEDGE18] = new FeH1LagrangeWedge18();

          // Use explicit tensor type pyramids
          refElems_[region][Elem::ET_PYRA5]  = new FeH1LagrangePyra14();
          refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra14();
          refElems_[region][Elem::ET_PYRA14] = new FeH1LagrangePyra14();
        }
      }

      if( isoOrder > 2 && order.GetPolyMapType() == ApproxOrder::TENSOR_TYPE) {
        // Here we have now arbitrary order elements wit tensor-type
        // polynomial order.
        refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
        refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLineVar();

        refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
        refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
        refElems_[region][Elem::ET_QUAD9]  = new FeH1LagrangeQuadVar();

        refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
        refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHexVar();
        refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHexVar();

        std::map<Elem::FEType, FeH1* >::iterator i = refElems_[region].begin();
        for( ; i != refElems_[region].end(); ++i ) {
          FeH1LagrangeVar * ptFe = dynamic_cast<FeH1LagrangeVar*>(i->second);
          ptFe->SetIsoOrder(isoOrder+orderOffset_);
        }
      }
    }// IF: MappingType

    // Store mapping type for this region
    mappingType_[region] = mType;

  }

  void FeSpaceH1Nodal::CheckConsistency(){
    //just set lobatto integration with element order for each spectral region
    std::set<RegionIdType>::iterator spIt = spectralRegions_.begin();

    while(spIt != spectralRegions_.end()){

      // get region node
      std::string regionName = ptGrid_->regionData[*spIt].name;
      PtrParamNode regionNode = infoNode_->Get("regionList")->Get(regionName);
      IntegOrder order;
      // UInt test = *spIt;
      //every reference element has the same order
      order.SetIsoOrder( refElems_[*spIt][Elem::ET_LINE2]->GetIsoOrder() );
      SetRegionIntegration( *spIt,IntScheme::LOBATTO, order, INTEG_MODE_ABSOLUTE,
                            regionNode );
      // test = 0;
      spIt++;
    }

  }


  void FeSpaceH1Nodal::ReadCustomAttributes(PtrParamNode node, RegionIdType region){
    LOG_DBG(feSpaceH1Nodal) << "Processing polyRegioNode for region " << region;

    bool spectral = node->Get("spectral",ParamNode::EX)->As<bool>();

    if(spectral)
      spectralRegions_.insert(region);
  }

  void FeSpaceH1Nodal::SetDefaultElements(PtrParamNode infoNode ){
    //but it could be, that the PDE requires a minimum order of elements..
    ApproxOrder order (ptGrid_->GetDim());
    order.SetIsoOrder(1);

    if(orderOffset_>0){
      order.SetIsoOrder(orderOffset_);
      SetRegionElements(ALL_REGIONS,POLYNOMIAL,order,infoNode);
    }else{
      //now we are pretty sure that we need a grid mapping
      SetRegionElements(ALL_REGIONS,GRID,order,infoNode);
    }
  }

  //! sets the default integration scheme and order
  void FeSpaceH1Nodal::SetDefaultIntegration(PtrParamNode infoNode ){
    regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].order.SetIsoOrder(0);
    regionIntegration_[ALL_REGIONS].mode = INTEG_MODE_RELATIVE;
  }
  
  bool FeSpaceH1Nodal::IsSameEntityApproximation( shared_ptr<EntityList> list,
                                                  shared_ptr<FeSpace> space ) {
    if( this->GetSpaceType()  != space->GetSpaceType()  ) {
      return false;
    }
    if( this->IsHierarchical() != space->IsHierarchical()) {
      return false;
    }
    
    // Cast other space to same type
    shared_ptr<FeSpaceH1Nodal> otherSpace = dynamic_pointer_cast<FeSpaceH1Nodal>(space);
    
    EntityList::ListType actListType = list->GetType();
    if ( ! (actListType == EntityList::ELEM_LIST) &&
        ! (actListType == EntityList::SURF_ELEM_LIST) &&
        ! (actListType == EntityList::NC_ELEM_LIST))  {
      return true;
    }
    
    // Loop over all elements
    EntityIterator it = list->GetIterator();

    // switch depending on mapping type
    for( it.Begin(); !it.IsEnd(); it++) {
      if( mappingType_[it.GetElem()->regionId] == GRID ) {
        FeH1LagrangeExpl * myElem = static_cast<FeH1LagrangeExpl*>(this->GetFe(it));
        FeH1LagrangeExpl * otherElem = static_cast<FeH1LagrangeExpl*>(otherSpace->GetFe(it));
        if( !( *myElem == *otherElem) ) {
          return false;
        }
      } else {
        FeH1LagrangeVar * myElem = static_cast<FeH1LagrangeVar*>(this->GetFe(it));
        FeH1LagrangeVar * otherElem = static_cast<FeH1LagrangeVar*>(otherSpace->GetFe(it));
        if( !( *myElem == *otherElem) ) {
          return false;
        }
      }
    }
    return true;
  }

  void FeSpaceH1Nodal::MapNodalBCs(){
    LOG_TRACE(feSpaceH1Nodal) << "Mapping Nodal BCs";
    StdVector<UInt> actNodes;

    // check if feFunction was defined
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
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
          LOG_DBG(feSpaceH1Nodal) << "\tHDBC for vNode " << vNode << ", dof " << *dofIt;
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
        // check first, if this node was already processed
        // loop over all dofs
        std::set<UInt>::const_iterator dofIt = (*actIBC)->dofs.begin();
        for( ; dofIt != (*actIBC)->dofs.end(); ++dofIt) { 
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
      UInt masterDof = (*actConstr)->masterDof;
      UInt slaveDof = (*actConstr)->slaveDof;
      GetNodesOfEntities(slaveNodes,(*actConstr)->slaveEntities);
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

}//namespace
