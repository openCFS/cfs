// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FeSpaceL2Nodal.cc
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

#include <boost/tokenizer.hpp>
#include "FeSpaceL2Nodal.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(feSpaceL2Nodal, "feSpaceL2Nodal")

namespace CoupledField{

FeSpaceL2Nodal::FeSpaceL2Nodal(PtrParamNode aNode, PtrParamNode infoNode,
                               Grid* ptGrid)
    : FeSpaceL2(aNode, infoNode, ptGrid ){
  type_ = L2;
  isHierarchical_ = false;
  mapType_ = POLYNOMIAL;
  polyType_ = LAGRANGE;

  infoNode_ = infoNode->Get("l2Nodal");
}


FeSpaceL2Nodal::~FeSpaceL2Nodal(){
  std::map< RegionIdType, std::map<Elem::FEType, FeH1* > >::iterator regionIt;
  regionIt = refElems_.begin();
  for( ; regionIt != refElems_.end(); ++regionIt ) {
    std::map<Elem::FEType, FeH1* > & elems = regionIt->second;
    std::map<Elem::FEType, FeH1* >::iterator elemIt = elems.begin();
    for( ; elemIt != elems.end(); ++elemIt ) {
      delete elemIt->second;
    }
  }
}

void FeSpaceL2Nodal::Init( shared_ptr<SolStrategy> solStrat ) {
  solStrat_ = solStrat;
  ReadIntegList();
  ReadPolyList();
}

BaseFE* FeSpaceL2Nodal::GetFe( const EntityIterator ent ,
                               IntScheme::IntegMethod& method,
                               IntegOrder & order  ) {
  BaseFE * ret = GetFe(ent);

  // Set correct integration order
  RegionIdType eRegion;// =  ent.GetElem()->regionId;
  if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
    eRegion = ent.GetSurfElem()->ptVolElems[0]->regionId;
  } else {
    eRegion = ent.GetElem()->regionId;
  }

  this->GetIntegration(ret, eRegion, method, order);
  // Note: The order is currently more or less hard-coded for isotropic order
  //intScheme->SetOrder( method, order[0][0] );

  return ret;

}

BaseFE* FeSpaceL2Nodal::GetFe( const EntityIterator ent ){

  if(ent.GetType() != EntityList::ELEM_LIST &&
      ent.GetType() != EntityList::SURF_ELEM_LIST &&
      ent.GetType() != EntityList::NC_ELEM_LIST){
    EXCEPTION("This version of GetFe expects a element iterator")
  }


  // Note: if the element is a surface element, we must omit the regionId
  // and look for the neigbor. Which one to take? Well, we had the
  // discussion already ....
  RegionIdType eRegion = NO_REGION_ID;
  if( ent.GetType() == EntityList::SURF_ELEM_LIST ) {
    eRegion = ent.GetSurfElem()->ptVolElems[0]->regionId;
  } else {
    eRegion = ent.GetElem()->regionId;
  }

  //Check if the region is there, otherwise fall back to default
  if(refElems_.find(eRegion) == refElems_.end()){
    eRegion = ALL_REGIONS;
  }

  if(refElems_[eRegion].find(ent.GetElem()->type) == refElems_[eRegion].end()){
    EXCEPTION("fespacel2::getfe( const entityiterator): requested fetype which is noch supported by space");
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

  // No need to set the order here, as this is already done once and for all in the
  // SetMapType() method. For higher order spaces with non-uniform polynomial order, this necessary.
  // myFe->SetIsoOrder( isoOrder_);

  return myFe;
}

BaseFE* FeSpaceL2Nodal::GetFe( UInt elemNum ){
  shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
  assert(feFct);
  const Elem * ptElem = feFct->GetGrid()->GetElem(elemNum);
  RegionIdType eRegion = ptElem->regionId;

  //Check if the region is there, otherwise fall back to default
  if(refElems_.find(eRegion) == refElems_.end()){
    eRegion = ALL_REGIONS;
  }

  if(refElems_[eRegion].find(ptElem->type) == refElems_[eRegion].end()){
    EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
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

void FeSpaceL2Nodal::PreCalcShapeFncs(){
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
void FeSpaceL2Nodal::Finalize(){
  /* Basic idea:
   * 1. Set correct element order at elements
   * 2. Create the VirtualNode Array
   * 3. Map boundary conditions
   * 4. Map equations only based on the virtualNodeArray
   */

  //      UInt numEqns_ = 0;
  //      UInt numFreeEquations_ = 0;
  //the damn thing is that we need to map the edges right here if we have internal surfaces
  //ok, should not hurt anyway
  shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
  assert(feFct);
  if(interiorElemMap_.size() > 0){
    feFct->GetGrid()->MapEdges();
    feFct->GetGrid()->MapFaces();
    CreateSurfaceElems();
  }

  CreateVirtualNodes();
  //Determine boundary Unknowns
  MapNodalBCs();
  MapNodalEqns(1);
  MapNodalEqns(2);


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

void FeSpaceL2Nodal::SetRegionElements(RegionIdType region,
                                       MappingType mType,
                                       const ApproxOrder& order,
                                       PtrParamNode infoNode ){

  LOG_DBG(feSpaceL2Nodal) << "Setting region elements";
  LOG_DBG2(feSpaceL2Nodal) << "\tegion: "
      << ptGrid_->GetRegion().ToString(region);
  LOG_DBG2(feSpaceL2Nodal) << "\tmappingType: " << mType;
  LOG_DBG2(feSpaceL2Nodal) << "\torder: " << order.ToString();

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
    refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLine1();
    refElems_[region][Elem::ET_TRIA3]  = new FeH1LagrangeTria1();
    refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuad1();
    refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHex1();
    refElems_[region][Elem::ET_WEDGE6] = new FeH1LagrangeWedge1();
    refElems_[region][Elem::ET_PYRA5]  = new FeH1LagrangePyra1();
    refElems_[region][Elem::ET_TET4]  = new FeH1LagrangeTet1();
    refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLine2();
    refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
    refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex2();
    refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
    refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra2();
    refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet2();

    UInt gridOrder = ptGrid_->IsQuadratic() ? 2 : 1;
    infoNode->Get("order")->SetValue(gridOrder);

  } else if (mType == POLYNOMIAL) {
    refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
    refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
    refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
    refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLineVar();
    refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
    refElems_[region][Elem::ET_QUAD9]  = new FeH1LagrangeQuadVar();
    refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHexVar();
    refElems_[region][Elem::ET_HEXA27] = new FeH1LagrangeHexVar();

    UInt isoOrder = order.GetIsoOrder();
    
    std::map<Elem::FEType, FeH1* >::iterator i = refElems_[region].begin();
    for( ; i != refElems_[region].end(); ++i ) {
      FeH1LagrangeVar * ptFe = dynamic_cast<FeH1LagrangeVar*>(i->second);
      ptFe->SetIsoOrder(isoOrder+orderOffset_);
    }
    mapType_ = POLYNOMIAL;
    infoNode->Get("order")->SetValue(isoOrder);

    switch(isoOrder)
    {
    case 1:
      refElems_[region][Elem::ET_TRIA3]  = new FeH1LagrangeTria1();
      refElems_[region][Elem::ET_TRIA6]  = new FeH1LagrangeTria1();

      refElems_[region][Elem::ET_TET4]  = new FeH1LagrangeTet1();
      refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet1();

      refElems_[region][Elem::ET_WEDGE6] = new FeH1LagrangeWedge1();
      refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge1();
      refElems_[region][Elem::ET_WEDGE18] = new FeH1LagrangeWedge1();

      refElems_[region][Elem::ET_PYRA5]  = new FeH1LagrangePyra1();
      refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra1();
      refElems_[region][Elem::ET_PYRA14] = new FeH1LagrangePyra1();
      break;
    case 2:
      refElems_[region][Elem::ET_TRIA6]  = new FeH1LagrangeTria2();

      refElems_[region][Elem::ET_TET10]  = new FeH1LagrangeTet2();

      // ET_WEDGE15 elements are not compatible with tensor product hexas.
      //refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
      refElems_[region][Elem::ET_WEDGE18] = new FeH1LagrangeWedge18();

      // ET_PYRA13 elements are not compatible with tensor product hexas.
      //refElems_[region][Elem::ET_PYRA13] = new FeH1LagrangePyra2();
      refElems_[region][Elem::ET_PYRA14] = new FeH1LagrangePyra14();
      break;
    default:
      break;
    }

  }

  // Store mapping type for this region
  mappingType_[region] = mType;

}

void FeSpaceL2Nodal::CheckConsistency(){
  //just set lobatto integration with element order for each spectral region
  std::set<RegionIdType>::iterator spIt = spectralRegions_.begin();

  while(spIt != spectralRegions_.end()){

    // get region node
    std::string regionName = ptGrid_->regionData[*spIt].name;
    PtrParamNode regionNode = infoNode_->Get("regionList")->Get(regionName);
    IntegOrder order;
    //every reference element has the same order
    order.SetIsoOrder( refElems_[*spIt][Elem::ET_LINE2]->GetIsoOrder() );
    SetRegionIntegration( *spIt,IntScheme::LOBATTO, order, INTEG_MODE_ABSOLUTE,
                          regionNode );
    spIt++;
  }

}


void FeSpaceL2Nodal::ReadCustomAttributes(PtrParamNode node, RegionIdType region){
  LOG_DBG(feSpaceL2Nodal) << "Processing polyRegioNode for region " << region;

  bool spectral = node->Get("spectral",ParamNode::EX)->As<bool>();

  if(spectral)
    spectralRegions_.insert(region);
}

void FeSpaceL2Nodal::SetDefaultElements(PtrParamNode infoNode ){
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
void FeSpaceL2Nodal::SetDefaultIntegration(PtrParamNode infoNode ){
  regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
  regionIntegration_[ALL_REGIONS].order.SetIsoOrder( 0 );
  regionIntegration_[ALL_REGIONS].mode = INTEG_MODE_RELATIVE;
}

bool FeSpaceL2Nodal::IsSameEntityApproximation( shared_ptr<EntityList> list,
                                             shared_ptr<FeSpace> space ) {
  
  if( this->GetSpaceType()  != space->GetSpaceType()  ) {
    return false;
  }
  if( this->IsHierarchical() != space->IsHierarchical()) {
    return false;
  }
  
  // Cast other space to same type
  shared_ptr<FeSpaceL2Nodal> otherSpace = dynamic_pointer_cast<FeSpaceL2Nodal>(space);
  
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

}

