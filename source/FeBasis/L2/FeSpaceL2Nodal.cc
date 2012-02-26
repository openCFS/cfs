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

#include "FeSpaceL2Nodal.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

DECLARE_LOG(feSpaceL2Nodal)
DEFINE_LOG(feSpaceL2Nodal, "feSpaceL2Nodal")


namespace CoupledField{

FeSpaceL2Nodal::FeSpaceL2Nodal(PtrParamNode aNode, PtrParamNode infoNode )
    : FeSpaceL2(aNode, infoNode){
  type_ = L2;
  isHierarchical_ = false;
  mapType_ = POLYNOMIAL;
  polyType_ = LAGRANGE;

  infoNode_ = infoNode->Get("l2Nodal");
}


FeSpaceL2Nodal::~FeSpaceL2Nodal(){

}

void FeSpaceL2Nodal::Init( shared_ptr<SolStrategy> solStrat ) {
  solStrat_ = solStrat;
  ReadIntegList();
  ReadPolyList();
}

BaseFE* FeSpaceL2Nodal::GetFe( const EntityIterator ent,
                               shared_ptr<IntScheme>& intScheme ) {
  BaseFE * ret = GetFe(ent);

  // Set correct integration order
  RegionIdType eRegion;// =  ent.GetElem()->regionId;
  if( ent.GetType() == EntityList::SURF_ELEM_LIST) {
    eRegion = ent.GetSurfElem()->ptVolElems[0]->regionId;
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
    EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
  }
  BaseFE * myFe = refElems_[eRegion][ent.GetElem()->type];

  // No need to set the order here, as this is already done once and for all in the
  // SetMapType() method. For higher order spaces with non-uniform polynomial order, this necessary.
  // myFe->SetIsoOrder( isoOrder_);

  return myFe;
}

BaseFE* FeSpaceL2Nodal::GetFe( UInt elemNum ){
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

  return myFe;

}



UInt FeSpaceL2Nodal::GetEntityOrder( UInt elemNum, BaseFE::EntityType type,
                                     UInt entityNum, UInt comp  ) {

  // For the Lagrangian space this is a trivial implementation, as we only have nodal unknowns
  // which have just 1 unkown per node.
  if( type == BaseFE::NODE ) {
    return 1;
  } else {
    return 0;
  }
}

void FeSpaceL2Nodal::PreCalcShapeFncs(){
  //now pre-calculate shape functions for the
  // desired element orders
  //get grip of the integrationScheme object

  // leave, if element space is hierarchical
  if (isHierarchical_)return;
  shared_ptr<IntScheme> integScheme = feFunction_->GetGrid()->GetIntegrationScheme();

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
        Matrix<Integer> curOrder;
        IntScheme::IntegMethod curMethod;
        GetIntegration(elemIt->second,*integIter,curMethod,curOrder);
        integScheme->GetIntegrationPoints(integPoints,shape,curOrder[0][0],curMethod);
        dynamic_cast<FeH1*>(elemIt->second)->SetFunctionsAtIp(integPoints);
        integIter++;
      }
      elemIt++;
    }
    regIt++;
  }
}


UInt FeSpaceL2Nodal::GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type,
                                        UInt entityNum  ) {

  // For the Lagrangian space this is a trivial implementation, as we only have nodal unknowns
  // which have just 1 unkown per node.
  if( type == BaseFE::NODE ) {
    return 1;
  } else {
    return 0;
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
  if(interiorElemMap_.size() > 0){
    this->feFunction_->GetGrid()->MapEdges();
    this->feFunction_->GetGrid()->MapFaces();
    CreateSurfaceElems();
  }

  CreateVirtualNodes();
  //Determine boundary Unknowns
  MapNodalBCs();
  MapNodalEqns(1);
  MapNodalEqns(2);


  CheckConsistency();
  isFinalized_ = true;
}

void FeSpaceL2Nodal::SetRegionElements(RegionIdType region,
                                       MappingType mType,
                                       const Matrix<Integer>& order,
                                       PtrParamNode infoNode ){

  LOG_DBG(feSpaceL2Nodal) << "Setting region elements";
  LOG_DBG2(feSpaceL2Nodal) << "\tegion: "
      << domain->GetGrid()->GetRegion().ToString(region);
  LOG_DBG2(feSpaceL2Nodal) << "\tmappingType: " << mType;
  LOG_DBG2(feSpaceL2Nodal) << "\torder: " << order[0][0];

  //This method may not be called after the space is finalized!
  if(isFinalized_){
    Exception("FeSpace::SetRegionMapping is called after finalization");
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

    UInt gridOrder = domain->GetGrid()->IsQuadratic() ? 2 : 1;
    infoNode->Get("order")->SetValue(gridOrder);

  } else if (mType == POLYNOMIAL) {
    refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
    refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
    refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
    refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLineVar();
    refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
    refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHexVar();

    //now set the order
    if(order.GetNumCols() != 1 || order.GetNumRows() != 1){
      Exception("feSpaceL2Nodal::SetRegionMapping : The order matrix may have only one entry for lagrange elements");
    }
    std::map<Elem::FEType, FeH1* >::iterator i = refElems_[region].begin();
    for( ; i != refElems_[region].end(); ++i ) {
      FeH1LagrangeVar * ptFe = dynamic_cast<FeH1LagrangeVar*>(i->second);
      ptFe->SetIsoOrder(order[0][0]+orderOffset_);
    }
    mapType_ = POLYNOMIAL;
    infoNode->Get("order")->SetValue(order[0][0]);
  }

  // print information to info.xml


}

void FeSpaceL2Nodal::CheckConsistency(){
  //just set lobatto integration with element order for each spectral region
  std::set<RegionIdType>::iterator spIt = spectralRegions_.begin();

  while(spIt != spectralRegions_.end()){

    // get region node
    std::string regionName = domain->GetGrid()->regionData[*spIt].name;
    PtrParamNode regionNode = infoNode_->Get("regionList")->Get(regionName);
    Matrix<Integer> order(1,1);
    //every reference element has the same order
    order[0][0] = refElems_[*spIt][Elem::ET_LINE2]->GetIsoOrder();
    SetRegionIntegration( *spIt,IntScheme::LOBATTO, order, ABSOLUTE,
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
  Matrix<Integer> order(1,1);
  if(orderOffset_>0){
    order[0][0] = orderOffset_;
    SetRegionElements(ALL_REGIONS,POLYNOMIAL,order,infoNode);
  }else{
    //now we are pretty sure that we need a grid mapping
    SetRegionElements(ALL_REGIONS,GRID,order,infoNode);
  }
}

//! sets the default integration scheme and order
void FeSpaceL2Nodal::SetDefaultIntegration(PtrParamNode infoNode ){
  regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
  regionIntegration_[ALL_REGIONS].order = Matrix<Integer>(1,1);
  regionIntegration_[ALL_REGIONS].order[0][0] = 0;
  regionIntegration_[ALL_REGIONS].mode = RELATIVE;
}

}

