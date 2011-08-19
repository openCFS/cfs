// =====================================================================================
// 
//       Filename:  fespaceH1Lagrange.cc
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
#include "Elements/fespaceH1Lagrange.hh"
//#include "Domain/elem.hh"
//#include "Elements/integrationScheme.hh"

namespace CoupledField{

    //! Constructor
    FeSpaceH1Lagrange::FeSpaceH1Lagrange(PtrParamNode aNode) 
    : FeSpaceH1(aNode) {
      type_ = H1;
      isHierarchical_ = false;
      mapType_ = GRID;
      polyType_ = LAGRANGE;
    }

    //! Destructor
    FeSpaceH1Lagrange::~FeSpaceH1Lagrange(){
    }
    
    void FeSpaceH1Lagrange::Init() {
      // read order of function space
      // read, if map type should be isotropic
      //TODO: Specify default integration method
      //       special case of spectral elements:
      //       If we have Gauss-LOBATTO, we set it to the element order
      //       (if the user did not specify something else).
      //       If we do not have Gauss-Lobatto we set to Gauss and order to 2
      //       (if the user did not specify something else).
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
      CreateRefElems();
    }

    BaseFE* FeSpaceH1Lagrange::GetFe( const EntityIterator ent ){

      if(ent.GetType() != EntityList::ELEM_LIST){
        EXCEPTION("This version of GetFe expects a element iterator")
      }

      RegionIdType eRegion = ent.GetElem()->regionId;

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

    BaseFE* FeSpaceH1Lagrange::GetFe( UInt elemNum ){
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
    
    UInt FeSpaceH1Lagrange::GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                              UInt entityNum, UInt comp  ) {
      
      // For the Lagrangian space this is a trivial implementation, as we only have nodal unknowns
      // which have just 1 unkown per node.
      if( type == BaseFE::NODE ) {
        return 1;
      } else {
        return 0;
      }
    }
    
    UInt FeSpaceH1Lagrange::GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
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
    void FeSpaceH1Lagrange::Finalize(){
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
      //PrintEqnMap();
      CheckConsistency();
      isFinalized_ = true;
    }

    void FeSpaceH1Lagrange::SetRegionElements(RegionIdType region, MappingType mType,Matrix<Integer> order){
      //This method may not be called after the space is finalized!
      if(isFinalized_){
        Exception("FeSpace::SetRegionMapping is called after finalization");
      }

      if(mType == GRID){
        refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLine1();
        refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuad1();
        refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHex1();
        refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLine2();
        refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
        refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex2();
      } else if (mType == POLYNOMIAL) {
        refElems_[region][Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
        refElems_[region][Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
        refElems_[region][Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
        refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLineVar();
        refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
        refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHexVar();

        //now set the order
        if(order.GetNumCols() != 1 || order.GetNumRows() != 1){
          Exception("FeSpaceH1Lagrange::SetRegionMapping : The order matrix may have only one entry for lagrange elements");
        }
        std::map<Elem::FEType, BaseFE* >::iterator i = refElems_[region].begin();
        for( ; i != refElems_[region].end(); ++i ) {
          i->second->SetIsoOrder(order[0][0]+orderOffset_);
        }
        mapType_ = POLYNOMIAL;
      }
    }


    void FeSpaceH1Lagrange::SetRegionIntegration(RegionIdType region, IntScheme::IntegMethod method, Matrix<Integer> order){
      //TODO:Implementation of defaults (ALL_REGIONS) and XML
      regionIntegration_[region].first = method;
      regionIntegration_[region].second = order;
    }

    void FeSpaceH1Lagrange::CheckConsistency(){
      //just set lobatto integration with element order for each spectral region
      std::map<RegionIdType,bool>::iterator spIt = spectralRegions_.begin();

      while(spIt != spectralRegions_.end()){
        Matrix<Integer> order(1,1);
        //every reference element has the same order
        order[0][0] = refElems_[spIt->first][Elem::ET_LINE2]->GetIsoOrder();
        SetRegionIntegration(spIt->first,IntScheme::LOBATTO,order);
        spIt++;
      }

    }


    void FeSpaceH1Lagrange::ProcessPolyRegionNode(PtrParamNode node, RegionIdType region){
      Matrix<Integer> order(1,1);
      order[0][0] = -1;
      PtrParamNode isoOrderNode = node->Get("isoOrder", ParamNode::PASS );
      bool spectral = node->Get("spectral",ParamNode::EX)->As<bool>();
      bool grid = node->Get("useGridOrder",ParamNode::EX)->As<bool>();

      spectralRegions_[region] = spectral;
      //determine mapping type
      MappingType curMap = POLYNOMIAL;
      if(!spectral && grid){
        curMap = GRID;
      }
      if(isoOrderNode){
        Integer isoOrder = isoOrderNode->As<Integer>();
        order[0][0] = isoOrder;
      }else{
        if(spectral){
          WARN("No order specified for spectral element region. setting it to 1");
        }
        order[0][0] = 1;
      }
      SetRegionElements(region,curMap,order);
    }

    void FeSpaceH1Lagrange::CreateDefaultElements(){
      //but it could be, that the PDE requires a minimum order of elements..
      Matrix<Integer> order(1,1);
      if(orderOffset_>0){
        order[0][0] = orderOffset_;
        SetRegionElements(ALL_REGIONS,POLYNOMIAL,order);
      }else{
        //now we are pretty sure that we need a grid mapping
        SetRegionElements(ALL_REGIONS,GRID,order);
      }
    }

    //! sets the default integration scheme and order
    void FeSpaceH1Lagrange::SetDefaultIntegration(){
      regionIntegration_[ALL_REGIONS].first = IntScheme::GAUSS;
      regionIntegration_[ALL_REGIONS].second = Matrix<Integer>(1,1);
      regionIntegration_[ALL_REGIONS].second[0][0] = -1;
    }

 }//namespace
