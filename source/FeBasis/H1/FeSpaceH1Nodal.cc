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


  DECLARE_LOG(feSpaceH1Nodal)
  DEFINE_LOG(feSpaceH1Nodal, "feSpaceH1Nodal")
namespace CoupledField{

    //! Constructor
    FeSpaceH1Nodal::FeSpaceH1Nodal(PtrParamNode aNode, 
                                   PtrParamNode infoNode) 
    : FeSpaceH1(aNode, infoNode) {
      type_ = H1;
      isHierarchical_ = false;
      mapType_ = GRID;
      polyType_ = LAGRANGE;
      
      infoNode_ = infoNode->Get("h1Nodal");
    }

    //! Destructor
    FeSpaceH1Nodal::~FeSpaceH1Nodal(){
    }
    
    void FeSpaceH1Nodal::Init() {
      // read order of function space
      // read, if map type should be isotropic
      //TODO: Specify default integration method
      //       special case of spectral elements:
      //       If we have Gauss-LOBATTO, we set it to the element order
      //       (if the user did not specify something else).
      //       If we do not have Gauss-Lobatto we set to Gauss and order to 2
      //       (if the user did not specify something else).
      //read in polyLists and integLists for easier access later
      ReadIntegList();
      ReadPolyList();
    }
    
    BaseFE* FeSpaceH1Nodal::GetFe( const EntityIterator ent, 
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

    BaseFE* FeSpaceH1Nodal::GetFe( const EntityIterator ent ){

      if(ent.GetType() != EntityList::ELEM_LIST &&
          ent.GetType() != EntityList::SURF_ELEM_LIST){
        EXCEPTION("This version of GetFe expects a element iterator")
      }

      
      // Note: if the element is a surface element, we must omit the regionId
      // and look for the neigbor. Which one to take? Well, we had the 
      // discussion already ....
      RegionIdType eRegion = NO_REGION_ID;
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

      // No need to set the order here, as this is already done once and for all in the 
      // SetMapType() method. For higher order spaces with non-uniform polynomial order, this necessary.
      // myFe->SetIsoOrder( isoOrder_);

      return myFe;
    }

    BaseFE* FeSpaceH1Nodal::GetFe( UInt elemNum ){
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
    
    
    
    UInt FeSpaceH1Nodal::GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                              UInt entityNum, UInt comp  ) {
      
      // For the Lagrangian space this is a trivial implementation, as we only have nodal unknowns
      // which have just 1 unkown per node.
      if( type == BaseFE::NODE ) {
        return 1;
      } else {
        return 0;
      }
    }
    
    UInt FeSpaceH1Nodal::GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                                 UInt entityNum  ) {

      // For the Lagrangian space this is a trivial implementation, as we only have nodal unknowns
      // which have just 1 unkown per node.
      if( type == BaseFE::NODE ) {
        return 1;
      } else {
        return 0;
      }
    }

    void FeSpaceH1Nodal::PreCalcShapeFncs(){
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
          Matrix<Integer> curOrder;
          IntScheme::IntegMethod curMethod;
          GetIntegration(elemIt->second,regIt->first,curMethod,curOrder);
          integScheme->GetIntegrationPoints(integPoints,shape,curOrder[0][0],curMethod);
          dynamic_cast<FeH1*>(elemIt->second)->SetFunctionsAtIp(integPoints);
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
      //PrintEqnMap();
      CheckConsistency();
      isFinalized_ = true;
    }

    void FeSpaceH1Nodal::SetRegionElements(RegionIdType region, 
                                              MappingType mType,
                                              const Matrix<Integer>& order,
                                              PtrParamNode infoNode ){
      
      LOG_DBG(feSpaceH1Nodal) << "Setting region elements";
      LOG_DBG2(feSpaceH1Nodal) << "\tegion: "
          << domain->GetGrid()->GetRegion().ToString(region);
      LOG_DBG2(feSpaceH1Nodal) << "\tmappingType: " << mType;
      LOG_DBG2(feSpaceH1Nodal) << "\torder: " << order[0][0];
      
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
        refElems_[region][Elem::ET_LINE3]  = new FeH1LagrangeLine2();
        refElems_[region][Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
        refElems_[region][Elem::ET_HEXA20] = new FeH1LagrangeHex2();
        refElems_[region][Elem::ET_WEDGE15] = new FeH1LagrangeWedge2();
        
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
          Exception("FeSpaceH1Nodal::SetRegionMapping : The order matrix may have only one entry for lagrange elements");
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

    void FeSpaceH1Nodal::CheckConsistency(){
      //just set lobatto integration with element order for each spectral region
      std::set<RegionIdType>::iterator spIt = spectralRegions_.begin();

      while(spIt != spectralRegions_.end()){
      
        // get region node
        std::string regionName = domain->GetGrid()->GetRegion().ToString(*spIt);
        PtrParamNode regionNode = infoNode_->Get("regionList")->Get(regionName);
        Matrix<Integer> order(1,1);
        //every reference element has the same order
        order[0][0] = refElems_[*spIt][Elem::ET_LINE2]->GetIsoOrder();
        SetRegionIntegration( *spIt,IntScheme::LOBATTO, order, ABSOLUTE,
                              regionNode );
        spIt++;
      }

    }


    void FeSpaceH1Nodal::ReadCustomAttributes(PtrParamNode node, RegionIdType region){
      LOG_DBG(feSpaceH1Nodal) << "Processing polyRegioNode for region " << region;
      
      bool spectral = node->Get("spectral",ParamNode::EX)->As<bool>();
      
      if(spectral)
        spectralRegions_.insert(spectral);
    }

    void FeSpaceH1Nodal::SetDefaultElements(PtrParamNode infoNode ){
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
    void FeSpaceH1Nodal::SetDefaultIntegration(PtrParamNode infoNode ){
      regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
      regionIntegration_[ALL_REGIONS].order = Matrix<Integer>(1,1);
      regionIntegration_[ALL_REGIONS].order[0][0] = 2;
      regionIntegration_[ALL_REGIONS].mode = RELATIVE;
    }

 }//namespace
