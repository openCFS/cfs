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
    FeSpaceH1Lagrange::FeSpaceH1Lagrange(ParamNode* aNode) 
    : FeSpaceH1(aNode) {
      mapType_ = GRID;
      type_ = H1;
      isHierarchical_ = false;
    }

    //! Destructor
    FeSpaceH1Lagrange::~FeSpaceH1Lagrange(){
    }
    
    void FeSpaceH1Lagrange::Init() {
      // read order of function space
      // read, if map type should be isotropic
      
      ParamNode * orderNode = NULL;
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
      }
    }

    void FeSpaceH1Lagrange::SetMapType( MappingType mapT){
      mapType_ = mapT;
      //build up the pointerMap
      
      // Note: at the moment this is not implemented in a 
      // clean fashion. The equation mapping currently relies
      // on the global node numbers of the element.
      // Instead, it should use the GetNumFncs(AnsatzFct::NODE)
      // method, which would deliver the correct number of
      // unknowns for the calculation element instead of the 
      // geometric element.

      if( mapType_ == GRID ) {
        refElems_[Elem::ET_LINE2]  = new FeH1LagrangeLine1();
        refElems_[Elem::ET_QUAD4]  = new FeH1LagrangeQuad1();
        refElems_[Elem::ET_HEXA8]  = new FeH1LagrangeHex1();
        refElems_[Elem::ET_LINE3]  = new FeH1LagrangeLine2();
        refElems_[Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
        refElems_[Elem::ET_HEXA20] = new FeH1LagrangeHex2();
      } else if (mapType_ == POLYNOMIAL) {
        refElems_[Elem::ET_LINE2]  = new FeH1LagrangeLineVar();
        refElems_[Elem::ET_QUAD4]  = new FeH1LagrangeQuadVar();
        refElems_[Elem::ET_HEXA8]  = new FeH1LagrangeHexVar();
        refElems_[Elem::ET_LINE3]  = new FeH1LagrangeLineVar();
        refElems_[Elem::ET_QUAD8]  = new FeH1LagrangeQuadVar();
        refElems_[Elem::ET_HEXA20] = new FeH1LagrangeHexVar();
      }
    }
 
    BaseFE* FeSpaceH1Lagrange::GetFe( const EntityIterator ent ){
      if(refElems_.find(ent.GetElem()->type) == refElems_.end()){
        EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
      }
      BaseFE * myFe = refElems_[ent.GetElem()->type];

      // No need to set the order here, as this is already done once and for all in the 
      // SetMapType() method. For higher order spaces with non-uniform polynomial order, this necessary.
      // myFe->SetIsoOrder( isoOrder_);

      return myFe;
    }

    BaseFE* FeSpaceH1Lagrange::GetFe( UInt elemNum ){
      const Elem * ptElem = feFunction_->GetGrid()->GetElem(elemNum); 
      if(refElems_.find(ptElem->type) == refElems_.end()){
        EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
      }
      BaseFE * myFe = refElems_[ptElem->type];

      // No need to set the order here, as this is already done once and for all in the 
      // SetMapType() method. For higher order spaces with non-uniform polynomial order, this necessary.
      // myFe->SetIsoOrder( isoOrder_);

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


      if (mapType_ != GRID ) {
        // for lagrangian elements only isoparametric element orders are supported
        std::map<Elem::FEType,BaseFE*>::iterator i = refElems_.begin();
        for( ; i != refElems_.end(); ++i ) {
          i->second->SetIsoOrder(isoOrder_);
        }
      }
      
      UInt numEqns_ = 0;
      UInt numFreeEquations_ = 0;
      CreateVirtualNodes();
      //Determine boundary Unknowns
      MapNodalBCs();
      MapNodalEqns(1);
      MapNodalEqns(2);
      
      // TEMPORARY: print information 
      //PrintEqnMap();
      isFinalized_ = true;
    }


    }
