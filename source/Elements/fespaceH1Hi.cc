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
  FeSpaceH1Hi::FeSpaceH1Hi(ParamNode* aNode)
  : FeSpaceH1(aNode) {
    mapType_ = POLYNOMIAL;
    type_ = H1;
    isHierarchical_ = true;
  }

  //! Destructor
  FeSpaceH1Hi::~FeSpaceH1Hi(){
  }
  
  //! Initialize class
  void FeSpaceH1Hi::Init() {
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

  void FeSpaceH1Hi::SetMapType( MappingType mapT){

    // check, that the correct map type is set:
    // For hierarchical elements we always need the 
    mapType_ = mapT;
    refElems_[Elem::ET_LINE2]  = new FeH1HiLine();
    refElems_[Elem::ET_QUAD4]  = new FeH1HiQuad();
    refElems_[Elem::ET_HEXA8]  = new FeH1HiHex();

    // Do we need the association of geometric "higher" 
    // elements to the ones here?
    //      refElems_[Elem::ET_LINE3]  = new FeH1LagrangeLine2();
    //      refElems_[Elem::ET_QUAD8]  = new FeH1LagrangeQuad2();
    //      refElems_[Elem::ET_HEXA20] = new FeH1LagrangeHex2();
  }


  UInt FeSpaceH1Hi::GetEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                      UInt entityNum, UInt comp ) {
    if( type == BaseFE::NODE) {
      return 1;
    } else if( type == BaseFE::EDGE ) {
      // Currently we support just isotropic order
      //StdVector<UInt> & orders = edgeOrder_[elemNum];
      
      return isoOrder_;
    } else if( type == BaseFE::FACE ) {
      // IMPLEMENT ME
    }
  }

  UInt FeSpaceH1Hi::GetMaxEntityOrder( UInt elemNum, BaseFE::EntityType type, 
                                         UInt entityNum ) {
    if( type == BaseFE::NODE) {
      return 1;
    } else if( type == BaseFE::EDGE ) {
      // Currently we support just isotropic order
      //StdVector<UInt> & orders = edgeOrder_[elemNum];
      //return *(std::max_element( orders.Begin(), orders.End() ));
      return isoOrder_;
    } else if( type == BaseFE::FACE ) {
      // IMPLEMENT ME
    }
  }

  //! Map equations i.e. intialize object
  void FeSpaceH1Hi::Finalize(){
    /* Basic idea:
     * 1. Create the VirtualNode Array
     * 2. Make the polynomial order of edges / faces consistent
     * 3. Map boundary conditions
     * 4. Map equations only based on the virtualNodeArray
     */

    UInt numEqns_ = 0;
    UInt numFreeEquations_ = 0;
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
    if(refElems_.find(ent.GetElem()->type) == refElems_.end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    BaseFE * myFe = refElems_[ent.GetElem()->type];


    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entitiy orders.
    myFe->SetIsoOrder( isoOrder_);

    return myFe;
  }
  
  BaseFE* FeSpaceH1Hi::GetFe( UInt elemNum ){
    const Elem * ptElem = feFunction_->GetGrid()->GetElem(elemNum); 
    if(refElems_.find(ptElem->type) == refElems_.end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    BaseFE * myFe = refElems_[ptElem->type];


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
} // end of namespace
