#include "FeSpaceHCurlHi.hh"

// boost graph related stuff
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


#include "HCurlElemsHi.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/algsys/SolStrategy.hh"



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

// declare class specific logging stream
DECLARE_LOG(feSpaceHCurlHi)
DEFINE_LOG(feSpaceHCurlHi, "feSpaceHCurlHi")
namespace CoupledField{

  //! Constructor
  FeSpaceHCurlHi::FeSpaceHCurlHi( PtrParamNode aNode, 
                                  PtrParamNode infoNode,
                                  Grid* ptGrid )
  : FeSpaceHi(aNode, infoNode, ptGrid ) {
    mapType_ = POLYNOMIAL;
    type_ = HCURL;
    isHierarchical_ = true;
    polyType_ = LEGENDRE;
    onlyLowestOrder_ = false;
    groupAnisoEdges_ = false;
    maxAspectRatio_ = 0.0;
    
    infoNode_ = infoNode->Get("hCurlHierarchical");
    
    // important: trigger mapping of edges / faces
    ptGrid_->MapEdges();
    ptGrid_->MapFaces();
  }

  //! Destructor
  FeSpaceHCurlHi::~FeSpaceHCurlHi(){
    std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > >::iterator regionIt;
    regionIt = refElems_.begin();
    for( ; regionIt != refElems_.end(); ++regionIt ) {
      std::map<Elem::FEType, FeHCurlHi* > & elems = regionIt->second;
      std::map<Elem::FEType, FeHCurlHi* >::iterator elemIt = elems.begin();
      for( ; elemIt != elems.end(); ++elemIt ) {
        delete elemIt->second;
      }
    }
    
    
    // delete also reference elements of 1st order (in case of 2 level approach)
    regionIt = refElems1St_.begin();
    for( ; regionIt != refElems1St_.end(); ++regionIt ) {
      std::map<Elem::FEType, FeHCurlHi* > & elems = regionIt->second;
      std::map<Elem::FEType, FeHCurlHi* >::iterator elemIt = elems.begin();
      for( ; elemIt != elems.end(); ++elemIt ) {
        delete elemIt->second;
      }
    }
  }

  void FeSpaceHCurlHi::Init( shared_ptr<SolStrategy> solStrat ) {
    
    solStrat_ = solStrat;
    // read order of function space
    // read, if map type should be isotropic

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
    //read in polyLists and integLists for easier access later
    ReadIntegList();
    ReadPolyList();
  }
  
//  void FeSpaceHCurlHi::SetStrategy( SolStrategyType strategy, UInt step ) {
//
//    solStrategy_ = strategy;
//    solStep_ = step;

//    // Check: If we have TWO_LEVEL strategy and we are in the 2nd step we 
//    // have to do the following:
//    if(solStrategy_ == STRAT_TWO_LEVEL && solStep_ == 1 ) {
//      
//      // in the 1st step we use hard-codedo order 0
//      //SetMapType(POLYNOMIAL);
//      //isoOrder_ = 0;
//      mapType_ = POLYNOMIAL;
//    }
//    
//    if( solStrategy_ == STRAT_TWO_LEVEL && solStep_ == 2 ) {
//      isFinalized_ = false;
//      numEqns_ = 0;
//      numFreeEquations_ = 0;
//      bcCounter_.clear();
//      gridToVirtualNodes_.clear();
//      nodes_.Clear();
//      nodesType_.clear();
//      virtualNodes_.clear();
//      nodeMap_.eqns.clear();
//      nodeMap_.BcKeys.clear();
//      nodeMap_.constraintNodes.clear();
//    
//    // Re-initialize structure
//      this->Init();
//      this->Finalize();
//    }
//  }


  void FeSpaceHCurlHi::SetRegionElements(RegionIdType region, 
                                         MappingType mType,
                                         const ApproxOrder& order,
                                         PtrParamNode infoNode ){
    
    LOG_DBG(feSpaceHCurlHi) 
        << "Set region elements for region " << ptGrid_->GetRegion().ToString(region)
        << ", order is " << order.ToString();
        
    //This method may not be called after the space is finalized!
    if(isFinalized_){
      Exception("FeSpace::SetRegionMapping is called after finalization");
    }

    //ToDo: save the information...
    // QUERY FOR USER PARAMS IS STILL TO COME
    refElems_[region][Elem::ET_QUAD4]  = new FeHCurlHiQuad();
    refElems_[region][Elem::ET_HEXA8]  = new FeHCurlHiHex();
    
    refElems_[region][Elem::ET_QUAD8]  = new FeHCurlHiQuad();
    refElems_[region][Elem::ET_HEXA20]  = new FeHCurlHiHex();


    refElems1St_[region][Elem::ET_QUAD4]  = new FeHCurlHiQuad();
    refElems1St_[region][Elem::ET_HEXA8]  = new FeHCurlHiHex();

    refElems1St_[region][Elem::ET_QUAD8]  = new FeHCurlHiQuad();
    refElems1St_[region][Elem::ET_HEXA20]  = new FeHCurlHiHex();

    
    // set order for every region
    SetRegionOrder( region, order );
    
    // set explicit first order of all elements
    std::map<Elem::FEType, FeHCurlHi* >::iterator i;
     i = refElems1St_[region].begin();
    for( ; i != refElems1St_[region].end(); ++i ) {
      i->second->SetIsoOrder(0);
    }
    
    infoNode->Get("order")->SetValue(order.ToString());

  }


  void FeSpaceHCurlHi::CheckConsistency(){

  }

  //! sets the default integration scheme and order
  void FeSpaceHCurlHi::SetDefaultIntegration(PtrParamNode infoNode ){
    regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].order.SetIsoOrder( 0 );
    regionIntegration_[ALL_REGIONS].mode = RELATIVE;
  }

  //! Map equations i.e. intialize object
  void FeSpaceHCurlHi::Finalize(){
    /* Basic idea:
     * 1. Create the VirtualNode Array
     * 2. Make the polynomial order of edges / faces consistent
     * 3. Map boundary conditions
     * 4. Map equations only based on the virtualNodeArray
     */
    AdjustEntityOrder();
    CreateVirtualNodes();

    // Map normal Bcs
    MapNodalBCs();
    
    // Fix higher order dofs in anisotropic case
    FixHigherOrderAnisoDofs();    
    MapNodalEqns(1);
    MapNodalEqns(2);
    
    isFinalized_ = true;
  }
  
  void FeSpaceHCurlHi::UpdateToSolStrategy() {
    if( solStrat_->GetType() == SolStrategy::TWO_LEVEL_STRATEGY &&
        solStrat_->GetNumSolSteps() == 2 &&
        solStrat_->GetActSolStep() == 1 ) {
      onlyLowestOrder_ = true;
      std::cerr << "=> only lowest order\n";
    } else {
      onlyLowestOrder_ = false;
    }
         
  }
  
  void FeSpaceHCurlHi::GetNodesOfElement( StdVector<UInt>& nodes,
                           const Elem* ptElem,
                           BaseFE::EntityType entType){
     UInt elemNum = ptElem->elemNum;
     if(virtualNodes_.find(elemNum) ==virtualNodes_.end()){

       EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element #"
           << ptElem->elemNum << " of region " 
           << ptGrid_->GetRegion().ToString(ptElem->regionId));
     }
     if(entType == BaseFE::ALL){
       
       if( onlyLowestOrder_) {
         nodes.Resize(virtualNodes_[elemNum][BaseFE::EDGE].vNodes.GetSize());
         const StdVector<UInt>& entNodes =  virtualNodes_[elemNum][BaseFE::EDGE].vNodes;
         for (UInt i = 0; i < entNodes.GetSize(); i++ ){
           nodes[i] =  entNodes[i];
         }
       } else {
         nodes.Resize(virtualNodes_[elemNum][BaseFE::VERTEX].vNodes.GetSize()+
                      virtualNodes_[elemNum][BaseFE::EDGE].vNodes.GetSize()+
                      virtualNodes_[elemNum][BaseFE::FACE].vNodes.GetSize()+
                      virtualNodes_[elemNum][BaseFE::INTERIOR].vNodes.GetSize());
         ElemVirtualNodes::iterator nodeIt = virtualNodes_[elemNum].begin();
         UInt c = 0;
         while(nodeIt !=virtualNodes_[elemNum].end()){

           StdVector<UInt> & entNodes =  nodeIt->second.vNodes;
           for (UInt i = 0; i < entNodes.GetSize(); i++ ){
             nodes[c++] =  entNodes[i];
           }
           nodeIt++;
         }
       } // lowest order
     }else{
       nodes.Resize(virtualNodes_[elemNum][entType].vNodes.GetSize());
       const StdVector<UInt>& entNodes =  virtualNodes_[elemNum][entType].vNodes;
       for (UInt i = 0; i < entNodes.GetSize(); i++ ){
         nodes[i] =  entNodes[i];
       }
     }
   }
  
  BaseFE* FeSpaceHCurlHi::GetFe( const EntityIterator ent ,
                                 IntScheme::IntegMethod& method,
                                 IntegOrder & order  ) {
    BaseFE * ret = GetFe(ent);

    // Set correct integration order
    RegionIdType eRegion = GetVolElem(ent.GetElem())->regionId;
       
    this->GetIntegration(ret, eRegion, method, order);
    // Note: The order is currently more or less hard-coded for isotropic order
    
    if( onlyLowestOrder_) {
      order = IntegOrder();
      order.SetIsoOrder( 2 );
    }

    return ret;

  }
  
  BaseFE* FeSpaceHCurlHi::GetFe( const EntityIterator ent ){
    RegionIdType eRegion = GetVolElem(ent.GetElem())->regionId;
    
    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ent.GetElem()->type) == refElems_[eRegion].end()){
      EXCEPTION("FeSpaceHCurlHi: requested fetype which is not supported by space");
    }

    FeHCurlHi * myFe = NULL;
    if( onlyLowestOrder_) {
      ApproxOrder order;
      order.SetIsoOrder(0);
      myFe = refElems1St_[eRegion][ent.GetElem()->type];
      // attention: here we do NOT apply the max/min rule, as we assume constant
      // element order for all elements
      SetElemOrder( ent.GetElem(), myFe, order, false );
    } else {
      // Fetch reference element and set correct order
      myFe = refElems_[eRegion][ent.GetElem()->type];
      std::map<RegionIdType,ApproxOrder>::iterator it = regionOrder_.find(eRegion);
      SetElemOrder( ent.GetElem(), myFe, it->second, true );
      myFe = refElems_[eRegion][ent.GetElem()->type];
    }

    // ToDo: Currently hard coded to isotropic order. Here we should generalize the 
    // setting of entity orders.
    assert (myFe);
    return myFe;
  }
  
  BaseFE* FeSpaceHCurlHi::GetFe( UInt elemNum ){
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    const Elem * ptElem = feFct->GetGrid()->GetElem(elemNum); 
    RegionIdType eRegion = ptElem->regionId;
    

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    if(refElems_[eRegion].find(ptElem->type) == refElems_[eRegion].end()){
      EXCEPTION("FeSpaceHCurlHi::getfe( const entityiterator): requested fetype which is noch supported by space");
    }
    FeHCurlHi * myFe = NULL;
    if( onlyLowestOrder_) {
      ApproxOrder order;
      order.SetIsoOrder(0);
      myFe = refElems1St_[eRegion][ptElem->type];
      // attention: here we do NOT apply the max/min rule, as we assume constant
      // element order for all elements
      SetElemOrder( ptElem, myFe, order, false );
    } else {
      // Fetch reference element and set correct order
      myFe = refElems_[eRegion][ptElem->type];
      std::map<RegionIdType,ApproxOrder>::iterator it = regionOrder_.find(eRegion);
      SetElemOrder( ptElem, myFe, it->second, true );
      myFe = refElems_[eRegion][ptElem->type];
    }
    return myFe;
  }

  void FeSpaceHCurlHi::GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                        std::map<UInt,StdVector<std::set<Integer> > >&
                                        minorBlocks ) {
    
    LOG_DBG(feSpaceHCurlHi) << "Performing OLAS mappings ...";
    
    // Check: If we have a "standard" solution strategy, just call the 
    // method of the base class
    if( solStrat_->GetType() == SolStrategy::STD_STRATEGY ) {
      FeSpace::GetOlasMappings(sbmBlocks, minorBlocks );
      return;
    }
    
    // check in addition, if we need the two-level approach
    if( solStrat_->GetType() != SolStrategy::TWO_LEVEL_STRATEGY ) 
      EXCEPTION("Solution strategy of type ' " 
                << SolStrategy::strategyType.ToString(solStrat_->GetType()) 
                << "' not implemented for HCurl of higher order.");
    
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    FeFctIdType fctId = feFct->GetFctId();
    
    // maintain of already used entities
    std::set<UInt> edges, faces, elems;
    
    // Resize sbm-blocks
    sbmBlocks.Resize(3);
    
    // ==========================================
    // I) SPECIAL TREATMENT OF ANISOTROPIC FACES
    // ==========================================
    
    std::map<UInt, UInt> faceElems;
    StdVector<StdVector<UInt> > faceGroups;

    // obtain list of faces which have to be grouped together
    GetThinFaceGroups( faceElems, faceGroups );
    
    PtrParamNode thinNode = infoNode_->Get("thinElements");
    
    UInt numGroups = faceGroups.GetSize();
    UInt totalEntries = 0;
    for( UInt i = 0; i < numGroups; ++i) {
      totalEntries += faceGroups[i].GetSize();
    }
    UInt avgGroupSize = 0;
    if( numGroups > 0 )
      avgGroupSize = totalEntries / numGroups;
    thinNode->Get("active")->SetValue(groupAnisoEdges_);
    thinNode->Get("maxAspectRatio")->SetValue(maxAspectRatio_);
    thinNode->Get("numGroups")->SetValue(numGroups);
    thinNode->Get("avgGroupSize")->SetValue(avgGroupSize);
    
    
    // Loop over all faceGroups
    for(UInt iGroup = 0; iGroup < numGroups; ++iGroup ) {
      std::set<Integer>  faceEqns;

//      std::cerr << "====== Face Group #" << iGroup << " =========\n";
      UInt numFacePerGroup = faceGroups[iGroup].GetSize();
      const StdVector<UInt> & faceNums = faceGroups[iGroup]; 
//      std::cerr << "got " << faceNums.GetSize() << " faces\n";
//      
      // Loop over all faces in specific group
      for( UInt iFace = 0; iFace < numFacePerGroup; ++iFace ) {

        // obtain element
        UInt faceNum = faceNums[iFace];
        UInt elemNum = faceElems[faceNum];
//        std::cerr << "faceNum is " << faceNum << std::endl;
//                std::cerr << "elemNum is " << elemNum << std::endl;
        const Elem * ptEl = ptGrid_->GetElem(elemNum);
        
        ElemVirtualNodes & vn = virtualNodes_[elemNum];
        const StdVector<UInt> & faceNodes = vn[BaseFE::FACE].vNodes;
        const StdVector<UInt> & faceOffsets = vn[BaseFE::FACE].offset;
        
        
        // look for the correct face in the array
        Integer locFace = ptEl->faces.Find( faceNum );
        assert(locFace > -1 );
        UInt pos = 0;
        for( Integer i = 0; i < locFace; ++i ) {
          pos+= faceOffsets[UInt(i)];
        }
//        std::cerr << "pos is " << pos << std::endl;
        // Loop over all faceNodes
        for(UInt iNode = pos; iNode < faceOffsets[locFace]+pos; ++iNode ) {
          LOG_DBG(feSpaceHCurlHi) << "treating facenode " << iNode << std::endl;
          // check, if face was already numbered
          if( faces.find(faceNodes[iNode]) == faces.end()) {
            LOG_DBG(feSpaceHCurlHi) << "=> inserting " 
                << nodeMap_[faceNodes[iNode]].GetSize()
                << " equations\n";
            faceEqns.insert(nodeMap_[faceNodes[iNode]].Begin(),
                            nodeMap_[faceNodes[iNode]].End() );
            faces.insert(faceNodes[iNode]);
          }
        } // loop over face nodes
      } // loop over faces within group
      
      // Now define group by all collected faceEqns 
      if( faceEqns.size()) {
        LOG_DBG(feSpaceHCurlHi) << "faceEqns has size " << faceEqns.size() 
                                                       << std::endl;
//        std::set<Integer>::const_iterator fIt = faceEqns.begin();
//        std::cerr << "For face group #" << iGroup << " eqns:\n\t";
//        for( ; fIt != faceEqns.end(); ++fIt) {
//          std::cerr << *fIt << ", ";
//        }
//        std::cerr << "\n";
        minorBlocks[1].Push_back(faceEqns);
      }
    } // loop over groups
         
    
    // ==========================================
    // II) NORMAL MAPPING OF EDGES / FACES, ETC.
    // ==========================================
    
    // Loop over all elements
    boost::unordered_map< UInt, ElemVirtualNodes >::iterator elemIt = virtualNodes_.begin();
    for( ; elemIt != virtualNodes_.end(); ++elemIt ) {
      const UInt elemNum = elemIt->first;
      const Elem * elem = ptGrid_->GetElem(elemNum);
      UInt dim = Elem::shapes[elem->type].dim;;
      LOG_DBG(feSpaceHCurlHi) << "\nDim of elem #" 
          << elemNum << ": " << dim << std::endl;
//      if (dim != 3 ) continue;
//          
      
      
      ElemVirtualNodes& vn = elemIt->second;

      // ===================
      //  1) SBM-Definition
      // ===================
      
      // loop over all edges -> block #0 (if gradients are disabled)
      const StdVector<UInt> & edgeNodes = vn[BaseFE::EDGE].vNodes;
      for(UInt i = 0; i < edgeNodes.GetSize(); ++i ) {
        sbmBlocks[0][fctId].insert(nodeMap_[edgeNodes[i]].Begin(),
                                   nodeMap_[edgeNodes[i]].End());
      } 

      // loop over all faces -> block #1
      StdVector<UInt> & faceNodes = vn[BaseFE::FACE].vNodes;
      LOG_DBG(feSpaceHCurlHi) << "\n\nfaceNodes has size "
                              << faceNodes.GetSize() << std::endl;
      for(UInt i = 0; i < faceNodes.GetSize(); ++i ) {
        sbmBlocks[1][fctId].insert(nodeMap_[faceNodes[i]].Begin(),
                                   nodeMap_[faceNodes[i]].End());
      }

      // collect all inner nodes -> block #2
      StdVector<UInt> & innerNodes = vn[BaseFE::INTERIOR].vNodes;
      LOG_DBG(feSpaceHCurlHi) << "innerNode has size " << innerNodes.GetSize() 
                              << std::endl;
      for(UInt i = 0; i < innerNodes.GetSize(); ++i ) {
        sbmBlocks[2][fctId].insert(nodeMap_[innerNodes[i]].Begin(),
                                   nodeMap_[innerNodes[i]].End());
      }

      // ======================
      //  2) Matrix Sub-Blocks
      // ======================

      StdVector<UInt> & faceOffsets = vn[BaseFE::FACE].offset;
      
      
      // b) number remaining faces
      // SBM block #1: Group per face

      LOG_DBG(feSpaceHCurlHi)<< "faceOffset of elem #" << elemIt->first 
                             << " is "<< faceOffsets.ToString() << std::endl;
      UInt pos = 0;
      // loop over faces
      for( UInt iFace = 0; iFace < faceOffsets.GetSize(); ++iFace ) {
        LOG_DBG(feSpaceHCurlHi) << "treating face " << iFace << std::endl;
        
        std::set<Integer>  faceEqns;
        
        // loop over faceNodes
        LOG_DBG(feSpaceHCurlHi) << "looping from " << pos << " to " 
                                << faceOffsets[iFace]+pos << std::endl;
        for(UInt iNode = pos; iNode < faceOffsets[iFace]+pos; ++iNode ) {
          LOG_DBG(feSpaceHCurlHi) << "treating facenode " << iNode << std::endl;
          // check, if face was already numbered
          if( faces.find(faceNodes[iNode]) == faces.end()) {
            LOG_DBG(feSpaceHCurlHi) << "=> inserting " 
                                    << nodeMap_[faceNodes[iNode]].GetSize()
                                    << " equations\n";
            faceEqns.insert(nodeMap_[faceNodes[iNode]].Begin(),
                            nodeMap_[faceNodes[iNode]].End() );
            faces.insert(faceNodes[iNode]);
          }
        }
        pos += faceOffsets[iFace];
        if( faceEqns.size()) {
          LOG_DBG(feSpaceHCurlHi) << "faceEqns has size " << faceEqns.size() 
                                  << std::endl;
          minorBlocks[1].Push_back(faceEqns);
        }
        
      }
      
      // SBM block #2: interior equations per element
      std::set<Integer> innerEqns;
      LOG_DBG(feSpaceHCurlHi) << "size of innerNodes is " 
                              << innerNodes.GetSize() << std::endl;
      for(UInt i = 0; i < innerNodes.GetSize(); ++i ) {
        innerEqns.insert(nodeMap_[innerNodes[i]].Begin(),
                         nodeMap_[innerNodes[i]].End());
      }
      if( innerEqns.size() ) {
        minorBlocks[2].Push_back(innerEqns);
        LOG_DBG(feSpaceHCurlHi) << "innerEqns has size " 
                                << innerEqns.size() << std::endl;
      }
      
    }


  }

  void FeSpaceHCurlHi::GetThinFaceGroups( std::map<UInt,UInt>& faceElems,
                                          StdVector<StdVector<UInt> >&  fg ) {

    
    // only treat thin faces, if activated
    if( !groupAnisoEdges_)
      return;
    
    // use boost namespace to shorten thins a little
    using namespace boost;
    
    // define mapping from face-index (key) to face-number (value)
    std::map<UInt,UInt> i2f;
    std::map<UInt,UInt> f2i;
    
    // define graph type, describing the connectivity of the faces in terms
    // of their indices
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph faceGraph;

    // list of regions to look at
    StdVector<RegionIdType> regions;
    ptGrid_->GetVolRegionIds( regions );

    // maximum allowed aspect ratio
    Double minEdge, maxEdge;
    
    // Loop over all regions
    for( UInt iRegion = 0; iRegion < regions.GetSize(); ++iRegion ) {

      StdVector<Elem*> elems;
      ptGrid_->GetElems( elems, regions[iRegion] );
      UInt numElems = elems.GetSize();

      // Loop over all elements
      for( UInt iEl = 0; iEl < numElems; ++iEl) {
        Elem* ptEl = elems[iEl];
        shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap( ptEl );
        const ElemShape & shape = Elem::shapes[ptEl->type];
        
        // Determine aspect ratio, if within tolerance -> continue
        esm->GetMaxMinEdgeLength(maxEdge,minEdge);
        if( maxEdge/minEdge < maxAspectRatio_ ) 
          continue;

        // Determine extension of element in local directions
        Vector<Double> extension;
        esm->GetExtensionLocalDir( extension );
        UInt size = extension.GetSize();
        
        // sort extension vector in descending order, such that
        // the largest extension comes at first
        StdVector<UInt> indices(size);
        // initialize indices array
        for( UInt i = 0; i < size; i++ ) {
          indices[i] = i;
        }

        // -------------------------
        // Insertion sort algorithm
        // ------------------------
        UInt j;
        Double  comp;

        for( UInt i = 1; i < size; i++ ) {
          comp = extension[i];
          j = i;
          while( ( j > 0 ) && ( extension[j - 1] < comp ) ) {
            extension[j] = extension[j - 1];
            indices[j] = indices[j - 1];
            j = j - 1;
          }
          extension[j] = comp;
          indices[j] = i;
        }
        // -----------------------
        
        // determine direction, which get treated:
        // if we are here, the edge with the shortest direction (at position 
        // indices[2]) is the shortest and gets smoothed anyway.
        // But we also check the the second largest direction, which could be
        // also smaller than the aspect ratio.
        StdVector<UInt> minDirs;
        if( extension[0]/extension[2] > maxAspectRatio_ )
          minDirs.Push_back(indices[2]);
        if( extension[0]/extension[1] > maxAspectRatio_ )
          minDirs.Push_back(indices[1]);

        std::set<UInt> elemFaces;      
        for( UInt iDir = 0; iDir < minDirs.GetSize(); ++iDir ) {
          Integer minDir = minDirs[iDir];

          // Loop over all faces of this element orthogonal to "short"-direction
          for( UInt iFace = 0; iFace < ptEl->faces.GetSize(); ++iFace ) {

            if( shape.faceLocDirs[iFace][0] != minDir &&
                shape.faceLocDirs[iFace][1] != minDir ) {
              //            if( shape.faceLocDirs[iFace][0] == minDir ||
              //                shape.faceLocDirs[iFace][1] == minDir ) {
              elemFaces.insert( ptEl->faces[iFace]);
              faceElems[ptEl->faces[iFace]] = ptEl->elemNum;
            }
          } // loop element faces
        } // loop local directions

        
        std::set<UInt>::const_iterator faceIt1 = elemFaces.begin();
        std::set<UInt>::const_iterator faceIt2 = elemFaces.begin();
        // loop over faces (we need all pairs)
        Integer id1, id2; // continuous indices for faces
        
        for( ; faceIt1 != elemFaces.end(); ++faceIt1 ) {
          UInt face1 = *faceIt1;
          // Perform mapping of face (faceNum->index)
          if( f2i.find( face1) == f2i.end() ){
            id1 = i2f.size(); 
            i2f[id1] = face1;
            f2i[face1] = id1; 
          } else {
            id1 = f2i[face1];
          }

          for( faceIt2 = faceIt1; faceIt2 != elemFaces.end(); ++faceIt2 ) {
            UInt face2 = *faceIt2;  
            // Perform mapping of face (faceNum->index)
            if( f2i.find( face2 ) == f2i.end() ){
              id2 = i2f.size(); 
              i2f[id2] = face2;
              f2i[face2] = id2; 
            } else {
              id2 = f2i[face2];
            }
            // put every faceIndex-pair into the graph
            add_edge( id1, id2, faceGraph );
          } // loop elem faces (inner)
        } // loop elem faces (outer)
      } // loop elements
    } // loop over regions

    
    // Here we are finished setting up the graph. Now get all "connected"
    // faces using boost graph algorithms
    std::vector<int> groupNums(num_vertices(faceGraph));
    int numGroups = connected_components(faceGraph, &groupNums[0]);
    fg.Resize( numGroups );
    std::vector<int>::size_type i;
    //    std::cout << "Total number of components: " << numGroups << std::endl;
    for (i = 0; i != groupNums.size(); ++i) {
      Integer groupNum = groupNums[i];
      UInt faceNum = i2f[i]; 
      //      std::cout << "face " << faceNum <<" is in group " << groupNum << std::endl;
      fg[groupNum].Push_back(faceNum);
    }
    //    std::cout << std::endl;


    
//    // ----------------------------------------------------------------------
//    // Version a): Create named volume elements for connected faces
//    // ----------------------------------------------------------------------
//    // define named elements based on the groups
//    for( UInt i = 0; i < UInt(numGroups); ++i ) {
//      const StdVector<UInt> & faceNums = fg[i];
//      
//      StdVector<UInt> elemNums;
//      for( UInt iEl = 0; iEl < faceNums.GetSize(); ++ iEl ) {
//        UInt faceNum = faceNums[iEl];
//        elemNums.Push_back( faceElems[faceNum] );
//      }
//      std::string name = "fg_";
//      name += lexical_cast<std::string>(i+1);
//      ptGrid_->AddNamedElems(name, elemNums);
//      std::cerr << "creating named elems '" << name << "' on " << elemNums.Serialize() << std::endl;
//    }
    
    
//    // ----------------------------------------------------------------------
//    // Version b): Create surface elements for connected faces
//    // ----------------------------------------------------------------------
//    // define named surface elements based on the groups
//    
//    
//    UInt actElem = ptGrid_->GetNumElems()+1;
//    for( UInt i = 0; i < UInt(numGroups); ++i ) {
//      const StdVector<UInt> & faceNums = fg[i];
//      StdVector<UInt> elemNums;
//      for( UInt iFace = 0; iFace < faceNums.GetSize(); ++ iFace ) {
//        UInt faceNum = faceNums[iFace];
//        const Face & actFace = ptGrid_->GetFace( faceNum );
////        const Elem * ptEl = ptGrid_->GetElem( faceElems[faceNum] );
////        UInt locFace = ptEl->faces.Find(faceNum);
//        
//        StdVector<UInt> nodes = actFace.nodes;
//        //Face::GetSortedIndices( nodes, actFace.nodes, 4, ptEl->faceFlags[locFace]);
//        
//        
//        
//        std::cerr << "nodes are " << nodes.Serialize() << std::endl;
//        ptGrid_->AddElems(1);
//        ptGrid_->SetElemData( actElem, Elem::ET_QUAD4, NO_REGION_ID, &nodes[0] );
//        elemNums.Push_back(actElem++);
//      }
//      std::string name = "fg_";
//      name += lexical_cast<std::string>(i+1);
//      ptGrid_->AddNamedElems(name, elemNums);
//      std::cerr << "creating named elems '" << name << "' on " << elemNums.Serialize() << std::endl;
  }

  
  void FeSpaceHCurlHi::SetDefaultElements(PtrParamNode infoNode ){
    //but it could be, that the PDE requires a minimum order of elements...
    ApproxOrder order (ptGrid_->GetDim());
    order.SetIsoOrder(0);
    if(orderOffset_>0){
      order.SetIsoOrder(orderOffset_);
    }
    SetRegionElements(ALL_REGIONS,POLYNOMIAL,order, infoNode );
  }
  
  UInt FeSpaceHCurlHi::GetNumDofs() const {
    // As we have already vectorial basis functions, every
    // virtual "node" is basically just a scalar, so we
    // always return 1.
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    return 1;
  }
  
  void FeSpaceHCurlHi::TreatThinElements(Double maxAspectRatio ) {
    groupAnisoEdges_ = true;
    maxAspectRatio_ = maxAspectRatio;
    
  }
  
  FeHi* FeSpaceHCurlHi::GetFeHi( RegionIdType region, Elem::FEType type ) {
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
