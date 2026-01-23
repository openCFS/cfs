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

#include "def_use_openmp.hh"

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

  void FeSpaceHCurlHi::SetUseGradients(RegionIdType region) {
    LOG_DBG(feSpaceHCurlHi) << "Setting gradient for region"
        << ptGrid_->GetRegion().ToString(region) << " to true\n";
    useGradients_[region] = true;
  }
  

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
    
    // default: use no gradients
    useGradients_[region] = false;

    //ToDo: save the information...
    // QUERY FOR USER PARAMS IS STILL TO COME
    refElems_[region][Elem::ET_LINE2]  = new FeHCurlHiLine();
    refElems_[region][Elem::ET_TRIA3]  = new FeHCurlHiTria();
    refElems_[region][Elem::ET_QUAD4]  = new FeHCurlHiQuad();
    refElems_[region][Elem::ET_TET4]  = new FeHCurlHiTet();
    refElems_[region][Elem::ET_WEDGE6]  = new FeHCurlHiWedge();
    refElems_[region][Elem::ET_HEXA8]  = new FeHCurlHiHex();
    refElems_[region][Elem::ET_PYRA5]  = new FeHCurlHiPyra();
    
    refElems_[region][Elem::ET_TRIA6]  = new FeHCurlHiTria();
    refElems_[region][Elem::ET_QUAD8]  = new FeHCurlHiQuad();
    refElems_[region][Elem::ET_TET10]  = new FeHCurlHiTet();
    refElems_[region][Elem::ET_WEDGE15]  = new FeHCurlHiWedge();
    refElems_[region][Elem::ET_HEXA20]  = new FeHCurlHiHex();
    refElems_[region][Elem::ET_PYRA13]  = new FeHCurlHiPyra();

    refElems1St_[region][Elem::ET_TRIA3]  = new FeHCurlHiTria();
    refElems1St_[region][Elem::ET_QUAD4]  = new FeHCurlHiQuad();
    refElems1St_[region][Elem::ET_TET4]  = new FeHCurlHiTet();
    refElems1St_[region][Elem::ET_WEDGE6]  = new FeHCurlHiWedge();
    refElems1St_[region][Elem::ET_HEXA8]  = new FeHCurlHiHex();
    refElems1St_[region][Elem::ET_PYRA5]  = new FeHCurlHiPyra();

    refElems1St_[region][Elem::ET_TRIA6]  = new FeHCurlHiTria();
    refElems1St_[region][Elem::ET_QUAD8]  = new FeHCurlHiQuad();
    refElems1St_[region][Elem::ET_TET10]  = new FeHCurlHiTet();
    refElems1St_[region][Elem::ET_WEDGE15]  = new FeHCurlHiWedge();;
    refElems1St_[region][Elem::ET_HEXA20]  = new FeHCurlHiHex();
    refElems1St_[region][Elem::ET_PYRA13]  = new FeHCurlHiPyra();


    
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
  
  void FeSpaceHCurlHi::AdjustGradients( ) {
    
    LOG_DBG(feSpaceHCurlHi) << "Setting Region Gradients";

    // Loop over all regions and fix edge / face on
    // shared edges /faces
    std::map< RegionIdType, 
               std::map<Elem::FEType, FeHCurlHi* > >::const_iterator it;
    it = refElems_.begin();
    
    for( ; it != refElems_.end(); ++it ){
      StdVector<Elem*> elems;
      RegionIdType regionId = it->first;

      ptGrid_->GetElems( elems, regionId );
      UInt numElems = elems.GetSize();
      for( UInt i = 0; i < numElems; ++i ) {

        const Elem& el = *(elems[i]);
        LOG_DBG3(feSpaceHCurlHi) << "Treating element #" << el.elemNum;

        // Fetch reference element and set correct order
        FeHCurlHi * myFe = static_cast<FeHCurlHi*>(GetFeHi( regionId, el.type));

        // important: do not adjust the entity order here, as the
        // edge / face order is not yet initialized
        SetElemGrad( elems[i], myFe, regionId, false );

        // a) loop over all edges
        // -----------------------
        UInt numEdges = el.extended->edges.GetSize();
        LOG_DBG3(feSpaceHCurlHi) << "Checking " << numEdges << " edges ";
        const StdVector<bool>& edgeGradients = myFe->GetEdgeGradient();
        for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ){
          UInt edgeNum = std::abs(el.extended->edges[iEdge]);
          bool gradient = edgeGradients[iEdge];

          // check if edge got already mapped
          if( gradEdges_.find(edgeNum) == gradEdges_.end() ) {
            gradEdges_[edgeNum] = gradient;
            LOG_DBG3(feSpaceHCurlHi) << "\tedge #" << edgeNum << ", grad: " 
                << gradient ;
          } else {
            bool oldGrad = gradEdges_[edgeNum];
            LOG_DBG3(feSpaceHCurlHi) << "\tedge #" << edgeNum 
                << " -> Already set to " << oldGrad; 
            // check if previous set order is different from current
            if( gradient != oldGrad) {
              gradEdges_[edgeNum] = true;
              LOG_DBG3(feSpaceHCurlHi) << "\t\t=> Re-Set to " 
                  << gradEdges_[edgeNum]<< " (max. rule)";
              adjustedGradEdges_.insert( edgeNum );
            }
          }
        } // loop over edges


        // b) loop over all faces (only in 3D case )
        // -----------------------------------------
        if( ptGrid_->GetDim() == 3 ) {
          UInt numFaces = el.extended->faces.GetSize();
          LOG_DBG3(feSpaceHCurlHi) << "Checking " << numFaces << " faces ";
          const StdVector<bool>& faceGradients = myFe->GetFaceGradient();
          for( UInt iFace = 0; iFace < numFaces; ++iFace ){
            UInt faceNum = el.extended->faces[iFace];
            bool gradient = faceGradients[iFace];

            // check if face got already mapped
            if( gradFaces_.find(faceNum) == gradFaces_.end() ) {
              gradFaces_[faceNum] = gradient;
              LOG_DBG3(feSpaceHCurlHi) << "\tface #" << faceNum << ", grad: " 
                  << gradient;
            } else {
              bool oldGrad = gradFaces_[faceNum];
              LOG_DBG3(feSpaceHCurlHi) << "\tface #" << faceNum 
                  << " -> Already set to " << oldGrad;
              // check if previous set order is different from current
              if( gradient != oldGrad) {
                gradFaces_[faceNum] = true;
                LOG_DBG3(feSpaceHCurlHi) << "\t\t=> Re-Set to " 
                    << gradFaces_[faceNum] << " (max. rule)";
                adjustedGradFaces_.insert( faceNum );
              }
            }
          }  // loop over faces
        } // check for 3D case 

      } // loop over elements
    }
    
    // This part is called in the end.

    // The idea is as follows:
    // Up to now the polynomial orders of all edges and faces are stored
    // explicitly. However, we need to store only those, which
    // have been adjusted explicitly due to the maximum rule.
    // Thus we store only explicitly the edge and face order of those
    // edges / faces adjusted.

    // Create temporary map for edges
    boost::unordered_map<UInt, bool> edgeGrads;
    
    // Loop over all adjusted edges
    boost::unordered_set<UInt>::const_iterator edgeIt = adjustedGradEdges_.begin();
    for( ; edgeIt != adjustedGradEdges_.end(); ++edgeIt ) {
      edgeGrads[*edgeIt] = gradEdges_[*edgeIt];
    } // loop over edges

    // store only adjusted edges back
    gradEdges_ = edgeGrads;
    
    boost::unordered_map<UInt, bool> faceGrads;
    boost::unordered_set<UInt>::const_iterator faceIt = adjustedGradFaces_.begin();
      for( ; faceIt != adjustedGradFaces_.end(); ++faceIt ) {
        faceGrads[*faceIt] = gradFaces_[*faceIt];
      } // loop over edges

      // store only adjusted edges back
      gradFaces_ = faceGrads;
    
  }
  
  void FeSpaceHCurlHi::SetElemGrad( const Elem* ptEl, FeHCurlHi* ptFe,
                                    RegionIdType regionId,
                                    bool applyMaxRule ) {
    LOG_DBG(feSpaceHCurlHi) << "In SetElemGrad for elem " << ptEl->elemNum 
                        << ", maxRule: " << applyMaxRule;
     

    // Stage 1: Set gradient as given from the region template.
    ptFe->SetUseGradients(useGradients_[regionId]);
    LOG_DBG(feSpaceHCurlHi) << "\t->Gradient: " << 
        (useGradients_[ptEl->regionId] ? "true" : "false");

     // Stage 2: Adjust edge / face order according to minimum /
     // maximum rule, in case the element has neighboring elements with
     // different order. We will use the maximum rule
     if( applyMaxRule && gradEdges_.size() ) {
       // loop over all edges
       UInt numEdges = ptEl->extended->edges.GetSize();
       for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ) {
         UInt edgeNum = std::abs( ptEl->extended->edges[iEdge] );
         // check if edge got adjusted
         if( gradEdges_.find(edgeNum) != gradEdges_.end() ) {
           LOG_DBG3(feSpaceHCurlHi) << "Setting grad edge " << edgeNum
               << " to " << gradEdges_[edgeNum];
           ptFe->SetEdgeGradient(iEdge, gradEdges_[edgeNum]);
         }
       }

     }

     if( applyMaxRule && gradFaces_.size() && ptGrid_->GetDim() == 3  ) {
       // loop over all faces
       UInt numFaces = ptEl->extended->faces.GetSize();
       for( UInt iFace = 0; iFace < numFaces; ++iFace ) {
         UInt faceNum = ptEl->extended->faces[iFace];
         // check if face got adjusted
         if( gradFaces_.find(faceNum) != gradFaces_.end() ) {
           LOG_DBG3(feSpaceHCurlHi) << "Setting grad face " << faceNum
               << " to order " << gradFaces_[faceNum];
           ptFe->SetFaceGradient( iFace, gradFaces_[faceNum] );
         }
       }

     }
  }


  void FeSpaceHCurlHi::CheckConsistency(){

  }

  //! sets the default integration scheme and order
  void FeSpaceHCurlHi::SetDefaultIntegration(PtrParamNode infoNode ){
    regionIntegration_[ALL_REGIONS].method = IntScheme::GAUSS;
    regionIntegration_[ALL_REGIONS].order.SetIsoOrder( 0 );
    regionIntegration_[ALL_REGIONS].mode = INTEG_MODE_RELATIVE;
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
    
    // adjust the gradients
    AdjustGradients();
    
    CreateVirtualNodes();

    // Map normal Bcs
    MapNodalBCs();
    
    // Fix higher order dofs in anisotropic case
    FixHigherOrderAnisoDofs();    
    MapNodalEqns(1);
    MapNodalEqns(2);
    
#ifdef USE_OPENMP
    std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > >::iterator regIt = refElems_.begin();
    while(regIt != refElems_.end()){
      TL_refElems_[regIt->first] = regIt->second;
      ++regIt;
    }
    std::map< RegionIdType, std::map<Elem::FEType, FeHCurlHi* > >::iterator regIt1St = refElems1St_.begin();
    while(regIt1St != refElems1St_.end()){
      TL_refElems1St_[regIt1St->first] = regIt1St->second;
      ++regIt1St;
    }
#endif
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
    nodes.Clear();
    nodes.Reserve(30);
    // Use const references and .at() for thread-safe read-only access
    const EntityNodesType& eNodes = vNodesCont_[BaseFE::EDGE];
    const EntityNodesType& fNodes = vNodesCont_[BaseFE::FACE];
    const EntityNodesType& iNodes = vNodesCont_[BaseFE::INTERIOR];

    // Collect edge nodes
    {
      UInt numEdges = ptElem->extended->edges.GetSize();
      if( entType == BaseFE::EDGE || entType == BaseFE::ALL ) {

        for( UInt i = 0; i < numEdges; ++i ) {
          auto it = eNodes.find(std::abs(ptElem->extended->edges[i]));
          if( it != eNodes.end() ) {
            const StdVector<UInt>& edgeNodes = it->second;
            for( UInt j = 0; j < edgeNodes.GetSize(); ++j ) {
              nodes.Push_back(edgeNodes[j]);
            }
          }
        }
      }
    }

    // The inclusion of higher order terms might be deactivated
    // (as indicated by the onlyLowestOrder_ flag).
    if( !onlyLowestOrder_) {
      // Collect face nodes
      {
        UInt numFaces = ptElem->extended->faces.GetSize();
        if( entType == BaseFE::FACE || entType == BaseFE::ALL ) {

          for( UInt i = 0; i < numFaces; ++i ) {
            auto it = fNodes.find(std::abs(ptElem->extended->faces[i]));
            if( it != fNodes.end() ) {
              const StdVector<UInt>& faceNodes = it->second;
              for( UInt j = 0; j < faceNodes.GetSize(); ++j ) {
                nodes.Push_back(faceNodes[j]);
              }
            }
          }
        }
      }

      // Collect interior nodes
      {
        if( iNodes.size() ) {
          if( entType == BaseFE::INTERIOR || entType == BaseFE::ALL ) {
            auto it = iNodes.find(ptElem->elemNum);
            if( it != iNodes.end() ) {
              const StdVector<UInt>& intNodes = it->second;
              for( UInt j = 0; j < intNodes.GetSize(); ++j ) {
                nodes.Push_back(intNodes[j]);
              }
            }
          }
        }
      }
    } // if: !lowestOrder

    // Ensure, that at least one virtual node is present
    if( nodes.GetSize() == 0 ) { 
      EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element #"
          << ptElem->elemNum << " of region " 
          <<  ptGrid_->GetRegion().ToString(ptElem->regionId));
   }
  }
  
  BaseFE* FeSpaceHCurlHi::GetFe( const EntityIterator ent ,
                                 IntScheme::IntegMethod& method,
                                 IntegOrder & order  ) {
    BaseFE* ret = GetFe(ent);

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

    // Use .at() for thread-safe read-only access (throws if key missing)
    const auto& refElemsRegion = refElems_.at(eRegion);
    if(refElemsRegion.find(ent.GetElem()->type) == refElemsRegion.end()){
      EXCEPTION("FeSpaceHCurlHi: requested fetype which is not supported by space");
    }
#ifdef USE_OPENMP
    std::map<Elem::FEType, FeHCurlHi* >&  refMap = (isFinalized_ && omp_get_num_threads()>1)? TL_refElems_.at(eRegion).Mine() : refElems_.at(eRegion);
    std::map<Elem::FEType, FeHCurlHi* >&  refMap1St = (isFinalized_ && omp_get_num_threads()>1)? TL_refElems_.at(eRegion).Mine() : refElems1St_.at(eRegion);
#else
    std::map<Elem::FEType, FeHCurlHi* >&  refMap = refElems_.at(eRegion);
    std::map<Elem::FEType, FeHCurlHi* >&  refMap1St = refElems1St_.at(eRegion);
#endif

    FeHCurlHi * myFe = NULL;
    if( onlyLowestOrder_) {
      ApproxOrder order;
      order.SetIsoOrder(0);
      myFe = refMap1St[ent.GetElem()->type];
      // attention: here we do NOT apply the max/min rule, as we assume constant
      // element order for all elements
      SetElemOrder( ent.GetElem(), myFe, order, false );
    } else {
      // Fetch reference element and set correct order
      myFe = refMap[ent.GetElem()->type];
      std::map<RegionIdType,ApproxOrder>::iterator it = regionOrder_.find(eRegion);
      SetElemOrder( ent.GetElem(), myFe, it->second, true );
      SetElemGrad( ent.GetElem(), myFe, eRegion, true );
      myFe = refMap[ent.GetElem()->type];
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

    // Note: if the element is a surface element, we must omit the regionId
    // and look for the neigbor
    RegionIdType eRegion = GetVolElem(ptElem)->regionId;

    //Check if the region is there, otherwise fall back to default
    if(refElems_.find(eRegion) == refElems_.end()){
      eRegion = ALL_REGIONS;
    }

    // Use .at() for thread-safe read-only access (throws if key missing)
    const auto& refElemsRegion = refElems_.at(eRegion);
    if(refElemsRegion.find(ptElem->type) == refElemsRegion.end()){
      EXCEPTION("FeSpaceHCurlHi::getfe( const entityiterator): requested fetype which is noch supported by space");
    }

#ifdef USE_OPENMP
    std::map<Elem::FEType, FeHCurlHi* >&  refMap = (isFinalized_ && omp_get_num_threads()>1)? TL_refElems_.at(eRegion).Mine() : refElems_.at(eRegion);
    std::map<Elem::FEType, FeHCurlHi* >&  refMap1St = (isFinalized_ && omp_get_num_threads()>1)? TL_refElems_.at(eRegion).Mine() : refElems1St_.at(eRegion);
#else
    std::map<Elem::FEType, FeHCurlHi* >&  refMap = refElems_.at(eRegion);
    std::map<Elem::FEType, FeHCurlHi* >&  refMap1St = refElems1St_.at(eRegion);
#endif

    FeHCurlHi * myFe = NULL;
    if( onlyLowestOrder_) {
      ApproxOrder order;
      order.SetIsoOrder(0);
      myFe = refMap1St[ptElem->type];
      // attention: here we do NOT apply the max/min rule, as we assume constant
      // element order for all elements
      SetElemOrder( ptElem, myFe, order, false );
    } else {
      // Fetch reference element and set correct order
      myFe = refMap[ptElem->type];
      std::map<RegionIdType,ApproxOrder>::iterator it = regionOrder_.find(eRegion);
      SetElemOrder( ptElem, myFe, it->second, true );
      SetElemGrad( ptElem, myFe, eRegion, true );
      myFe = refMap[ptElem->type];
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
    
    // maintain of already mapped faces
    std::set<UInt> faces;
    
    // Resize sbm-blocks
    sbmBlocks.Resize(3);
    
    // ==========================================
    // I) SPECIAL TREATMENT OF ANISOTROPIC FACES
    // ==========================================

    // obtain list of faces which have to be grouped together
    StdVector<StdVector<UInt> > faceGroups;
    GetThinFaceGroups( faceGroups );
    
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
      LOG_DBG(feSpaceHCurlHi) << "treating group " << iGroup <<  std::endl;
      UInt numFacePerGroup = faceGroups[iGroup].GetSize();
      const StdVector<UInt> & faceNums = faceGroups[iGroup]; 

      // Loop over all faces in specific group
      for( UInt iFace = 0; iFace < numFacePerGroup; ++iFace ) {
        UInt faceNum = faceNums[iFace];
        LOG_DBG(feSpaceHCurlHi) << "\tFace # " << faceNum << std::endl;
        StdVector<UInt>&  faceNodes = vNodesCont_[BaseFE::FACE][faceNum];
        UInt numFaceNodes = faceNodes.GetSize();
        LOG_DBG(feSpaceHCurlHi) << "\tFaceNodes: " << faceNodes.Serialize() << std::endl;
        if( faces.find(faceNum) == faces.end()) {
          // Loop over all faceNodes
          for(UInt iNode = 0; iNode < numFaceNodes; ++iNode ) {
            UInt virtualNodeNum = faceNodes[iNode];
            LOG_DBG(feSpaceHCurlHi) << "treating facenode " << iNode << std::endl;
            // check, if face was already numbered
            LOG_DBG(feSpaceHCurlHi) << "=> inserting " 
                << nodeMap_[faceNodes[iNode]].GetSize()
                << " equations\n";
            faceEqns.insert(nodeMap_[virtualNodeNum].Begin(),
                            nodeMap_[virtualNodeNum].End() );
          } // loop over face nodes
          faces.insert(faceNum);
        }
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


    // ===================
    //  1) SBM-Definition
    // ===================
    EntityNodesType::const_iterator entIt;

    // loop over all edges -> block #0 (if gradients are disabled)
    const EntityNodesType& eNodes = vNodesCont_[BaseFE::EDGE];
    for(entIt = eNodes.begin(); entIt != eNodes.end(); ++entIt ) {
      UInt numNodes = entIt->second.GetSize();
      for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
        UInt virtNode = entIt->second[iNode];
        sbmBlocks[0][fctId].insert(nodeMap_[virtNode].Begin(),
                                   nodeMap_[virtNode].End() );
      }
    }

    // loop over all faces -> block #1
    const EntityNodesType& fNodes = vNodesCont_[BaseFE::FACE];
    for(entIt = fNodes.begin(); entIt != fNodes.end(); ++entIt ) {
      std::set<Integer>  faceEqns;
      UInt faceNum = entIt->first;

      UInt numNodes = entIt->second.GetSize();
      for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
        UInt virtNode = entIt->second[iNode];
        sbmBlocks[1][fctId].insert(nodeMap_[virtNode].Begin(),
                                   nodeMap_[virtNode].End() );
        faceEqns.insert(nodeMap_[virtNode].Begin(),
                        nodeMap_[virtNode].End() );
      } // loop: faceNodes

      // Make sure to use each face in only one minor block
      if( faces.find(faceNum) == faces.end() ) {
        if( faceEqns.size())
          minorBlocks[1].Push_back(faceEqns);
      } // check: face mapped
    } // loop: faces 

    
    // collect all inner nodes -> block #2
    const EntityNodesType& iNodes = vNodesCont_[BaseFE::INTERIOR];
    for(entIt = iNodes.begin(); entIt != iNodes.end(); ++entIt ) {
      std::set<Integer>  innerEqns;
      UInt numNodes = entIt->second.GetSize();
      for( UInt iNode = 0; iNode < numNodes; ++iNode ) {
        UInt virtNode = entIt->second[iNode];
        sbmBlocks[2][fctId].insert(nodeMap_[virtNode].Begin(),
                                   nodeMap_[virtNode].End() );
        innerEqns.insert(nodeMap_[virtNode].Begin(), nodeMap_[virtNode].End());
      }

      if( innerEqns.size())
        minorBlocks[2].Push_back(innerEqns);
    } 
  }

  void FeSpaceHCurlHi::GetThinFaceGroups( StdVector<StdVector<UInt> >&  fg ) {


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
          for( UInt iFace = 0; iFace < ptEl->extended->faces.GetSize(); ++iFace ) {

            if( shape.faceLocDirs[iFace][0] != minDir &&
                shape.faceLocDirs[iFace][1] != minDir ) {
              //            if( shape.faceLocDirs[iFace][0] == minDir ||
              //                shape.faceLocDirs[iFace][1] == minDir ) {
              elemFaces.insert( ptEl->extended->faces[iFace]);
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

  bool FeSpaceHCurlHi::IsSameEntityApproximation( shared_ptr<EntityList> list,
                                                  shared_ptr<FeSpace> space ) {
    if( this->GetSpaceType()  != space->GetSpaceType()  ) {
      return false;
    }
    if( this->IsHierarchical() != space->IsHierarchical()) {
      return false;
    }

    // Cast other space to same type
    shared_ptr<FeSpaceHCurlHi> otherSpace = dynamic_pointer_cast<FeSpaceHCurlHi>(space);

    EntityList::ListType actListType = list->GetType();
    if ( ! (actListType == EntityList::ELEM_LIST) &&
        ! (actListType == EntityList::SURF_ELEM_LIST) &&
        ! (actListType == EntityList::NC_ELEM_LIST))  {
      return true;
    }

    // Loop over all elements
    EntityIterator it = list->GetIterator();
    for( it.Begin(); !it.IsEnd(); it++) {
      FeHCurlHi * myElem = static_cast<FeHCurlHi*>(this->GetFe(it));
      FeHCurlHi * otherElem = static_cast<FeHCurlHi*>(otherSpace->GetFe(it));
      if( !( *myElem == *otherElem) ) {
        return false;
      } else {
      }
    }
    return true;
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

    // Use .at() for thread-safe read-only access (throws if key missing)
    const auto& refElemsRegion = refElems_.at(region);
    if(refElemsRegion.find(type) == refElemsRegion.end()){
      EXCEPTION("fespaceh1::getfe( const entityiterator): requested fetype which is noch supported by space");
    }

#ifdef USE_OPENMP
    if(isFinalized_ && omp_get_num_threads()>1)
      ret = TL_refElems_.at(region).Mine().at(type);
    else
      ret = refElems_.at(region).at(type);
#else
    ret = refElems_.at(region).at(type);
#endif
    return ret;
  }

} // end of namespace
