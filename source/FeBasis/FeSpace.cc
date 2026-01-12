#include "FeSpace.hh"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/algsys/SolStrategy.hh"

// include headers of sub-classes for factory method
#include "H1/FeSpaceH1Nodal.hh"
#include "L2/FeSpaceL2Nodal.hh"
#include "H1/FeSpaceH1Hi.hh"
#include "HCurl/FeSpaceHCurlHi.hh"
#include "FeSpaceConst.hh"


#ifdef _OPENMP
#include <omp.h>
#endif

namespace CoupledField {


// ==========================================================================
//   A P P R O X O R D E R    C L A S S  
// ==========================================================================

ApproxOrder::ApproxOrder( ) {
  isoOrder_ = 0;
  maxOrder_= 0;
  dim_ = 0;
  isIsotropic_ = false;
  polyMapType_ = UNDEF;
}

ApproxOrder::ApproxOrder(UInt dim ) {
  isoOrder_ = 0;
  maxOrder_= 0;
  dim_ = dim;
  isIsotropic_ = false;
  polyMapType_ = UNDEF;
}
  
  void ApproxOrder::SetIsoOrder( UInt isoOrder ) {
    isoOrder_ = isoOrder;
    isIsotropic_ = true;
  }
  
  void ApproxOrder::SetAnisoOrder( const Matrix<UInt>& anisoOrder ) {
    
    // the structure of the anisoOrder matrix has the following 
    // structure:
    
    //         dofs ->
    //      ux  uy uz
    //   xi
    //  eta
    // zeta
    
    assert( anisoOrder.GetNumRows() == dim_ );
    
    anisoOrder_ = anisoOrder;
    isIsotropic_ = false;
    maxOrder_ = 0;
    // determine maximum w.r.t. to local coordinate directions
    maxOrderLocDir_.Resize( anisoOrder_.GetNumRows() );
    maxOrderLocDir_.Init( 0 );
    for( UInt iLoc = 0; iLoc < anisoOrder_.GetNumRows(); iLoc++ ) {
      for( UInt iDof = 0; iDof < anisoOrder_.GetNumCols(); iDof++ ) {
        if( anisoOrder_[iLoc][iDof] > maxOrderLocDir_[iLoc] ) {
          maxOrderLocDir_[iLoc] = anisoOrder_[iLoc][iDof];
        }
        if( anisoOrder_[iLoc][iDof] > maxOrder_ ) {
          maxOrder_ = anisoOrder_[iLoc][iDof];
        }
      }
    } 
  }
  
  void ApproxOrder::SetPolyMapping(PolyMapType type) {
    polyMapType_ = type;
  }
  
  bool ApproxOrder::IsIsotropic() const {
    return isIsotropic_;
  }
  
  ApproxOrder::PolyMapType ApproxOrder::GetPolyMapType() const {
    return polyMapType_;
  }
  
  UInt ApproxOrder::GetIsoOrder() const {
    if( !isIsotropic_) {
      WARN(("Approximation is not isotropic!"))
    }
    return isoOrder_;
    
  }
  
  void ApproxOrder::GetAnisoOrder( Matrix<UInt>& order ) const {
    order = anisoOrder_;
  }
  
  UInt ApproxOrder::GetMaxOrder() const {
    return maxOrder_;
  }
  
  const StdVector<UInt>& ApproxOrder::GetMaxOrderLocDir( ) const {
    return maxOrderLocDir_;
  
  }
  
  StdVector<UInt> ApproxOrder::GetDofOrder( UInt dof ) const {
    StdVector<UInt> order;
    if( isIsotropic_) {
      assert(dim_ != 0 );
      order.Resize( dim_ );
      order.Init( isoOrder_ );
    } else {
      assert(dof < anisoOrder_.GetNumCols() );
      order.Resize( anisoOrder_.GetNumRows() );
      for( UInt i = 0; i < order.GetSize(); ++i ) {
        order[i] = anisoOrder_[i][dof];
      }
    }
    return order;
  }
  
  std::string ApproxOrder::ToString() const {
    if( isIsotropic_) {
      return lexical_cast<std::string>(isoOrder_);
    } else {
      return anisoOrder_.ToString();
    }
  }


// ==========================================================================
//   M A I N   F E S P A C E    C L A S S
// ==========================================================================

  // declare class specific logging stream
  DEFINE_LOG(feSpace, "feSpace")
  
  
  //! Constructor
  FeSpace::FeSpace( PtrParamNode paramNode, PtrParamNode infoNode, 
                    Grid* ptGrid ){
   
    lagrangeSurfSpace_ = false;
    myParam_ = paramNode;
    infoNode_ = infoNode;
    ptGrid_ = ptGrid;
    type_ = UNDEF_SPACE;
    mapType_ = GRID;
    polyType_ = UNDEF_POLY;
    isFinalized_ = false;
    isContinuous_ = true;
    isHierarchical_ = false;
    numEqns_ = 0;
    numFreeEquations_ = 0;
    solStep_ = 1;
    orderOffset_ = 0;

    bcCounter_[NOBC] = 0;
    bcCounter_[HDBC] = 0;
    bcCounter_[IDBC] = 0;
    bcCounter_[CS] = 0;
    
    lastElemNum_.Set(0);

    // get already integrationScheme from grid
    // In the future we could also create our own instance,
    // where only the maximum integration order defined by 
    // the user gets initialized
    intScheme_ = ptGrid_->GetIntegrationScheme();
  }

  FeSpace::~FeSpace(){

  }

  shared_ptr<FeSpace> FeSpace::CreateInstance( PtrParamNode aNode, 
                                               PtrParamNode infoNode,
                                               SpaceType reqType,
                                               Grid* ptGrid,
                                               bool isAVExc) {
     

    /* One big Problem> Due to the splitting of spaces in Hi and Lo/Lagrange, it is not
     * possible to deal with different Polynomial types in different regions
     * e.g. acoustic pressure is approximated by 2nd order Spectral elements in region1 and
     * anisotropic Legendre elements in region2 coupled by NcInterfaces.
     * This is to be changed in the future!
     * Already possible:
     *    - The combination of Grid Mapping and Lagrange higher order is ok.
     *    - The combination of Lagrange in different orders is ok
     *    - The combination of Legendre in different anisotropic orders is ok
     *    - Differnt Integration types for different regions is not affected
     so we follow the idea:
      1. obtain the fePolynomialList from ParamNode
      2. loop over the list and check the following things
         a. if the useGrid Attribute is set, continue
         b. else we store the given Polytype (e.g. Lagrange)
         c. if we find later on a incompatible Polytype we throw an error
      3. Create the Correct space and pass the parameter node to the constructor to
         determine the correct reference elements.
    */

    LOG_TRACE(feSpace) << "Creating FeSpace instance";
    LOG_DBG(feSpace) << "\trequired type: " << SpaceTypeEnum.ToString(reqType);
    
    //= PolyTypeEnum.Parse(polyStr);
    PolyType polyType = UNDEF_POLY;
    //obtain the fePolynomialList
    PtrParamNode polyNode = aNode->GetRoot()->Get("fePolynomialList", ParamNode::PASS );
    if(!polyNode){
      polyType = LAGRANGE;
      LOG_DBG(feSpace) << "No explicit definition available, using default";
    }else{
      
      
      // loop over all polynomial definitions, get the id and query, if the
      // polynomial gets also referenced in the <reegionList> of this PDE.
      PtrParamNode regList = aNode->Get("regionList", ParamNode::PASS);
      if( !regList)
        regList = aNode->Get("surfRegionList", ParamNode::PASS);
      if( !regList) {
        EXCEPTION("Could not find <regionList> or <surfaceRegionList> element");
      }
      ParamNodeList polys = polyNode->GetChildren();
      std::set<shared_ptr<ParamNode> > usedPolys;
      std::set<std::string> usedIds; // check for multiple definitions of Ids
      for( UInt i = 0; i < polys.GetSize(); ++i ) {
        std::string id = polys[i]->Get("id")->As<std::string>();
        // check, if id was already used
        if( usedIds.find( id ) != usedIds.end() ) {
          EXCEPTION("Polynomial ID '" << id << "' defined multiple times" );
        } else {
          usedIds.insert( id );
        }
        if( regList->HasByVal("region", "polyId", id ) ) {
          usedPolys.insert(polys[i]);
        }
        if( regList->HasByVal("surfRegion", "polyId", id ) ) {
          usedPolys.insert(polys[i]);
        }
        // check for "velPolyId" and "presPolyId" in the PDE definition (used in LinFlowPDE)
        std::vector<std::string> polyIds = {"velPolyId", "presPolyId"};
        for (const std::string & polyId :polyIds)
        {
          if ( aNode->HasByVal(polyId,id) ) {
            usedPolys.insert(polys[i]);
          }
        }
      }
      
      // check consistency of referenced polynomials
      // a) if no polynomial gets referenced, use the default lagrange one
      if( usedPolys.size() == 0 ) {
        polyType = LAGRANGE;
      }
      
      // b) if one polynomial gets referenced, use this one
      if( usedPolys.size() == 1 ) {
        polyType = PolyTypeEnum.Parse((*usedPolys.begin())->GetName());
      }
      
      // c) if several ones get referenced, check for same type
      if( usedPolys.size() > 1 ) {
        polyType = PolyTypeEnum.Parse((*usedPolys.begin())->GetName());
        std::set<shared_ptr<ParamNode> >::const_iterator it = usedPolys.begin();
        it++;
        while (it != usedPolys.end()) {
          PolyType actPoly = PolyTypeEnum.Parse((*it)->GetName());
          if( actPoly != polyType ) {
            EXCEPTION( "Can not mix '" << PolyTypeEnum.ToString(polyType)
                       << "' and '" << PolyTypeEnum.ToString(actPoly) 
                       << "' spaces in one PDE.");
          }
          ++it;
        }
      }
      
      
      
//      ParamNodeList lagList  = polyNode->GetList("Lagrange");
//      ParamNodeList legList  = polyNode->GetList("Legendre");
//      if(lagList.GetSize() > 0 && legList.GetSize() > 0){
//        EXCEPTION("Different polynomial types for different regions not supported!")
//      }else if(lagList.GetSize() > 0 ){
//        polyType = LAGRANGE;
//      }else if(legList.GetSize() > 0 ){
//        polyType = LEGENDRE;
//      }else{
//        EXCEPTION("An error occurred while reading the XML file.Specify at least one fePolynomial");
//      }
    }
    // switch depending on space type
    shared_ptr<FeSpace> ret;
    switch( reqType ) {
      case H1:
        LOG_DBG(feSpace) << "Creating H1 space";
        if( polyType == LAGRANGE )  {
          ret.reset(new FeSpaceH1Nodal(aNode, infoNode, ptGrid));
        } else {  
          ret.reset(new FeSpaceH1Hi(aNode, infoNode, ptGrid));
        }
        break;
      case HCURL:
        LOG_DBG(feSpace) << "Creating HCurl space";
        //if( polyType == LAGRANGE )  {
        //  // create explicit lower order HCurl space
        //  EXCEPTION("Explicit lower order space H-Curl not defined yet");
        //} else {
        ret.reset(new FeSpaceHCurlHi(aNode, infoNode, ptGrid));
        //}
        break;
      case CONSTANT:
        LOG_DBG(feSpace) << "Creating CONST space";
        ret.reset(new FeSpaceConst(aNode, infoNode, ptGrid, isAVExc));
        break;
      case HDIV:
      case L2:
        LOG_DBG(feSpace) << "Creating L2 space";
        if( polyType == LAGRANGE )  {
          ret.reset(new FeSpaceL2Nodal(aNode, infoNode, ptGrid));
        } else {
          EXCEPTION("Higher order L2 space not defined yet");
        }
        break;
        break;
      case UNDEF_SPACE:
        EXCEPTION("Can not create unknown function space");
        break;
    }
    return ret;    
  }


  void FeSpace::GetNodesOfEntities( StdVector<UInt>& nodes,
                                    shared_ptr<EntityList> ent,
                                    BaseFE::EntityType entType ) {

    // Case 1: No "virtual" nodes are present. Therefore we can directly use
    //          the grid-nodes of the geometric element
    if(mapType_ == GRID){
      // get name of entitylist
      StdVector<UInt> tempNodes;
      std::string name= ent->GetName();
      shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
      assert(feFct);
      
      // If we have a nodal entity list, we can directly access those nodes.
      // Otherwise the grid class is queried
      if( ent->GetType() == EntityList::NODE_LIST ) {
        tempNodes = static_cast<NodeList&>(*ent).GetNodes();
      } else {
        feFct->GetGrid()->GetNodesByName( tempNodes, name );
      }
      
      
      // careful: not all nodes of the entity list are necessarily mapped for
      // this space. So only select those nodes, also contained in the nodesType_ array.
      nodes.Clear();
      UInt numNodes = tempNodes.GetSize();
      for( UInt i = 0; i < numNodes; ++i ) {
        if( nodesType_.find(tempNodes[i]) != nodesType_.end() ) {
          nodes.Push_back( tempNodes[i] );
        }
      }
      
      

      // Case 2: We rely on the the information of the virtualNodes_ array.
      //         which gets filled depending on the number of "unknown" nodes
      //         initially.
    } else {
      EntityList::ListType defType = ent->GetType();
      EntityIterator entIt = ent->GetIterator();
      nodes.Resize(0);
      // UInt nCount = 0;
      switch (defType){
        case EntityList::REGION_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<Elem*> eList;
            shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
            assert(feFct);
            feFct->GetGrid()->GetElems(eList,entIt.GetRegion());
            StdVector<UInt> eNodes;
            for (UInt i = 0; i < eList.GetSize(); i ++ ){

              this->GetNodesOfElement(eNodes,eList[i],entType);
              for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
                nodes.Push_back(eNodes[aNode]);
              }
            }
          }
        break;
        case EntityList::ELEM_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            this->GetNodesOfElement(eNodes,entIt.GetElem(),entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::NC_ELEM_LIST:
        case EntityList::SURF_ELEM_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            this->GetNodesOfElement(eNodes,  entIt.GetSurfElem() ,entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::NODE_LIST:
          // nCount = 0;
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            if(gridToVirtualNodes_.find(entIt.GetNode()) != gridToVirtualNodes_.end()){
              for(UInt i = 0;i < gridToVirtualNodes_[entIt.GetNode()].GetSize();i++){
                nodes.Push_back( gridToVirtualNodes_[entIt.GetNode()][i]);
              }
            }
          }
        break;
        default:
          EXCEPTION("FeSpace::GetNodesOfEntities: Called with a list which is not supported ");
          break;
      }
    }
  }

  void FeSpace::GetNodesOfElement( StdVector<UInt>& nodes,
                                   const Elem* ptElem,
                                   BaseFE::EntityType entType){
    nodes.Clear(true);

    // Use const references for thread-safe read-only access
    const EntityNodesType& vNodes = vNodesCont_[BaseFE::VERTEX];
    const EntityNodesType& eNodes = vNodesCont_[BaseFE::EDGE];
    const EntityNodesType& fNodes = vNodesCont_[BaseFE::FACE];
    const EntityNodesType& iNodes = vNodesCont_[BaseFE::INTERIOR];

    // Get hold of element
    BaseFE * ptFe = this->GetFe(ptElem->elemNum);

    // Collect vertex nodes
    {
      UInt numNodes = ptElem->connect.GetSize();
      if( entType == BaseFE::VERTEX || entType == BaseFE::ALL ) {
        for( UInt i = 0; i < numNodes; ++i ) {
          auto it = vNodes.find(ptElem->connect[i]);
          if( it != vNodes.end() ) {
            const StdVector<UInt>& vertexNodes = it->second;
            for( UInt j = 0; j < vertexNodes.GetSize(); ++j ) {
              nodes.Push_back(vertexNodes[j]);
            }
          }
        }
      }
    }

    // Collect edge nodes
    if(eNodes.size() >0){
      UInt numEdges = ptElem->extended->edges.GetSize();

      if( entType == BaseFE::EDGE || entType == BaseFE::ALL ) {
        // Check for permutation
        if( ptFe->NeedsNodalPermutation() ) {
          StdVector<UInt> perm;
          for( UInt i = 0; i < numEdges; ++i ) {
            auto it = eNodes.find(std::abs(ptElem->extended->edges[i]));
            if( it != eNodes.end() ) {
              const StdVector<UInt>& edgeNodes = it->second;
              ptFe->GetNodalPermutation( perm, ptElem, BaseFE::EDGE, i);
              for( UInt j = 0; j < edgeNodes.GetSize(); ++j ) {
                nodes.Push_back(edgeNodes[perm[j]]);
              }
            }
          }
        } else {
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
    }

    // Collect face nodes
    if(fNodes.size() >0){
      UInt numFaces = ptElem->extended->faces.GetSize();
      if( entType == BaseFE::FACE || entType == BaseFE::ALL ) {
        if( ptFe->NeedsNodalPermutation() ) {
          StdVector<UInt> perm;
          for( UInt i = 0; i < numFaces; ++i ) {
            auto it = fNodes.find(std::abs(ptElem->extended->faces[i]));
            if( it != fNodes.end() ) {
              ptFe->GetNodalPermutation( perm, ptElem, BaseFE::FACE, i);
              const StdVector<UInt>& faceNodes = it->second;
              for( UInt j = 0; j < faceNodes.GetSize(); ++j ) {
                nodes.Push_back(faceNodes[perm[j]]);
              }
            }
          }
        } else {
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

    // Ensure, that at least one virtual node is present
    if( nodes.GetSize() == 0 ) { 
      EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element #"
          << ptElem->elemNum << " of region " 
          <<  ptGrid_->GetRegion().ToString(ptElem->regionId));
    }
  }

  void FeSpace::GetIntegration( BaseFE * fe, 
                                RegionIdType region,
                                IntScheme::IntegMethod & method,
                                IntegOrder & order ){
    
    
    // find integration definition for this region 
    std::map<RegionIdType, IntegDefinition>::iterator it;
    it = regionIntegration_.find(region);
    if( it == regionIntegration_.end() ) {
      it = regionIntegration_.find(ALL_REGIONS);
    }
    IntegDefinition &id = it->second;
    LOG_DBG(feSpace) << "returning integScheme for region " << region;
    
    // We distinguish the following cases:
    // 1) order is given ABSOLUTE: 
    //    in this case we use directly this value
    // 2) order is given RELATIVE:
    //    distinguish if element is
    //    a) ISOTROPIC: get order p of element and set integration
    //                  order to 2*(p+1) + relativeOrder-1.
    //                  This leads e.g. to integration order 3
    //                  for p = 1 and order 5 for p = 2.
    //    b) ANISOTROPIC: get maximum order in each direction
    //                    and set integration oder to
    //                    2*( p_maxLicDir) + relativeOrder_locDir
    method = id.method;
    LOG_DBG2(feSpace) << "\tmethod: " 
        << IntScheme::IntegMethodEnum.ToString(method);
    if( id.mode == INTEG_MODE_ABSOLUTE ) {

      // ==== ABSOLUTE order given ===
      order = id.order;
      LOG_DBG2(feSpace) << "\tuser order (absolute):" << order.ToString();
    } else {
      // === RELATIVE order given ===

      LOG_DBG2(feSpace) << "\tuser order (relative): " << id.order.ToString();
      // Check: If method is GAUSS and orider is anisotropic,
      //        we take the anisotropic version

      IntegOrder p;
      if( method == IntScheme::GAUSS &&
          !fe->IsIsotropic() ) {
        // --- anisotropic version ---
        StdVector<UInt> anisoOrder;
        fe->GetAnisoOrder(anisoOrder);
        LOG_DBG2(feSpace) << "\tanisotropic element order:" << anisoOrder.ToString();
        p.SetAnisoOrder(anisoOrder);
      } else {
        // --- isotropic version ---
        LOG_DBG2(feSpace) << "\tisotropic element order:" << fe->GetMaxOrder();
        p.SetIsoOrder(fe->GetMaxOrder());
      }
      // Note: in case of H-curl elements, we have to add an order of 1
      // as we use here also the polynomial order of 0 for lowest order Nedelec elements
      if( type_ == FeSpace::HCURL )
        p += 1;
      
      // -------------------------
      //  Final integration order
      // -------------------------
      // In Bilinearforms like the mass-matrix we have a product of two polynomials
      // of order p, so the integration order must be at least 2*p. However, this holds
      // true exactly only on the reference element. The Jacobian-mapping introduces
      // and additional influence (e.g. for highly distorted elements). Thus we
      // increase the order to 2*(p+1). In the end we decrease the order by -1
      // to get efficient integration rules. mostly compatible with values in literature.
      // This yields the following mapping
      // p = 1 => IntOrder: 3 => Gauss 1D: 2 points
      // p = 2 => IntOrder: 5 => Gauss 1D: 3 points
      // ....
      order = id.order;
      order += (p+1)*2 - 1; 
      LOG_DBG2(feSpace) << "\torder (final): " << order.ToString();
    }
  }

  void FeSpace::CreateVirtualNodes(){
    
    
    // store all regions the space is defined on
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    regions_ = feFct->GetRegions();
    
    //the function checks for the special case if we got only one element
    //in the region map and if this mapping is of type GRID
    //if so, we call a specialized but efficient function
    // if not, we get to the general case
    if(mapType_ == GRID && isContinuous_){
      this->CreateGridNodes();
    }else{
      this->CreatePolynomialNodes();
      //we have the general case so we iterate over all elements
    }
  }

  void FeSpace::CreateGridNodes(){

    LOG_TRACE(feSpace) << "starting to create virtual nodes based on GRID";

    StdVector< shared_ptr<EntityList> > fctEntList;
    StdVector<UInt> curNodes;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    fctEntList = feFct->GetEntityList();
    EntityNodesType& vNodes = vNodesCont_[BaseFE::VERTEX];

    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){
      LOG_DBG(feSpace) << "\tMapping entity list '" << fctEntList[actList]->GetName();
      EntityList::ListType actListType = fctEntList[actList]->GetType();
      if ( ! (actListType == EntityList::ELEM_LIST ||
          actListType == EntityList::SURF_ELEM_LIST ||
          actListType == EntityList::NC_ELEM_LIST) )  {
        continue;
      }
      shared_ptr<EntityList> actElemList =  fctEntList[actList];

      std::string name= actElemList->GetName();
      feFct->GetGrid()->GetNodesByName( curNodes, name );
      for ( UInt aNode= 0; aNode < curNodes.GetSize(); aNode++ ) {
        const UInt nodeNum = curNodes[aNode];
        if(gridToVirtualNodes_.find(nodeNum) == gridToVirtualNodes_.end()){
          gridToVirtualNodes_[nodeNum].Push_back(nodeNum);
          vNodes[nodeNum] = nodeNum;
        }
      }
      
    } //loop: lists

    for(std::map<UInt,StdVector<UInt> >::const_iterator it = gridToVirtualNodes_.begin(); it != gridToVirtualNodes_.end(); ++it) {
      //here we can hardcode to zero as we assume continuous approximation
      nodesType_[it->second[0]] = BaseFE::VERTEX;
    } // loop: lists


  }



  void FeSpace::CreatePolynomialNodes(){
    //follow the following algorithm
    // - loop over every element
    //  - get vertices and if not already mapped assign new virtual nodes
    //  - get edges and if not already mapped assign new virtual nodes
    //  - get faces and if not already mapped assign new virtual nodes
    //  - assign interior node to the element
    //  - now fill the Virtual node map with the information according to element orientation
    // - finally delete all intermediate arrays
    
    StdVector< shared_ptr<EntityList> > fctEntList;
    EntityList::ListType actListType = EntityList::NO_LIST;

    LOG_TRACE(feSpace) << "starting to create virtual nodes based on POLYNOMIAL";
    
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    fctEntList = feFct->GetEntityList();

    // This is the counter for the virtual node number
    UInt offset =0;

    // loop over all entitylists (i.e. regions)
    for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){


      actListType = fctEntList[actList]->GetType();
      //determine mapping type and order of current entity list
      //if we do not find the name of the reigon in our map, we fall back to the default
      //BE CARFUL, if the user has specified something for the volume region but not for surface regions,
      // it could be that the fallback to default leads to errors!

      //check if we got what we expected

      if ( ! (actListType == EntityList::ELEM_LIST) &&
           ! (actListType == EntityList::SURF_ELEM_LIST) &&
           ! (actListType == EntityList::NC_ELEM_LIST))  {
        LOG_DBG(feSpace) << "\tLEAVING";
        continue;
      }
      LOG_DBG(feSpace) << "========================================";
      LOG_DBG(feSpace) << "treating entityList '" << fctEntList[actList]->GetName() << "'";
      LOG_DBG(feSpace) << "========================================";
      //cast down to element list
      EntityList* actElemList = fctEntList[actList].get();

      // Get iterator of current element list
      EntityIterator entIt = actElemList->GetIterator();
      
      // ------------------------
      //  loop over all elements
      // ------------------------
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        // Fetch current finite element. This is performed by the specialized 
        // version, so this element "knows" already about its order, unknowns etc.
        BaseFE* ptFe = GetFe( entIt );
        const Elem* actEl = entIt.GetElem();  

        LOG_DBG2(feSpace) << "treating element #" << entIt.GetElem()->elemNum;
        LOG_DBG2(feSpace) << "mapped volume element #" << actEl->elemNum;

        StdVector<UInt> permutations; // initially size 0
        
        //===========================================================
        //Assign the BaseFE::VERTEX node numbers
        //===========================================================
        {
          LOG_DBG2(feSpace) << "mapping vertex nodes";
          UInt numVert = Elem::shapes[actEl->type].numVertices;
          EntityNodesType & vtn =  vNodesCont_[BaseFE::VERTEX];
          const StdVector<UInt> & elemNodes = actEl->connect;
          StdVector<UInt> numFncs;
          ptFe->GetNumFncs( numFncs, BaseFE::VERTEX );
          
          // loop over vertices
          for ( UInt iVert= 0; iVert< numVert; iVert++ ) {
            UInt vertexNum = elemNodes[iVert];
            StdVector<UInt> & en = vtn[vertexNum];
            if( numFncs.GetSize() > 0 ) {
              if ( en.GetSize() == 0 ) {
                for( UInt iFunc = 0; iFunc < numFncs[iVert]; ++iFunc ) {
                  en.Push_back(++offset);
                  nodesType_[offset] = BaseFE::VERTEX;
                  LOG_DBG3(feSpace) << "adding " << offset << " to node_";
                }
              }

              // Assign vertex also to virtualToGridNodes
              if(gridToVirtualNodes_.find(vertexNum) == gridToVirtualNodes_.end()){
                LOG_DBG3(feSpace) << "gridToVirtualNodes[" << vertexNum << "] = " << offset;
                gridToVirtualNodes_[vertexNum].Push_back(offset);
              } else {
                LOG_DBG3(feSpace) << "vertex " << vertexNum << " already mapped to virtualNode " << 
                    gridToVirtualNodes_[vertexNum];
              }
            }
          }  // loop: vertices
        }// mapping of vertex nodes
        feFct->GetGrid()->MapEdges();
        feFct->GetGrid()->MapFaces();

        ElemShape actShape = Elem::shapes[actEl->type];

        //===========================================================
        //Assign the Edge node numbers
        //===========================================================
        {
          LOG_DBG2(feSpace) << "mapping edge nodes";
          EntityNodesType & etn =  vNodesCont_[BaseFE::EDGE];
          StdVector<UInt> numFncs;
          ptFe->GetNumFncs( numFncs, BaseFE::EDGE );
          // loop over edges
          for ( UInt iEdge=0; iEdge < actShape.numEdges; iEdge++) {
            UInt edgeNum = std::abs(actEl->extended->edges[iEdge]);
            StdVector<UInt> & en = etn[edgeNum];

            // check if edge got already numbered
            if( numFncs.GetSize() > 0 ) {
              if ( en.GetSize() == 0 ) {
                for( UInt iFunc = 0; iFunc < numFncs[iEdge]; ++iFunc ) {

                  en.Push_back(++offset);
                  nodesType_[offset] = BaseFE::EDGE;
                  
                  LOG_DBG3(feSpace) << "adding " << offset << " to node_";
                } // loop: fncs
              } else {
                if( en.GetSize() != numFncs[iEdge] ) {
                  EXCEPTION("Edge " << iEdge+1 << ", edge #" << edgeNum << 
                            " of element #" << actEl->elemNum << " got already "
                            << en.GetSize() << " virtual nodes and is now assigned "
                            << numFncs[iEdge] << " virtual nodes" );
                }
              }
            }
          }  // loop: edges
        } // mapping of edge nodes

        //===========================================================
        //Assign the Face node numbers
        //===========================================================
        {
          LOG_DBG2(feSpace) << "mapping face nodes";
          EntityNodesType & ftn =  vNodesCont_[BaseFE::FACE];
          StdVector<UInt> numFncs;
          ptFe->GetNumFncs( numFncs, BaseFE::FACE );

          // loop over faces
          for ( UInt iFace=0; iFace < actShape.numFaces; iFace++) {
            UInt faceNum = actEl->extended->faces[iFace];
            StdVector<UInt> & fn = ftn[faceNum];
            // check if face got already numbered
            if( numFncs.GetSize() > 0 ) {
              if ( fn.GetSize() == 0 ) {
                for( UInt iFunc = 0; iFunc < numFncs[iFace]; ++iFunc ) {
                  fn.Push_back(++offset);
                  nodesType_[offset] = BaseFE::FACE;
                  LOG_DBG3(feSpace) << "adding " << offset << " to node_";
                } // loop: fncs
              } else {
                if( fn.GetSize() != numFncs[iFace] ) {
                  EXCEPTION("Face " << iFace+1 << 
                            " of element #" << actEl->elemNum << " got already "
                            << fn.GetSize() << " virtual nodes and is now assigned "
                            << numFncs[iFace] << " virtual nodes" );
                }
              }
            }
          } // loop: faces
        } // mapping of face nodes

        //===========================================================
        //Assign the Interior node numbers
        //===========================================================
        {
          LOG_DBG2(feSpace) << "mapping interior nodes";

        
          StdVector<UInt> numFncs;
          ptFe->GetNumFncs( numFncs, BaseFE::INTERIOR );

          // check, if element provides any interior functions at all
          if( numFncs.GetSize() > 0 ) {
            StdVector<UInt> & itn =  vNodesCont_[BaseFE::INTERIOR][actEl->elemNum];
            if( itn.GetSize() == 0 ) {
              for( UInt iFunc = 0; iFunc < numFncs[0]; ++iFunc ) {
                itn.Push_back(++offset);
                nodesType_[offset] = BaseFE::INTERIOR;
                LOG_DBG3(feSpace) << "adding " << offset << " to node_";
              } // loop: fncs
            } else {
              EXCEPTION("Interior of element #" << actEl->elemNum << " got already "
                        << itn.GetSize() << " virtual nodes and is now assigned "
                        << numFncs[0] << " virtual nodes" );
            }
          }
        } // mapping of interior nodes
      
      } // loop: elements
    } // loop: entity lists
    
    

    // --------------------------------
    // Print information to info block  
    // --------------------------------
    // Loop over entity types
    PtrParamNode vNodeInfo = infoNode_->Get("virtualNodes");
    
    boost::unordered_map<BaseFE::EntityType, EntityNodesType>::const_iterator typeIt;
    typeIt = vNodesCont_.begin();
    UInt numVirtNodes = 0;
    for( ; typeIt != vNodesCont_.end(); ++typeIt ) {
      const BaseFE::EntityType & type = typeIt->first;
      const EntityNodesType& nodesVec = typeIt->second;
      
      // Collect number of entities and total nodes
      UInt numEntities = nodesVec.size();
      UInt vNodesOfType = 0;
      EntityNodesType::const_iterator nodeIt = nodesVec.begin();
      for( ; nodeIt != nodesVec.end(); ++nodeIt ) {
        vNodesOfType += nodeIt->second.GetSize();
      }
      numVirtNodes += vNodesOfType;
      std::string entType = BaseFE::entityType.ToString(type);
      vNodeInfo->Get(entType)->Get("numEntitiies")->SetValue(numEntities);
      vNodeInfo->Get(entType)->Get("numVirtualNodes")->SetValue(vNodesOfType);
    }
    vNodeInfo->Get("numVirtualNodes")->SetValue(numVirtNodes);
    
    
    // trim all vectors to get unused memory back
    boost::unordered_map<BaseFE::EntityType, EntityNodesType>::iterator it;
    it = vNodesCont_.begin();
    for( ; it != vNodesCont_.end(); ++it ) {
      EntityNodesType & evn = it->second;
      EntityNodesType::iterator evnIt = evn.begin();
      for( ; evnIt != evn.end(); evnIt++ ) {
        evnIt->second.Trim();
      }
    }
    LOG_TRACE(feSpace) << "finished creation of virtual nodes";
  }
  
  // ************************************************************************
  // GENERATE REGION SPECIFIC DATA AND PROCESS USER INPUT
  // ************************************************************************
  void FeSpace::SetRegionApproximation(RegionIdType region, std::string polyId, std::string integId){
    // add region to _regions
    regions_.insert(region);

    std::string regionName = ptGrid_->GetRegion().ToString(region);
    PtrParamNode regionNode = infoNode_->Get("regionList")->Get(regionName);
    
    RegionIdType pReg,iReg;

    //SECTION1: Polynomials
    if(polyNodes_.find(polyId)!=polyNodes_.end()){
      //the user specified a valid polyId so we read its data and call create ref elems
      //additionally we let the space read specific attributes only valid for himself
      PtrParamNode pNode = polyNodes_[polyId];
      ReadCustomAttributes(pNode,region);
      ApproxOrder order(ptGrid_->GetDim() );
      MappingType curMap;
      ReadPolyNode(pNode,curMap,order);
      SetRegionElements(region,curMap,order,regionNode);
      pReg = region;
    }else if(polyId=="default"){
      //the user requested the default but did not specify it so we set it here
      SetDefaultElements(infoNode_->Get("regionList")->Get("default"));
      pReg = ALL_REGIONS;
    }else{
      EXCEPTION("The polynomial id does not match any in the fePolynomialList: " << polyId);
    }
    //Section2: Integral Definitions
    if(integNodes_.find(integId)!=integNodes_.end()){
      //the user specified a valid polyId so we read its data and set the region integration
      PtrParamNode iNode = integNodes_[integId];
      IntegOrder order;
      IntegOrderMode curMode;
      IntScheme::IntegMethod iMeth;
      ReadIntegNode(iNode,iMeth,order,curMode);
      SetRegionIntegration(region,iMeth,order,curMode, regionNode);
      iReg = region;
    }else if(integId=="default"){
      //the user requested the default but did not specify it so we set it here
      SetDefaultIntegration(infoNode_->Get("regionList")->Get("default"));
      iReg = ALL_REGIONS;
    }else{
      EXCEPTION("The integration id does not match any in the IntegrationSchemeList: " << integId);
    }
    polyToIntegMap[pReg].insert(iReg);
  }

  void FeSpace::SetDefaultRegionApproximation(){
    RegionIdType pReg=ALL_REGIONS,iReg=ALL_REGIONS;
    SetDefaultIntegration(infoNode_->Get("regionList")->Get("default"));
    SetDefaultElements(infoNode_->Get("regionList")->Get("default"));
    polyToIntegMap[pReg].insert(iReg);
  }
  
  const Elem* FeSpace::GetVolElem( const Elem* ptElem ) const {

    const Elem * ret = nullptr;
    
    // check dimension
    UInt elemDim = Elem::shapes[ptElem->type].dim;
    UInt gridDim = ptGrid_->GetDim(); 
    if( gridDim ==  elemDim ) {
      // 1) Volume case: simply return volume element
      ret = ptElem;
    } else if( gridDim-elemDim == 1 ) {
      // 2) Real surface case

      if(lagrangeSurfSpace_) 
      {
        return ptElem;
      }

      const SurfElem * ptSurfEl = dynamic_cast<const SurfElem*>(ptElem);
      boost::array<Elem*,2>::const_iterator it = ptSurfEl->ptVolElems.begin();
      for( ; it != ptSurfEl->ptVolElems.end(); ++it ) {
        // check, if element is set at all
        if( *it) {
          if(regions_.find( (*it)->regionId) != regions_.end()) {
            ret = *it; 
            break;
          }
        }
      } // loop over volume element neighbors

      // check, if element could be found
      if( !ret) {
        EXCEPTION("Could not find a suitable volume neighbor for surface element #"
            << ptSurfEl->elemNum << " in region "
            << ptGrid_->GetRegion().ToString( ptSurfEl->regionId ) << "." );
      }
    } else {
      // 3) 1D element in 3D simulation
      
      // In case the (only) edge is not numbered yet, we perform a 
      // mapping of edges now
      if(ptElem->extended->edges.GetSize() != 1) {
        ptGrid_->MapEdges();
      }
      UInt edgeNr = std::abs(ptElem->extended->edges[0]);
      const StdVector<Elem*> & neighbors = ptGrid_->GetEdge(edgeNr).neighbors; 
     
      // loop until element with highest dimension and acceptable regionId 
      // is found
      UInt ind = 0;
      bool found = false;
      while( ind < neighbors.GetSize() ) {
        if (Elem::shapes[neighbors[ind]->type].dim == gridDim  &&
            regions_.find(neighbors[ind]->regionId) != regions_.end() )  {
          found = true;
          break;
        }
        ind++;
      }
      if( !found )  {
        EXCEPTION("Could not find a suitable volume neighbor for  element #"
                   << ptElem->elemNum << ". " );
      }
      ret = neighbors[ind];
    }
    
    return ret;
  }

  void FeSpace::SetRegionIntegration(RegionIdType region, IntScheme::IntegMethod method ,
                                     const IntegOrder& order,
                                     IntegOrderMode mode, PtrParamNode infoNode ){
    regionIntegration_[region].method = method;
    regionIntegration_[region].order = order;
    regionIntegration_[region].mode = mode;
  }

  void FeSpace::ReadIntegNode(PtrParamNode node,  IntScheme::IntegMethod & iMeth,
                              IntegOrder& order, IntegOrderMode & mode){
    std::string methodStr = node->Get("method")->As<std::string>();
    iMeth = IntScheme::IntegMethodEnum.Parse(methodStr,IntScheme::UNDEFINED);
    if(iMeth == IntScheme::UNDEFINED){
      EXCEPTION("got undefined integration. This will lead to errors! Input string: " << methodStr);
    }
    std::string modeStr = node->Get("mode")->As<std::string>();
    mode = IntegOrderModeEnum.Parse(modeStr,INTEG_MODE_ABSOLUTE);

    //only isotropic order supported right now
    UInt isoOrder = node->Get("order")->As<UInt>();
    //create order Matrix, we only support isotropic interation right now
    order.SetIsoOrder( isoOrder );
  }

  void FeSpace::ReadPolyNode(PtrParamNode node, MappingType & mapType, ApproxOrder & order){

    // query for one of the following subchildren:
    //   - gridOrder
    //   - isoOrder
    //   - anisoOrder

    if(node->Has("gridOrder") ) {
      mapType= GRID;
      UInt gridOrder = ptGrid_->IsQuadratic() ? 2 : 1;
      order.SetIsoOrder( gridOrder );
      order.SetPolyMapping(ApproxOrder::GRID_TYPE);

    } else if(node->Has("isoOrder")) {
      PtrParamNode isoOrderNode = node->Get("isoOrder", ParamNode::PASS );
      mapType= POLYNOMIAL;
      Integer isoOrder = isoOrderNode->As<Integer>();
      order.SetIsoOrder(isoOrder);
      
      // Try to find out polynomial mapping type for lagrange type
      ApproxOrder::PolyMapType polyMapType = ApproxOrder::TENSOR_TYPE;
      if (polyType_ == LAGRANGE ) {
        if(node->Get("isoOrder")->Get("serendipity")->As<bool>() == true) 
          polyMapType = ApproxOrder::SERENDIPITY_TYPE;
      }
      order.SetPolyMapping(polyMapType);
     

    } else if(node->Has("anIsoOrder")) {
      PtrParamNode anIsoOrderNode = node->Get("anIsoOrder", ParamNode::PASS );
      mapType= POLYNOMIAL;

      UInt dim = ptGrid_->GetDim();
      Matrix<UInt> orderMat;
      orderMat.Resize(dim,3);
      StdVector<std::string> dofs(3);
      anIsoOrderNode->GetValue("dof1",dofs[0],ParamNode::PASS);
      anIsoOrderNode->GetValue("dof2",dofs[1],ParamNode::PASS);
      anIsoOrderNode->GetValue("dof3",dofs[2],ParamNode::PASS);
      if( dofs[2] == "" ) {
        orderMat.Resize(dim,2); 
      }
      if( dofs[1] == "" ) {
        orderMat.Resize(dim,1);
      }
      char_separator<char> sep(" ");

      for(UInt i = 0;i < orderMat.GetNumCols();i++){
        boost::tokenizer<char_separator<char> > tokens(dofs[i],sep);
        boost::tokenizer<char_separator<char> >::iterator tokIt=tokens.begin();
        for(UInt j = 0;j < orderMat.GetNumRows();j++){
          orderMat[j][i] = lexical_cast<UInt>(*tokIt);
          tokIt++;
        }
      }
      order.SetAnisoOrder(orderMat);
      // Anisotropic order makes only sense for tensorial functions.
      order.SetPolyMapping(ApproxOrder::TENSOR_TYPE);
      
    } else {
      EXCEPTION( "Could not find definition of polynomial order" );
    }
  }


  void FeSpace::ReadIntegList(){
    //obtain the fePolynomialList
    PtrParamNode integNode = myParam_->GetRoot()->Get("integrationSchemeList", ParamNode::PASS );
    if(integNode){
      ParamNodeList iList = integNode->GetList("scheme");
      for(UInt aInt = 0; aInt < iList.GetSize();aInt++){
        std::string curId = iList[aInt]->Get("id")->As<std::string>();
        if(integNodes_.find(curId) != integNodes_.end()){
          EXCEPTION("Found multiple IntegrationSchemes with the same id in XML.\n" << \
                    "Im confused and going to exit...");
        }
        integNodes_[curId] = iList[aInt];
      }
    }
  }

  void FeSpace::ReadPolyList(){
        PtrParamNode polyNode = myParam_->GetRoot()->Get("fePolynomialList", ParamNode::PASS );
        if(polyNode){
          std::string polyName = PolyTypeEnum.ToString(polyType_);
          ParamNodeList pList = polyNode->GetList(polyName);
          for(UInt aPol = 0; aPol < pList.GetSize();aPol++){
            std::string curId =  pList[aPol]->Get("id")->As<std::string>();
            if(polyNodes_.find(curId) != polyNodes_.end()){
              EXCEPTION("Found multiple fePolynomials with the same id in XML.\n" << \
                        "Im confused and going to exit...");
            }
            polyNodes_[curId] = pList[aPol];
          }
        }
      }

  
  // ========================================================================
  //  EQUATION HANDLING
  // ========================================================================
  //@{ \name Equation Handling

  void FeSpace::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent ){

    //Get result for the feFunction
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();

    //Get "dimension" of one Unknown
    UInt dofsPerUnknown = GetNumDofs();

    //This if clause should be avoided if the functionality for higher order entities
    //is implemented
    //First cover the nodal Case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(dofsPerUnknown);
      if (gridToVirtualNodes_.find(node) != gridToVirtualNodes_.end()){
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[iDof] = nodeMap_[gridToVirtualNodes_[node][0]][iDof];
        }
      } else {
        // In case a node was not found, we reset the eqnarray to
        // size 0
        eqns.Resize(0);
      }

    } else if( ent.GetType() == EntityList::ELEM_LIST ||
        ent.GetType() == EntityList::SURF_ELEM_LIST ||
        ent.GetType() == EntityList::NC_ELEM_LIST) {

      // check if last accessed element is the same as this one

      UInt& lastElemNum             = lastElemNum_.Mine();
      StdVector<Integer>& lastEqns  = lastEqns_.Mine();
      StdVector<UInt>& eqnNodes     = eqnNodes_.Mine();

      UInt elemNum = ent.GetElem()->elemNum;
      if( elemNum == lastElemNum ) {
        eqns = lastEqns;
      } else {
        GetNodesOfElement(eqnNodes, ent.GetElem());
        UInt nrNodes = eqnNodes.GetSize();
        eqns.Resize( nrNodes * dofsPerUnknown );
        //eqns.Init();
        UInt pos = 0;
        for (UInt iNode = 0; iNode < nrNodes; iNode++ ) {
          const StdVector<Integer>& nodeEqns = nodeMap_[eqnNodes[iNode]];
          for (UInt iDof = 0; iDof < dofsPerUnknown; iDof++) {
            eqns[pos + iDof] = nodeEqns[iDof];
          }
          pos += dofsPerUnknown;
        }
        lastElemNum = elemNum;
        lastEqns = eqns;
      } // if: lastElemNum == elemNum
    } else {
      EXCEPTION("In FeSpace::GetEqns(): Supplied an iterator which is not supported by FeSpace");
    }
  }

  void FeSpace::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                         , UInt dof ){ 
    //Get result for the feFunction
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();

    //First cover the nodal/grid case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(1);
      eqns.Init();
      //without the if clause we would have a segfault in case of higher order approximaton
      //for quadratic meshes
      if (gridToVirtualNodes_.find(node)!= gridToVirtualNodes_.end()) {
        eqns[0] = nodeMap_[gridToVirtualNodes_[node][0]][dof];
      } else {
        // In case a node was not found, we reset the eqnarray to
        // size 0
        eqns.Resize(0);
      }
    } else if( ent.GetType() == EntityList::ELEM_LIST ||
        ent.GetType() == EntityList::SURF_ELEM_LIST||
        ent.GetType() == EntityList::NC_ELEM_LIST){
      StdVector<UInt> nodes;
      this->GetNodesOfElement(nodes,ent.GetElem());
      eqns.Resize( nodes.GetSize() );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqns[iNode] = nodeMap_[nodes[iNode]][dof];
      }
    } else {
      EXCEPTION("In FeSpace::GetEqns(StdVector,EntityIterator,UInt):  Supplied an iterator which is not supported by FeSpace");
    }
    
    }
  void FeSpace::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent
                           , UInt dof, BaseFE::EntityType entityType){ 
    //Get result for the feFunction
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();

    //First cover the nodal/grid case
    if ( ent.GetType() == EntityList::NODE_LIST ) {
      UInt node = ent.GetNode();
      eqns.Resize(1);
      eqns.Init();
      if(gridToVirtualNodes_.find(node)!= gridToVirtualNodes_.end()) {
        eqns[0] = nodeMap_[gridToVirtualNodes_[node][0]][dof];
      } else {
        // In case a node was not found, we reset the eqnarray to
        // size 0
        eqns.Resize(0);
      }
    } else if( ent.GetType() == EntityList::ELEM_LIST ||
        ent.GetType() == EntityList::SURF_ELEM_LIST||
        ent.GetType() == EntityList::NC_ELEM_LIST){
      StdVector<UInt> nodes;
      this->GetNodesOfElement(nodes,ent.GetElem(), entityType);
      eqns.Resize( nodes.GetSize() );
      eqns.Init();
      for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
        eqns[iNode] = nodeMap_[nodes[iNode]][dof];
      }
    } else {
      EXCEPTION("In FeSpace::GetEqns(StdVector,EntityIterator,UInt):  Supplied an iterator which is not supported by FeSpace");
    }

  }
  
  void FeSpace::GetEqns( StdVector<Integer>& eqns, const EntityIterator ent,
                         BaseFE::EntityType entityType) { 
     //Get result for the feFunction
     shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
     assert(feFct);
     shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();

     //Get "dimension" of one Unknown
     UInt dofsPerUnknown = GetNumDofs();
     
     //First cover the nodal/grid case
     if ( ent.GetType() == EntityList::NODE_LIST ) {
       UInt node = ent.GetNode();
       eqns.Resize(dofsPerUnknown);
       eqns.Init();
       if(gridToVirtualNodes_.find(node)!= gridToVirtualNodes_.end()) {
         for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
           eqns[iDof] = nodeMap_[gridToVirtualNodes_[node][0]][iDof];
         }
       } else {
         // In case a node was not found, we reset the eqnarray to
         // size 0
         eqns.Resize(0);
       }
     } else if( ent.GetType() == EntityList::ELEM_LIST ||
         ent.GetType() == EntityList::SURF_ELEM_LIST||
         ent.GetType() == EntityList::NC_ELEM_LIST){
       StdVector<UInt> nodes;
       this->GetNodesOfElement(nodes,ent.GetElem(), entityType);
       eqns.Resize( nodes.GetSize() * dofsPerUnknown);
       eqns.Init();
       for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
         for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
           eqns[(iNode*dofsPerUnknown) + iDof] = nodeMap_[nodes[iNode]][iDof];
         }
       }
     } else {
       EXCEPTION("In FeSpace::GetEqns(StdVector,EntityIterator,UInt):  Supplied an iterator which is not supported by FeSpace");
     }

  }

    void FeSpace::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem){
      //Get result for the feFunction
      shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
      assert(feFct);
      shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();

      //Get "dimension" of one Unknown
      UInt dofsPerUnknown = GetNumDofs();

      // Thread-local caching: check if we already have results for this element
      UInt& lastElemNum             = lastElemNum_.Mine();
      StdVector<Integer>& lastEqns  = lastEqns_.Mine();
      StdVector<UInt>& eqnNodes     = eqnNodes_.Mine();

      UInt elemNum = elem->elemNum;
      if( elemNum == lastElemNum && lastEqns.GetSize() > 0 ) {
        eqns = lastEqns;
        return;
      }

      // GetNodesOfElement now uses .at() for thread-safe read-only access
      this->GetNodesOfElement(eqnNodes,elem);
      UInt nrNodes = eqnNodes.GetSize();
      eqns.Resize( nrNodes * dofsPerUnknown );
      UInt pos = 0;
      for (UInt iNode = 0; iNode < nrNodes; iNode++ ) {
        const StdVector<Integer>& nodeEqns = nodeMap_[eqnNodes[iNode]];
        for(UInt iDof = 0; iDof < dofsPerUnknown; iDof++){
          eqns[pos + iDof] = nodeEqns[iDof];
        }
        pos += dofsPerUnknown;
      }

      // Store in thread-local cache
      lastElemNum = elemNum;
      lastEqns = eqns;
    }
    
    void FeSpace::GetEntityListEqns(StdVector<Integer>& eqns, 
                                    shared_ptr<EntityList> ent,
                                    UInt dof,
                                    BaseFE::EntityType entType) {
      std::set<Integer> allEqns;
      StdVector<Integer> tmp;
      EntityIterator it = ent->GetIterator();
      for( it.Begin(); ! it.IsEnd(); it++ ) {
        tmp.Init();
        GetEqns( tmp, it, dof, entType );
        allEqns.insert(tmp.Begin(), tmp.End());
      }
      eqns.Clear();
      eqns.Resize(allEqns.size());
      std::copy(allEqns.begin(), allEqns.end(), eqns.Begin());
    }
    
    
    void FeSpace::GetEntityListEqns(StdVector<Integer>& eqns, 
                                    shared_ptr<EntityList> ent,
                                    BaseFE::EntityType entType) {
      std::set<Integer> allEqns;
      StdVector<Integer> tmp;
      EntityIterator it = ent->GetIterator();
      for( it.Begin(); ! it.IsEnd(); it++ ) {
        tmp.Init();
        GetEqns( tmp, it, entType);
        allEqns.insert(tmp.Begin(), tmp.End());
      }
      eqns.Clear();
      eqns.Resize(allEqns.size());
      eqns.Init();
      std::copy(allEqns.begin(), allEqns.end(), eqns.Begin());
    }
      

    //! Get Equation numbers for a specific element
    void FeSpace::GetElemEqns(StdVector<Integer>& eqns,const Elem* elem, UInt dof){
      EXCEPTION("Not implemented");
    }

    void FeSpace::AddFeFunction( shared_ptr<BaseFeFunction> fct ){
      feFunction_ = fct;
      return;
    }
  
   UInt FeSpace::GetNumFunctions( const EntityIterator ent ){
     BaseFE * fe = GetFe(ent);
     return fe->GetNumFncs();
   }
   
   //=======================================================
   //Perform Nodal equation Numbering
   //=======================================================
  void FeSpace::MapNodalEqns(UInt phase) {
    LOG_DBG(feSpace) << "Mapping nodal Eqns, Phase " << phase;
    UInt dofsPerUnknown = 0;
    shared_ptr<ResultInfo> feFctResult;
    StdVector< shared_ptr<EntityList> > fctEntList;
    StdVector< shared_ptr<EntityList> >::iterator entIt;
    boost::unordered_map< Integer , StdVector<BcType> >::iterator bcIt;
    UInt actNode = 0;
    NodeTypeMap::const_iterator it;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    
    switch(phase) {
      case 1:
        // number the nodal equations if they are not contained in the BcKeys
        // due to our virtual node array we have a valid structure

        // Get result for the feFunction
        feFctResult = feFct->GetResultInfo();
        
        // Get number of dofs
        dofsPerUnknown = GetNumDofs();
        
        // now loop over all nodes and assign an equation number
        it = nodesType_.begin();
        for( ; it != nodesType_.end(); ++it ) {
          actNode = it->first;
          LOG_DBG3(feSpace) << "actNode=" << actNode;
          // create entry if it does not exist
          if(nodeMap_.eqns.find(actNode) == nodeMap_.eqns.end()){
            LOG_DBG3(feSpace) << "not found -> create";
            nodeMap_[actNode] = StdVector<Integer>(dofsPerUnknown);
            // let's hope, we never have more nodes than INT_MAX...
            nodeMap_[actNode].Init(INT_MAX);
          }

          for(UInt iDof = 0; iDof < dofsPerUnknown; iDof ++){
            if(nodeMap_.BcKeys.find(actNode) != nodeMap_.BcKeys.end()){
              // node has no bc
              if(nodeMap_.BcKeys[actNode][iDof] == NOBC){
                // dof has no bc
                nodeMap_[actNode][iDof] = ++numEqns_;
                numFreeEquations_++;
                LOG_DBG3(feSpace) << "no bc -> eqn num=" << numEqns_;
              }
            } else if(nodeMap_[actNode][iDof] == INT_MAX) {
              // node is initialized but not set
              nodeMap_[actNode][iDof] = ++numEqns_;
              numFreeEquations_++;
              LOG_DBG3(feSpace) << "not set -> eqn num=" << numEqns_;
            }
          }
        }
        
        break;
        
      case 2:
        bcIt = nodeMap_.BcKeys.begin();
        while(bcIt != nodeMap_.BcKeys.end()){
          //TODO: take the kind of BC into account
          //OPTIONAL: Delete the BcKeys Map afterwards?!
          UInt vNode = bcIt->first;
          StdVector<BcType> & dofs = bcIt->second;
          
          LOG_DBG3(feSpace) << "vNode: " << vNode
                            << ", eqns: " << nodeMap_[vNode].ToString();
          // set dirichlet bc's
          for(UInt iDof =0; iDof < dofs.GetSize(); iDof++){
            if(dofs[iDof]!=NOBC){
               //nodeMap_[bcIt->first] = ++numEqns_;
              if(dofs[iDof] == IDBC){
                nodeMap_[vNode][iDof] = ++numEqns_;
                LOG_DBG3(feSpace) << "idbc -> eqn num=" << numEqns_;
              } else if(dofs[iDof] == HDBC){
                nodeMap_[vNode][iDof] = 0 ;
                LOG_DBG3(feSpace) << "hdbc -> eqn num=0";
              } 
              //else if(bcIt->second[iDof] == CONSTRAINT){
              //  Integer masterNode = nodeMap_.slaveMasterNodes[bcIt->first];
              //  nodeMap_[bcIt->first] = -1*masterNode;
              //}
            }
          }
          bcIt++;
        }
        {
          // set periodic bc's
          // the map is ordered by C++, so we do not know, in which order bc's are applied!
          std::queue<std::pair<Integer,Integer>> unprocessed;
          std::map<std::pair<Integer,Integer>, std::pair<Integer,Integer>  >::iterator  conNodes = nodeMap_.constraintNodes.begin();
          for(; conNodes !=nodeMap_.constraintNodes.end() ; conNodes++){

            std::pair<Integer,Integer> slaveP = conNodes->first;
            std::pair<Integer,Integer> masterP = conNodes->second;
            
            // attention: we can not rely on the fact, that the nodes mentioned in the
            // constraint nodes list are really contained in our region list
            if( nodeMap_.eqns.find(masterP.first) == nodeMap_.eqns.end() )  {
              EXCEPTION("Master node of constraints could not be found");
            }
            
            if( nodeMap_.eqns.find(slaveP.first) == nodeMap_.eqns.end() )  {
              EXCEPTION("Slave node of constraints could not be found");
            }

            Integer eqn = nodeMap_[masterP.first][masterP.second];

            if(eqn == INT_MAX) {
              // master is not set -> cannot set slave, put it into queue
              LOG_DBG3(feSpace) << "slave=" << slaveP.first
                                 << " | master=" << masterP.first
                                 << " -> not set";
              unprocessed.push(slaveP);
              continue;
            }

            if(eqn >= 0) {
              // master is regular node
              nodeMap_[slaveP.first][slaveP.second] = -1*eqn;
            } else {
              // master is itself a slave of another node
              nodeMap_[slaveP.first][slaveP.second] = eqn;
            }

            LOG_DBG3(feSpace) << "slave=" << slaveP.first << ", dof=" << slaveP.second
                               << " | master=" << masterP.first << ", dof=" << masterP.second
                               << " -> eqn num=" << nodeMap_[slaveP.first][slaveP.second];
          }

          while (!unprocessed.empty()) {
            // get first element from queue and remove it
            std::pair<Integer,Integer> slaveP = unprocessed.front();
            std::pair<Integer,Integer> masterP = nodeMap_.constraintNodes[slaveP];
            unprocessed.pop();

            Integer eqn = nodeMap_[masterP.first][masterP.second];

            if(eqn == INT_MAX) {
              // master is still not set -> cannot set slave,
              // put it at the end of queue
              unprocessed.push(slaveP);
            } else {
              // master is set and itself a slave, so eqn has to be negative
              assert(eqn < 0);
              nodeMap_[slaveP.first][slaveP.second] = eqn;

              LOG_DBG3(feSpace) << "slave=" << slaveP.first << ", dof=" << slaveP.second
                                << " | master=" << masterP.first << ", dof=" << masterP.second
                                << " -> eqn num=" << nodeMap_[slaveP.first][slaveP.second];
            }
          }
        }
        break;
      default:
        EXCEPTION("In FeSpace::MapNodalEqns(): Supplied wrong argument for the numbering phase");
        break;
    }
  }

  
  void FeSpace::GetOlasMappings( StdVector<AlgebraicSys::SBMBlockDef>& sbmBlocks,
                                   std::map<UInt,StdVector<std::set<Integer> > >&
                                   minorBlocks ) {
    // currently we only support "standard" solution strategy
    if( solStrat_->GetType() != SolStrategy::STD_STRATEGY ) {
      EXCEPTION( "Currently we just support the standard solution strategy for H1.");
    }
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    FeFctIdType fctId = feFct->GetFctId();
    
    // Check, if static condensation is to be performed
    bool statCond = solStrat_->UseStaticCondensation();
    if( !statCond ) {
      // -------------------------
      //  No Static Condensation
      // -------------------------
      // we have just one big matrix block, where all equations are put in
      sbmBlocks.Resize(1);
      std::set<Integer> & allEqns = sbmBlocks[0][fctId];
      for( UInt i = 0; i < numEqns_; ++i ) {
        allEqns.insert(i+1);
      }
    } else {

      // ---------------------
      //  Static Condensation
      // ---------------------
      
      minorBlocks[0]=(StdVector<std::set<Integer> >());
      // Create 2 SBM block sets:
      // 1st block: Contains all nodal, edge and face contributions
      // 2nd block: contains all element interior equations (in 2D: faces)
      
      // Preliminary first set: all equations
      // If we have later on no interior equations, we can directly return
      // this set
      std::set<Integer>  nonIntEqns;
      for( UInt i = 0; i < numEqns_; ++i ) {
        nonIntEqns.insert(i+1);
      }
      
      // Determine type of entity type to condensate
      // 2D: faces, 3D: interior
      BaseFE::EntityType intType;
      if( ptGrid_->GetDim() == 3) {
        intType = BaseFE::INTERIOR;
      } else {
        intType = BaseFE::FACE;
      }

      // list with all equation numbers for interior block
      std::set<Integer> intBlock;
      EntityNodesType intNodes = vNodesCont_[intType];
      
      // Loop over all elements / faces
      EntityNodesType::iterator intIt = intNodes.begin();
      for( ; intIt != intNodes.end(); ++intIt ) {
        StdVector<UInt> & innerNodes = intIt->second;
        LOG_DBG(feSpace) << "innerNode has size " << innerNodes.GetSize() 
                                                             << std::endl;

        // set for all inner eqns per element
        std::set<Integer> innerEqns;
        for(UInt i = 0; i < innerNodes.GetSize(); ++i ) {
          // collect entries for SBM block
          intBlock.insert( nodeMap_[innerNodes[i]].Begin(),
                           nodeMap_[innerNodes[i]].End() );
          // collect entries for minor block
          innerEqns.insert(nodeMap_[innerNodes[i]].Begin(),
                           nodeMap_[innerNodes[i]].End());
        } // loop over interior nodes
        
        if( innerEqns.size() ) {
          minorBlocks[1].Push_back(innerEqns);
        }
      } // loop: elements / faces (= "inner" entity type)

      // take set-difference:
      // so far we have collected all inner equations, now subtract them
      // from the set of "all" equations, to get non-inner equations
      std::set<Integer> diff;
      std::insert_iterator<std::set<Integer> > insert_it (diff, diff.begin());
      std::set_difference( nonIntEqns.begin(), nonIntEqns.end(),
                           intBlock.begin(), intBlock.end(),
                           insert_it );
      sbmBlocks.Resize(2);
      sbmBlocks[0][fctId] = diff;
      sbmBlocks[1][fctId] = intBlock;

    } // if static condensation

  }

  
  void FeSpace::GetAllElems(std::set<const Elem*>& allElems ) {
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock();
    StdVector< shared_ptr<EntityList> > entityLists = feFct->GetEntityList();
    for( UInt iList = 0; iList < entityLists.GetSize(); ++iList ) {
      EntityList&  actList = *(entityLists[iList]);
      // check, if entitylist contains elements
      if( actList.GetType() == EntityList::ELEM_LIST ||
          actList.GetType() == EntityList::SURF_ELEM_LIST || 
          actList.GetType() == EntityList::NC_ELEM_LIST ) {
        EntityIterator elemIt = actList.GetIterator();
        for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++) {
          allElems.insert(elemIt.GetElem());
        }
      }
    }
  }
    
  void FeSpace::PrintEqnMap(std::ostream* file)
  {
    std::ostream& out = file != NULL ? *file : std::cout;

    // obtain (fctId,eqnNr) -> (sbm,index) mapping from OLAS. Note that indices are 1-based!
    StdVector<UInt> blockNums, indices;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    FeFctIdType fctId = feFct->GetFctId();
    AlgebraicSys * algSys = feFct->GetSystem();
    if( algSys )
      algSys->MapCompleteFctIdToIndex(fctId, blockNums, indices);
    
    std::string resultName = SolutionTypeEnum.ToString(feFct->GetResultInfo()->resultType);
    
    out << "*****************************************\n";
    out << " F E - S P A C E - I N F O R M A T I O N \n";
    out << "*****************************************\n";
    out << "Physical Quantity: " << resultName << std::endl;
    out << "FeFunction Id: " << fctId << std::endl << std::endl;
    
    out << "size of indices: " << indices.GetSize() << std::endl;
    out << "indices=" << indices.ToString() << std::endl << std::endl;
    
    // =================================
    // 1) ELEMENT INFORMATION
    // =================================
    
    // iterate over all elements
    std::set<const Elem*> allElems;
    GetAllElems(allElems);
  
    // loop over all elements and print geometric information
    std::set<const Elem*>::iterator elemIt;
    for( elemIt = allElems.begin(); elemIt != allElems.end(); ++elemIt ) {
      // print element information (type, region, connect, edges, faces)
      const Elem* ptElem = *elemIt;
      
      out << "Elem #" << ptElem->elemNum << std::endl;
      out << "Type: " << Elem::feType.ToString( ptElem->type ) << std::endl;
      out << "Connect: " << ptElem->connect.ToString( 0 ) << std::endl;
      
      // Print edge  information
      out << "Edges: ";
      for( UInt i=0, numEdges = ptElem->extended->edges.GetSize(); i < numEdges; ++i ) {
        StdVector<UInt> edgeNodes;
        Integer edgeNum = ptElem->extended->edges[i];
        ptElem->GetEdgeNodes( std::abs(edgeNum) , edgeNodes );
        out << "E #" << edgeNum << " (" << edgeNodes[0] << "-> " << edgeNodes[1] << "), ";
      }
      
      // Print face  information
      out << std::endl << "Faces: ";
      for( UInt i=0, numFaces = ptElem->extended->faces.GetSize(); i < numFaces; ++i ) {
        StdVector<UInt> faceNodes;
        UInt faceNum = ptElem->extended->faces[i];
        ptElem->GetFaceNodes( faceNum, faceNodes );
        out << "F #" << faceNum << " (" << faceNodes.ToString( 0 ) << "), ";
      }
      out << "\n\n";

      // print header
      out << "\t\t#num\tvNodes\tEqnNrs\tSBM\tindex\n";
      out << "\t\t=================================================================\n";

      std::string prefix = "\t\t\t";
      
      // print vertex / edge / face / inner information
      StdVector<BaseFE::EntityType> entTypes;
      ElemShape& shape = Elem::shapes[ptElem->type];
      entTypes = BaseFE::VERTEX, BaseFE::EDGE, BaseFE::FACE, BaseFE::INTERIOR;
      for( UInt iType = 0; iType < entTypes.GetSize(); ++iType ) {
        
        BaseFE::EntityType type = entTypes[iType];
        
        // vector containing entity numbers
        StdVector<UInt> entNumbers;
        EntityNodesType &vNodes = vNodesCont_[type]; 
        
        // fill vector containing entity types
        if( type == BaseFE::VERTEX ) {
          entNumbers = ptElem->connect;
          // note: we may only take the number of nodes, which are 
          // also used for the calculation element
          if( mapType_ != GRID ) {
            entNumbers.Resize( shape.numVertices );
          }
        } else if( type == BaseFE::EDGE ) {
          entNumbers.Resize(ptElem->extended->edges.GetSize());
          for( UInt i=0; i < ptElem->extended->edges.GetSize(); ++i )
            entNumbers[i] = std::abs(ptElem->extended->edges[i]);
        } else if( type == BaseFE::FACE ) {
          entNumbers.Resize(ptElem->extended->faces.GetSize());
          for( UInt i=0; i < ptElem->extended->faces.GetSize(); ++i )
            entNumbers[i] = ptElem->extended->faces[i];
        } else if( type == BaseFE::INTERIOR ) {
          // only treat "interior", if we have any unknowns at all assigned
          // (= size of vNodes != 0)
          if( vNodes.size() ) {
            entNumbers.Resize(1);
            entNumbers.Init(ptElem->elemNum);
          }
        }

        // if any nodes of the given are available
        if( vNodes.size() ) {
          out << iType+1 << ") " << BaseFE::entityType.ToString(type) << std::endl;
          out << "========\n";
        } else {
          continue;
        }
        // loop over all entities
        for( UInt i = 0; i < entNumbers.GetSize(); ++i ) {
          UInt entNumber = entNumbers[i];
          out << "\t\t#" << entNumber << "\t";

          // Get hold of virtual nodes for given entity type
          StdVector<UInt>& vNodesEnt = vNodes[entNumber];
          UInt entNumNodes = vNodesEnt.GetSize();
          
          // check, if any virtual nodes are assigned at all
          if( vNodesEnt.GetSize() == 0 ) {
            out << "-\t-\t-\t-\n";
          }
          // leave, virtual node numbers are assigned
          for( UInt j = 0; j < entNumNodes; ++j ) {

            // print virtual node only for first entry            
            if( j > 0 ) {
              out << prefix;
            }
            // print virtual node
            out << vNodesEnt[j] << "\t";
            

            // equation numbers (loop)
            StdVector<Integer> & eqns = nodeMap_[vNodesEnt[j]];
            // immediately begin new line, if entity has no equations
            if( eqns.GetSize() == 0 ){
              out << "\n";
            }
            
            for( UInt iEqn = 0; iEqn < eqns.GetSize(); ++iEqn ) {
              Integer & eqn = eqns[iEqn];

              // indent succeeding equations correctly
              if( iEqn > 0 ) {
                out << prefix << "\t";
              }

              //equation number
              out << eqn << "\t";

              if( eqn > 0 && algSys) {
                // SBM-Block
                out << blockNums[eqn-1] << "\t";

                // index
                out << indices[eqn-1] <<  " (" << indices[eqn-1]-1 << ")\n";
              } else {
                out << "-\t-\n";
              }
            }
          } // loop over virtual nodes
        } // loop over entity numbers
        out << "\n";
      } // loop over entity types
      out << "\n\n";
    } // loop over elements


    // =================================
    // 2) Global Equation numbering
    // =================================
    shared_ptr<ResultInfo> feFctResult = feFct->GetResultInfo();
    // ---------------
    //  NODES
    // ---------------
    boost::unordered_map< Integer , StdVector<Integer> >::const_iterator nodeIt = nodeMap_.eqns.begin();
    boost::unordered_map< Integer , StdVector<BcType> >::iterator nodeBcIt;

    out << "EQUATION MAPPING" << std::endl << std::endl;
    out << "nodeNr | type | comp | eqnNr | SBM |  index 1 / 0    |  BC " << std::endl;
    out << "------------------------------------------------------------" << std::endl;
    while(nodeIt != nodeMap_.eqns.end())
    {
      nodeBcIt = nodeMap_.BcKeys.find(nodeIt->first);
      for(UInt iDof =0; iDof < nodeIt->second.GetSize(); iDof++){

        const int& eqn = nodeIt->second[iDof];

        // virtual node number (only once for all dofs) and type
        if(iDof == 0)
          out << std::setw(6) << nodeIt->first << " |   " << BaseFE::entityType.ToString(nodesType_[nodeIt->first])[0];
        else
          out << "       |    ";

        // component and equation
        out << "  |" << std::setw(5) << feFctResult->dofNames[iDof];
        out << " |" << std::setw(6) << eqn;

        // sbm and index in original 1-based plus 0-based
        if(eqn == 0 || !algSys)
          out << " |   - |        -       ";
        else
          out << " | " << std::setw(3) << blockNums[std::abs(eqn)-1] << " | "
               << std::setw(6) << indices[std::abs(eqn)-1] << " / "
               << std::setw(6) << indices[std::abs(eqn)-1]-1;

        // bc type
        out << " | ";
        if(nodeBcIt != nodeMap_.BcKeys.end() && nodeBcIt->second[iDof] != NOBC)
            out << BcTypeEnum.ToString(nodeBcIt->second[iDof]);
        out << std::endl;
      }
      nodeIt++;
    }
  }
  


  void FeSpace::CreateEquIndGeomMap(boost::unordered_map< Integer, BaseFeFunction::EqNodeGeom>& eqIndGeomMap,
                                    UInt& maxEqn, UInt& dim){
    // since it's only used for lowest order FE-functions, the node-entity is
    // whether vertex or edge (for edge-elements)

    // obtain (fctId,eqnNr) -> (sbm,index) mapping from OLAS
    StdVector<UInt> blockNums, indices;
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    FeFctIdType fctId = feFct->GetFctId();
    AlgebraicSys * algSys = feFct->GetSystem();
    if( algSys )
      algSys->MapCompleteFctIdToIndex(fctId, blockNums, indices);
    if( this->GetSpaceType() != this->HCURL ){
      /************** H1 Lagrange - Version *******************/
      // 1) loop over every node and store it's equations,
      //    index and grid-geometry
      boost::unordered_map< Integer , StdVector<Integer> >::const_iterator nodeIt = nodeMap_.eqns.begin();
      boost::unordered_map< Integer , StdVector<BcType> >::iterator nodeBcIt;
      maxEqn = 0;
      while(nodeIt != nodeMap_.eqns.end()){
        BaseFeFunction::EqNodeGeom nMap;
        dim = nodeIt->second.GetSize();
        Vector<Double> c;
        feFct->GetGrid()->GetNodeCoordinate(c, nodeIt->first);
        // only if equation number != 0 (they start with 1)
        // this would indicate that we have a Dirichlet (inh. or hom.)
        // but they are not handled inside AMG!!!

        switch (dim){
        case 1:
          if( nodeIt->second[0]){
            for(UInt iDof = 0; iDof < dim; iDof++){
              // index
              nMap.indexNum = indices[nodeIt->second[iDof] - 1 ]; // index[i] = eqn[i] - 1
              //nodeNum
              nMap.nodeNum = nodeIt->first;
              //coord
              nMap.coord = c;
              if(nodeIt->second[iDof] > (Integer)maxEqn) maxEqn = nodeIt->second[iDof];
              // insert it into map
              eqIndGeomMap[nodeIt->second[iDof]] = nMap; //nodeIt->second[iDof] is the eqn
            }
          }
          break;

        case 2:
          if( nodeIt->second[0] || nodeIt->second[1]){
            for(UInt iDof = 0; iDof < dim; iDof++){
              // index
              nMap.indexNum = indices[nodeIt->second[iDof] - 1 ]; // index[i] = eqn[i] - 1
              //nodeNum
              nMap.nodeNum = nodeIt->first;
              //coord
              nMap.coord = c;
              if(nodeIt->second[iDof] > (Integer)maxEqn) maxEqn = nodeIt->second[iDof];
              // insert it into map
              eqIndGeomMap[nodeIt->second[iDof]] = nMap; //nodeIt->second[iDof] is the eqn
            }
          }
          break;

        case 3:
          if( nodeIt->second[0] || nodeIt->second[1] || nodeIt->second[2] != 0
          ){
            for(UInt iDof = 0; iDof < dim; iDof++){
              // index
              nMap.indexNum = indices[nodeIt->second[iDof] - 1 ]; // index[i] = eqn[i] - 1
              //nodeNum
              nMap.nodeNum = nodeIt->first;
              //coord
              nMap.coord = c;
              if(nodeIt->second[iDof] > (Integer)maxEqn) maxEqn = nodeIt->second[iDof];
              // insert it into map
              eqIndGeomMap[nodeIt->second[iDof]] = nMap; //nodeIt->second[iDof] is the eqn
            }
          }
          break;
        }//switch dim
        nodeIt++;
      }//while nodeIt != nodeMap_.eqns.end()
    }else{
      if( this->GetSpaceType() == this->HCURL ){
        /************** HCurl Nedelec - Version *******************/
        dim = 1;
        std::map<UInt,StdVector<UInt> > indexGeomMap;
        // iterate over all elements
        std::set<const Elem*> allElems;
        GetAllElems(allElems);

        BaseFeFunction::EqNodeGeom nMap;
        // loop over all elements and store geometric information
        std::set<const Elem*>::iterator elemIt;
        for( elemIt = allElems.begin(); elemIt != allElems.end(); ++elemIt ) {

          const Elem * ptElem = *elemIt;

          BaseFE::EntityType type;
          type = BaseFE::EDGE;

          // vector containing entity numbers
          StdVector<UInt> entNumbers;
          EntityNodesType &vNodes = vNodesCont_[type];


          // loop over all edges
          for( UInt i = 0; i < ptElem->extended->edges.GetSize(); ++i ) {
            UInt edgeNum = std::abs(ptElem->extended->edges[i]);
            // direction of edge
            UInt dir = ptElem->extended->edges[i] < 0 ? 0 : 1;
            // Get hold of virtual nodes for given edge
            StdVector<UInt>& vNodesEdge = vNodes[edgeNum];
            UInt edgeNumNodes = vNodesEdge.GetSize();
            // check, if any virtual nodes are assigned at all
            if( vNodesEdge.GetSize() != 0 ) {
              for( UInt j = 0; j < edgeNumNodes; ++j ) {
                // equation numbers (loop)
                StdVector<Integer> & eqns = nodeMap_[vNodesEdge[j]];

                //geometry information
                StdVector<UInt> edgeNodes;
                StdVector<UInt> edgeNodesC; //for getting the coordinates
                ptElem->GetEdgeNodes( edgeNum , edgeNodes );

                //for 0th order edge-elements (=Nédélec elements) there
                //should only be one equation per edge
                if(eqns.GetSize() > 1){
                  EXCEPTION("FeSpace.cc: Index<->Geometry mapping for edges only works for 0th order EdgeElements!");
                }else{
                  //for( UInt iEqn = 0; iEqn < eqns.GetSize(); ++iEqn ) {
                  Integer & eqn = eqns[0];
                  if( eqn > 0 && algSys) {
                    // index
                    UInt index = indices[eqn-1];
                    // Assign index and edge-number to the map only if it is not already mapped
                    if(indexGeomMap.find(index) == indexGeomMap.end()){
                      nMap.indexNum = index;
                      StdVector< Vector<Double> > nCoords;
                      feFct->GetGrid()->GetNodeCoordinates(nCoords, edgeNodes, false);

                      nMap.eCoords.Resize(2);
                      nMap.eCoords[0] = nCoords[0];
                      nMap.eCoords[1] = nCoords[1];
                      nMap.eNodes.Resize(2); // 2 nodes per edge
                      switch(dir)
                      {
                      case 1:
                        indexGeomMap[index].Push_back(edgeNodes[0]);
                        indexGeomMap[index].Push_back(edgeNodes[1]);
                        nMap.eNodes[0] = edgeNodes[0];
                        nMap.eNodes[1] = edgeNodes[1];
                        break;
                      case 0:
                        indexGeomMap[index].Push_back(edgeNodes[1]);
                        indexGeomMap[index].Push_back(edgeNodes[0]);
                        nMap.eNodes[0] = edgeNodes[1];
                        nMap.eNodes[1] = edgeNodes[0];
                        break;
                      }
                      // insert it into final map
                      eqIndGeomMap[eqns[0]] = nMap; //nodeIt->second[iDof] is the eqn
                    }
                  }//if( eqn > 0 && algSys)
                } // loop over virtual nodes
              }// loop over nodes of edge i
            } // if vNodes.GetSize != 0
          } // loop over edges
        }// loop over elements
      }else EXCEPTION("AMG: FeSpace is not H1 nor HCurl !!!");
    }
  }


  UInt FeSpace::GetNumDofs() const {
    // Just return number of components
    shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
    assert(feFct);
    return feFct->GetResultInfo()->dofNames.GetSize();
  }
  

  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************

  // Definition of finite element space types
   static EnumTuple spaceTypeTuples[] = {
     EnumTuple(FeSpace::UNDEF_SPACE, "Undef"), 
     EnumTuple(FeSpace::CONSTANT,    "Constant"), 
     EnumTuple(FeSpace::H1,          "H1"),
     EnumTuple(FeSpace::HCURL,       "HCurl"),
     EnumTuple(FeSpace::HDIV,        "HDiv"),
     EnumTuple(FeSpace::L2,          "L2"),
   };
   Enum<FeSpace::SpaceType> FeSpace::SpaceTypeEnum = \
      Enum<FeSpace::SpaceType>("Types of FE Spaces",
          sizeof(spaceTypeTuples) / sizeof(EnumTuple),
          spaceTypeTuples);
   
   // Definition of polynomial types
    static EnumTuple polyTypeTuples[] = {
      EnumTuple(FeSpace::UNDEF_POLY, "Undef"), 
      EnumTuple(FeSpace::LAGRANGE,   "Lagrange"), 
      EnumTuple(FeSpace::LEGENDRE,   "Legendre"),
      EnumTuple(FeSpace::JACOBI,     "Jacobi"),
    };
    Enum<FeSpace::PolyType> FeSpace::PolyTypeEnum = \
       Enum<FeSpace::PolyType>("Types of FE Spaces",
           sizeof(polyTypeTuples) / sizeof(EnumTuple),
           polyTypeTuples);
   
    // Definition of integration oder mode types
    static EnumTuple integModeTuples[] = {
      EnumTuple(FeSpace::INTEG_MODE_ABSOLUTE, "absolute"),
      EnumTuple(FeSpace::INTEG_MODE_RELATIVE, "relative")
    };
    Enum<FeSpace::IntegOrderMode> FeSpace::IntegOrderModeEnum = \
       Enum<FeSpace::IntegOrderMode>("Types of FE Spaces",
           sizeof(integModeTuples) / sizeof(EnumTuple),
           integModeTuples);
  
  // Definition of types of boundary conditions
  static EnumTuple bcTypeTuples[] = {
    EnumTuple(FeSpace::NOBC,    "NOBC"), 
    EnumTuple(FeSpace::HDBC,    "HDBC"),
    EnumTuple(FeSpace::IDBC,    "IDBC"),
    EnumTuple(FeSpace::CS,    "CS")
  };
  
  Enum<FeSpace::BcType> FeSpace::BcTypeEnum = \
     Enum<FeSpace::BcType>("Types of boundary conditions",
         sizeof(bcTypeTuples) / sizeof(EnumTuple),
         bcTypeTuples);
} // end of namespace 
