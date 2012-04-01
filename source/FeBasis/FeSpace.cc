#include "FeSpace.hh"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "DataInOut/Logging/LogConfigurator.hh"

// include headers of sub-classes for factory method
#include "H1/FeSpaceH1Nodal.hh"
#include "L2/FeSpaceL2Nodal.hh"
#include "H1/FeSpaceH1Hi.hh"
#include "HCurl/FeSpaceHCurlHi.hh"


namespace CoupledField {


  // declare class specific logging stream
  DECLARE_LOG(feSpace)
  DEFINE_LOG(feSpace, "feSpace")
  
  
  //! Constructor
  FeSpace::FeSpace( PtrParamNode paramNode, PtrParamNode infoNode, 
                    Grid* ptGrid ){
    
    myParam_ = paramNode;
    infoNode_ = infoNode;
    ptGrid_ = ptGrid;
    isFinalized_ = false;
    isContinuous_ = true;
    numEqns_ = 0;
    numFreeEquations_ = 0;
    solStep_ = 1;
    orderOffset_ = 0;

    bcCounter_[NOBC] = 0;
    bcCounter_[HDBC] = 0;
    bcCounter_[IDBC] = 0;
    bcCounter_[CS] = 0;
    
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
                                               Grid* ptGrid ) {
     

    /* One big Problem> Due to the splitting of spaces in Hi and Lo/Lagrange, it is not
     * possible to deal with Diffent Polynomial types in different regions
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
    PtrParamNode polyNode = param->Get("fePolynomialList", ParamNode::PASS );
    if(!polyNode){
      WARN("No Polynomial specified falling back to Defaults");
      polyType = LAGRANGE;
      LOG_DBG(feSpace) << "No explicit definition available, using default";
    }else{
      
      
      // loop over all polynomial definitions, get the id and query, if the
      // polynomial gets also referenced in the <reegionList> of this PDE.
      PtrParamNode regList = aNode->Get("regionList");
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
        if( polyType == LAGRANGE )  {
          // create explicit lower order HCurl space
          EXCEPTION("Explicit lower order space H-Curl not defined yet");
        } else {
          ret.reset(new FeSpaceHCurlHi(aNode, infoNode, ptGrid));
        }
        break;
      case CONST:
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
      feFunction_->GetGrid()->GetNodesByName( tempNodes, name );
      
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
      //          initially.
    } else {
      EntityList::ListType defType = ent->GetType();
      EntityIterator entIt = ent->GetIterator();
      nodes.Resize(0);
      // UInt nCount = 0;
      switch (defType){
        case EntityList::REGION_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<Elem*> eList;
            feFunction_->GetGrid()->GetElems(eList,entIt.GetRegion());
            StdVector<UInt> eNodes;
            for (UInt i = 0; i < eList.GetSize(); i ++ ){

              GetNodesOfElement(eNodes,eList[i],entType);
              for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
                nodes.Push_back(eNodes[aNode]);
              }
            }
          }
        break;
        case EntityList::ELEM_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            GetNodesOfElement(eNodes,entIt.GetElem(),entType);
            for (UInt aNode = 0; aNode < eNodes.GetSize(); aNode += 1 ){
              nodes.Push_back(eNodes[aNode]);
            }
          }
        break;
        case EntityList::NC_ELEM_LIST:
        case EntityList::SURF_ELEM_LIST:
          for(entIt.Begin(); !entIt.IsEnd(); entIt++){
            StdVector<UInt> eNodes;
            GetNodesOfElement(eNodes,  entIt.GetSurfElem() ,entType);
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
    UInt elemNum = ptElem->elemNum;
    if(virtualNodes_.find(elemNum) ==virtualNodes_.end()){

      EXCEPTION("FeSpace::GetNodesOfElement: Could not find requested element #"
          << ptElem->elemNum << " of region " 
          <<  ptGrid_->GetRegion().ToString(ptElem->regionId));
    }
    if(entType == BaseFE::ALL){
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
    }else{
      nodes.Resize(virtualNodes_[elemNum][entType].vNodes.GetSize());
      const StdVector<UInt>& entNodes =  virtualNodes_[elemNum][entType].vNodes;
      for (UInt i = 0; i < entNodes.GetSize(); i++ ){
        nodes[i] =  entNodes[i];
      }
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
    //                  order to 2*(p+1) + relativeOrder
    //    b) ANISOTROPIC: get maximum order in each direction
    //                    and set integration oder to
    //                    2*(p_maxLicDir) + relativeOrder_locDir
    method = id.method;
    LOG_DBG2(feSpace) << "\tmethod: " 
                      << IntScheme::IntegMethodEnum.ToString(method);
    if( id.mode == ABSOLUTE ) {
      // ABSOLUTE order given
      order = id.order;
      LOG_DBG2(feSpace) << "\torder (absolute):" << order.ToString();
    } else {
      // RELATIVE order given
      LOG_DBG2(feSpace) << "\torder (relative): " << id.order.ToString();
      UInt p = fe->GetMaxOrder();
      order = id.order;
      order += 2*(p+1);
      LOG_DBG2(feSpace) << "\torder (final): " << order.ToString();
    }
  }

  void FeSpace::CreateVirtualNodes(){
    
    
    // store all regions the space is defined on
    regions_ = feFunction_->GetRegions();
    
    //the function checks for the special case if we got only one element
    //in the region map and if this mapping is of type GRID
    //if so, we call a specialized but efficient function
    // if not, we get to the general case
    if(mapType_ == GRID && isContinuous_){
      CreateGridNodes();
    }else{
      CreatePolynomialNodes();
      //we have the general case so we iterate over all elements
    }
  }

  void FeSpace::CreateGridNodes(){

    LOG_TRACE(feSpace) << "starting to create virtual nodes based on GRID";

    StdVector< shared_ptr<EntityList> > fctEntList;
    StdVector<UInt> curNodes;

    fctEntList = feFunction_->GetEntityList();

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
      feFunction_->GetGrid()->GetNodesByName( curNodes, name );
      //GetNodesOfEntities( curNodes, actElemList );
      for ( UInt aNode= 0; aNode < curNodes.GetSize(); aNode++ ) {
        if(gridToVirtualNodes_.find(curNodes[aNode]) == gridToVirtualNodes_.end()){
          gridToVirtualNodes_[curNodes[aNode]].Push_back(curNodes[aNode]);
        }
      }

      // loop over all elements a little overhead but makes our life easier later
      // otherwise one would have to distinguish everytime for this special case
      // for the future we should consider also to split the in vertex, edge, etc nodes
      // perhaps the meshing tools are getting better for higher order elements...
      EntityIterator entIt = actElemList->GetIterator();
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        const Elem* actEl = entIt.GetElem();
        const StdVector<UInt> & elemNodes = actEl->connect;
        virtualNodes_[actEl->elemNum][BaseFE::VERTEX].vNodes = elemNodes;
        virtualNodes_[actEl->elemNum][BaseFE::VERTEX].offset = StdVector<UInt>(elemNodes.GetSize());
        virtualNodes_[actEl->elemNum][BaseFE::VERTEX].offset.Init(1);
      }
    }

    for(std::map<UInt,StdVector<UInt> >::const_iterator it = gridToVirtualNodes_.begin(); it != gridToVirtualNodes_.end(); ++it) {
      //here we can hardcode to zero as we assume continuous approximation
      nodesType_[it->second[0]] = BaseFE::VERTEX;
    }
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

    //ok these data structures are a bit messy but this is the most obvious typ
    // and they are temporary anyway. furthermore we run here once in a lifetime
    std::map< UInt, std::map<UInt, StdVector<Integer> > > vertexNodes;
    std::map< UInt, std::map<UInt, StdVector<Integer> > > edgenodes;
    std::map< UInt, std::map<UInt, StdVector<Integer> > > facenodes;
    std::map<UInt, StdVector<Integer> > interiornodes;
    //Different lists have to be treated differently
    EntityList::ListType actListType = EntityList::NO_LIST;

    LOG_TRACE(feSpace) << "starting to create virtual nodes based on POLYNOMIAL";
    
    fctEntList = feFunction_->GetEntityList();

    //get the highest possible node number
    //UInt offset = feFunction_->GetGrid()->GetNumNodes();
    UInt offset =0;

    //Stores the current order
    // MappingType curMap = GRID;
    
    //changed algorithm first we add the volume elements and later on the surfaces
    // so we loop twice over the entities

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
      LOG_DBG(feSpace) << "treating entityList '" << fctEntList[actList]->GetName() << "'";
      //cast down to element list
      EntityList* actElemList = fctEntList[actList].get();
//      RegionIdType curReg = actElemList->GetRegion();

      // curMap = POLYNOMIAL;

      // Get iterator of current element list
      EntityIterator entIt = actElemList->GetIterator();
      
      // loop over all elements
      for(entIt.Begin(); !entIt.IsEnd();entIt++){

        // Fetch current finite element. This is performed by the specialized 
        // version, so this element "knows" already about its order, unknowns etc.
        BaseFE* ptFe = GetFe( entIt );
        const Elem* actEl = entIt.GetElem();  




        StdVector<UInt> permutations; // initially size 0
        UInt elemNum = actEl->elemNum;
        
        //===========================================================
        //Assign the BaseFE::VERTEX node numbers
        //===========================================================
        LOG_DBG2(feSpace) << "mapping vertex nodes";
        UInt numVertexNodes = 0;
        UInt numVert = Elem::shapes[actEl->type].numVertices;
        StdVector<UInt> elemNodes = actEl->connect;

        EntityTypeNodes & vtn =  virtualNodes_[actEl->elemNum][BaseFE::VERTEX];

        // check, if the vertices of this element were already numbered
        if( vtn.vNodes.GetSize() == 0 ) {
          for ( UInt iVert= 0; iVert< numVert; iVert++ ) {
            UInt vertexNum = elemNodes[iVert];
            ptFe->GetNodalPermutation(permutations,actEl,BaseFE::VERTEX,iVert);
            numVertexNodes = permutations.GetSize();

            if(isContinuous_){
              //in the continuous case we need to check if we already have an entry for
              //this vertex
              if(vertexNodes[vertexNum].size()>0){
                //so we need to redefine the volElemNumber
                elemNum = vertexNodes[vertexNum].begin()->first;
              }
            }

            // Check if the vertex is already numbered.
            if( vertexNodes[vertexNum][elemNum].GetSize() == 0 ) {

              vertexNodes[vertexNum][elemNum].Resize(numVertexNodes);
              vertexNodes[vertexNum][elemNum].Init();
              for( UInt vertNode = 0; vertNode < numVertexNodes; ++vertNode ) {
                vertexNodes[vertexNum][elemNum][vertNode] = ++offset;
                LOG_DBG3(feSpace) << "adding " << offset << " to node_";
                nodesType_[offset] = BaseFE::VERTEX;
              }
            }


            for( UInt i = 0; i < numVertexNodes; ++i ) {
              vtn.vNodes.Push_back(vertexNodes[vertexNum][elemNum][permutations[i] ]);
              LOG_DBG3(feSpace) << "adding " << vertexNodes[vertexNum][elemNum][permutations[i] ]
                                                                                   << " as virtual vertex node to element " << actEl->elemNum;
            }
            vtn.offset.Push_back( permutations.GetSize() );

            if(isContinuous_){
              if(gridToVirtualNodes_.find(vertexNum) == gridToVirtualNodes_.end()){
                LOG_DBG3(feSpace) << "gridToVirtualNodes[" << vertexNum << "] = " << offset;
                gridToVirtualNodes_[vertexNum].Push_back(offset);
              }
            }else{
              LOG_DBG3(feSpace) << "gridToVirtualNodes[" << vertexNum << "] = " << offset;
              gridToVirtualNodes_[vertexNum].Push_back(offset);
            }
          } // loop over vertices
        }
        
        feFunction_->GetGrid()->MapEdges();
        feFunction_->GetGrid()->MapFaces();
        
        //===========================================================
        //Assign the Edge node numbers
        //===========================================================
        UInt numEdgeNodes = 0;
        ElemShape actShape = Elem::shapes[actEl->type];
        EntityTypeNodes & etn =  virtualNodes_[actEl->elemNum][BaseFE::EDGE];
        
        // check if edges of this element were already numbered
        if( etn.vNodes.GetSize() == 0 ) {
          for ( UInt iEdge=0; iEdge < actShape.numEdges; iEdge++) {
            UInt edgeNum = std::abs(actEl->edges[iEdge]);
            //get the permutation Vector
            ptFe->GetNodalPermutation(permutations,actEl,BaseFE::EDGE,iEdge);
            numEdgeNodes = permutations.GetSize();
            if(isContinuous_){
              //in the continuous case we need to check if we already have an entry for
              //this vertex
              if(edgenodes[edgeNum].size()>0){
                //so we need to redefine the volElemNumber
                elemNum = edgenodes[edgeNum].begin()->first;
              }
            }
            // Check if the edge is already numbered.
            // Additionally, if we have the case of discontinuous approximation,
            // we number the nodes separately for every element anyway.
            if(edgenodes[edgeNum][elemNum].GetSize() == 0 ) {
              //here we assume spectral element approximation and we have
              //order-1 nodes on the edge
              edgenodes[edgeNum][elemNum].Resize(numEdgeNodes);
              edgenodes[edgeNum][elemNum].Init();
              for ( UInt edgeNode = 0;edgeNode < numEdgeNodes ;edgeNode++ ) {
                edgenodes[edgeNum][elemNum][edgeNode] = ++offset;
                nodesType_[offset] = BaseFE::EDGE;
              }
            }

            //fill the virtual Nodes in the correct ordering

            for ( UInt i = 0; i < numEdgeNodes ; i++ ) {
              etn.vNodes.Push_back(edgenodes[edgeNum][elemNum][ permutations[i] ]);
            }
            etn.offset.Push_back( permutations.GetSize() );
          } // loop over edges
        }

        //===========================================================
        //Assign the Face node numbers
        //===========================================================
        UInt numFaceNodes = 0;
        EntityTypeNodes & ftn =  virtualNodes_[actEl->elemNum][BaseFE::FACE];
        
        // check if faces of this element ware already numbered
        if( ftn.vNodes.GetSize() == 0 ) {
          for ( UInt iFace=0; iFace < actShape.numFaces; iFace++) {
            UInt faceNum = actEl->faces[iFace];
            //get the permutation Vector
            ptFe->GetNodalPermutation(permutations,actEl,BaseFE::FACE,iFace);
            numFaceNodes = permutations.GetSize();
            if(isContinuous_){
              //in the continuous case we need to check if we already have an entry for
              //this vertex
              if(facenodes[faceNum].size()>0){
                //so we need to redefine the volElemNumber
                elemNum = facenodes[faceNum].begin()->first;
              }
            }

            // Check if the face is already numbered.
            // Additionally, if we have the case of discontinuous approximation,
            // we number the nodes separately for every element separately anyway.
            if(facenodes[faceNum][elemNum].GetSize() == 0 ){
              facenodes[faceNum][elemNum].Resize(numFaceNodes);
              for ( UInt faceNode = 0;faceNode < numFaceNodes ;faceNode++ ) {
                facenodes[faceNum][elemNum][faceNode] = ++offset;
                nodesType_[offset] = BaseFE::FACE;
              }
            }
            //fill the virtual Nodes in the correct ordering

            for ( UInt i = 0; i < numFaceNodes ; i++ ) {
              ftn.vNodes.Push_back(facenodes[faceNum][elemNum][ permutations[i] ]);
            }
            ftn.offset.Push_back( permutations.GetSize() );
          }
        }
        //===========================================================
        //Assign the Interior node numbers
        //===========================================================
        //get the permutation Vector just for the number of nodes
        ptFe->GetNodalPermutation(permutations,actEl,BaseFE::INTERIOR,0);
        UInt numIntNodes = permutations.GetSize();
        EntityTypeNodes & itn =  virtualNodes_[actEl->elemNum][BaseFE::INTERIOR];
        
        //Check if the current element got already numbered
        if(interiornodes[actEl->elemNum].GetSize() == 0){
          interiornodes[actEl->elemNum].Resize(numIntNodes);
          for ( UInt intNode = 0;intNode < numIntNodes ;intNode++ ) {
            interiornodes[actEl->elemNum][intNode] = ++offset;
            nodesType_[offset] = BaseFE::INTERIOR;
          }
        }
        //fill the virtual Nodes in the correct ordering
        
        for ( UInt i = 0; i  < numIntNodes ; i++ ) {
          itn.vNodes.Push_back(interiornodes[actEl->elemNum][ permutations[i] ]);
        }
        itn.offset.Push_back( permutations.GetSize());
      } // loop elements 
    } // loop entity lists

    LOG_TRACE(feSpace) << "finished creation of virtual nodes";
  }
  
  // ************************************************************************
  // GENERATE REGION SPECIFIC DATA AND PROCESS USER INPUT
  // ************************************************************************
  void FeSpace::SetRegionApproximation(RegionIdType region, std::string polyId, std::string integId){
    
    std::string regionName = ptGrid_->GetRegion().ToString(region);
    PtrParamNode regionNode = infoNode_->Get("regionList")->Get(regionName);
    
    RegionIdType pReg,iReg;

    //SECTION1: Polynomials
    if(polyNodes_.find(polyId)!=polyNodes_.end()){
      //the user specified a valid polyId so we read its data and call create ref elems
      //additionally we let the space read specific attributes only valid for himself
      PtrParamNode pNode = polyNodes_[polyId];
      ReadCustomAttributes(pNode,region);
      Matrix<Integer> order;
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
      EXCEPTION("The integration id does not match any in the IntegratoinSchemeList: " << integId);
    }
    polyToIntegMap[pReg].insert(iReg);
  }

  void FeSpace::SetDefaultRegionApproximation(){
    RegionIdType pReg=ALL_REGIONS,iReg=ALL_REGIONS;
    SetDefaultIntegration(infoNode_->Get("regionList")->Get("default"));
    SetDefaultElements(infoNode_->Get("regionList")->Get("default"));
    polyToIntegMap[pReg].insert(iReg);
  }
  
  const Elem* FeSpace::GetVolElem( const SurfElem* surfElem ) const {

    Elem * ret = NULL;
    boost::array<Elem*,2>::const_iterator it = surfElem->ptVolElems.begin();
    for( ; it != surfElem->ptVolElems.end(); it++ ) {
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
          << surfElem->elemNum << ". " );
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
    mode = IntegOrderModeEnum.Parse(modeStr,ABSOLUTE);

    //only isotropic order supported right now
    UInt isoOrder = node->Get("order")->As<UInt>();
    //create order Matrix, we only support isotropic interation right now
    order.SetIsoOrder( isoOrder );
  }

  void FeSpace::ReadPolyNode(PtrParamNode node, MappingType & mapType, Matrix<Integer> & order){

    bool grid = node->Get("useGridOrder",ParamNode::EX)->As<bool>();
    //determine mapping type
    mapType = (grid)? GRID : POLYNOMIAL;

    //Read isoOrder or Aniso Order
    PtrParamNode isoOrderNode = node->Get("isoOrder", ParamNode::PASS );
    if(isoOrderNode){
      Integer isoOrder = isoOrderNode->As<Integer>();
      order.Resize(1,1);
      order[0][0] = isoOrder;
    }

    
    PtrParamNode anIsoOrderNode = node->Get("anIsoOrder", ParamNode::PASS );
    if(anIsoOrderNode){
      UInt dim = ptGrid_->GetDim();
      order.Resize(3,dim);
      StdVector<std::string> dofs(3);
       anIsoOrderNode->GetValue("dof1",dofs[0],ParamNode::PASS);
       anIsoOrderNode->GetValue("dof2",dofs[1],ParamNode::PASS);
       anIsoOrderNode->GetValue("dof3",dofs[2],ParamNode::PASS);
       if( dofs[2] == "" ) {
         order.Resize(2,dim); 
       }
       if( dofs[1] == "" ) {
         order.Resize(1,dim);
       }
      char_separator<char> sep(" ");

      for(UInt i = 0;i < order.GetNumRows();i++){
        boost::tokenizer<char_separator<char> > tokens(dofs[i],sep);
        boost::tokenizer<char_separator<char> >::iterator tokIt=tokens.begin();
        for(UInt j = 0;j < order.GetNumCols();j++){
          order[i][j] = lexical_cast<UInt>(*tokIt);
          tokIt++;
        }
      }
    }
  }


  void FeSpace::ReadIntegList(){
    //obtain the fePolynomialList
    PtrParamNode integNode = param->Get("integrationSchemeList", ParamNode::PASS );
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
        PtrParamNode polyNode = param->Get("fePolynomialList", ParamNode::PASS );
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

  // ************************************************************************
  // ENUM INITIALIZATION
  // ************************************************************************

  // Definition of finite element space types
   static EnumTuple spaceTypeTuples[] = {
     EnumTuple(FeSpace::UNDEF_SPACE, "Undef"), 
     EnumTuple(FeSpace::CONST,       "Const"), 
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
      EnumTuple(FeSpace::ABSOLUTE, "absolute"),
      EnumTuple(FeSpace::RELATIVE, "relative")
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
