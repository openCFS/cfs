#include "FeSpaceHi.hh"

#include "FeHi.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
namespace CoupledField {

// declare class specific logging stream
DECLARE_LOG(feSpaceHi)
DEFINE_LOG(feSpaceHi, "feSpaceHi")

FeSpaceHi ::FeSpaceHi (PtrParamNode paramNode, PtrParamNode infoNode, 
                       Grid* ptGrid )
: FeSpace( paramNode, infoNode, ptGrid ) {
  isHierarchical_ = true;
}

FeSpaceHi ::~FeSpaceHi() {
}

void FeSpaceHi::MapNodalBCs() {
  LOG_TRACE(feSpaceHi) << "Mapping Nodal BCs";
  StdVector<UInt> actNodes;

  // check if feFunction was defined
  if( !feFunction_ ) {
    EXCEPTION("No FeFunction set at FeSpace");
  }

  // check if feFunction has a result assigned
  if( !feFunction_->GetResultInfo()) {
    EXCEPTION("No resultinfo associated with feFunction");
  }

  //Get Grip of HdBC List for the fefunction
  const HdBcList hdbcs = feFunction_->GetHomDirichletBCs();
  HdBcList::const_iterator actHBC;
  UInt dofsPerUnknown = GetNumDofs();

  for(actHBC = hdbcs.Begin(); actHBC != hdbcs.End(); actHBC++) {
    // Get EntityIterator
    GetNodesOfEntities(actNodes,(*actHBC)->entities);
    for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){
      UInt vNode = actNodes[iNode];

      if( nodeMap_.BcKeys.find(vNode) == nodeMap_.BcKeys.end()){
        //nodeMap_.BcKeys[node] = StdVector<BaseFeFunction::BcType>(dofsPerUnknown,BaseFeFunction::NOBC);
        nodeMap_.BcKeys[vNode] = StdVector<BcType>(dofsPerUnknown);
        nodeMap_.BcKeys[vNode].Init(NOBC);
      }
      LOG_DBG(feSpaceHi) << "\tHDBC for vNode " << vNode << ", dof " <<  (*actHBC)->dof;
      nodeMap_.BcKeys[vNode][(*actHBC)->dof] = HDBC;
      bcCounter_[HDBC]++;
    }
  }
  /* 
   * In a hierarchical space, only the equations corresponding to a vertex
   * get fixed. The edge / face / inner equations are set to zero.
   */
  const IdBcList idbcs = feFunction_->GetInHomDirichletBCs();
  IdBcList::const_iterator actIBC;
  for(actIBC = idbcs.Begin(); actIBC != idbcs.End(); actIBC++) {
    // Get all (Virtual) Nodes of the list
    GetNodesOfEntities(actNodes,(*actIBC)->entities);
    for(UInt iNode = 0 ; iNode < actNodes.GetSize();iNode++){

      UInt vNode = actNodes[iNode];
      //TODO find the source
      //make it zero based
      if( nodeMap_.BcKeys.find(vNode) == nodeMap_.BcKeys.end()){
        nodeMap_.BcKeys[vNode] = StdVector<BcType>(dofsPerUnknown);
        nodeMap_.BcKeys[vNode].Init(NOBC);
      }

      // check, if this node belongs to a vertex
      // yes set "true" IDBC
      //  no: set HDBC
      if( nodesType_[vNode] == BaseFE::VERTEX ) {

        // check first, if this node was already processed
        if( nodeMap_.BcKeys[vNode][(*actIBC)->dof] != IDBC) {
          nodeMap_.BcKeys[vNode][(*actIBC)->dof] = IDBC;
          bcCounter_[IDBC]++;
        }
      } else {
        nodeMap_.BcKeys[vNode][(*actIBC)->dof] = HDBC ;
      }
    } // loop over virtual nodes
  } // loop over idbcs

  //Get Grip of constraint List for the fefunction
  const ConstraintList constraints = feFunction_->GetConstraints();
  ConstraintList::const_iterator actConstr;
  for(actConstr = constraints.Begin(); actConstr != constraints.End(); actConstr++) {
    StdVector<UInt> slaveNodes;
    StdVector<UInt> slaveVertexNodes;
    GetNodesOfEntities(slaveNodes,(*actConstr)->slaveEntities);
    GetNodesOfEntities(slaveVertexNodes,(*actConstr)->slaveEntities, 
                       BaseFE::VERTEX);
    UInt masterDof = (*actConstr)->masterDof;
    UInt slaveDof = (*actConstr)->slaveDof;
    UInt mNode = slaveNodes[0];

    for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
      UInt slaveNode = slaveNodes[iNode];
      if( nodeMap_.BcKeys.find(slaveNode) == nodeMap_.BcKeys.end()){
        nodeMap_.BcKeys[slaveNode] = StdVector<BcType>(dofsPerUnknown);
        nodeMap_.BcKeys[slaveNode].Init(NOBC);
      }

      // check, if this node belongs to a vertex
      // yes set "true" CS
      //  no: set HDBC
      if( slaveVertexNodes.Find( slaveNode ) != -1 ) {
        nodeMap_.BcKeys[slaveNode][slaveDof] = CS;
        nodeMap_.constraintNodes[std::pair<Integer,Integer>(slaveNode,slaveDof)] = 
            std::pair<Integer,Integer>(mNode,masterDof);
        bcCounter_[CS]++; 
      } else {
        nodeMap_.BcKeys[slaveNode][slaveDof] = HDBC;
      }
    } // loop over nodes
  } // loop over constraints
}



void FeSpaceHi::SetRegionOrder( RegionIdType region, 
                                const ApproxOrder& order ) {
  // store order for this region
  regionOrder_[region] = order;


  // loop over all elements of this region
  LOG_DBG(feSpaceHi) << "Assigning edge / face order";
  StdVector<Elem*> elems;
  ptGrid_->GetElems( elems, region );
  UInt numElems = elems.GetSize();
  for( UInt i = 0; i < numElems; ++i ) {

    const Elem& el = *(elems[i]);
    LOG_DBG3(feSpaceHi) << "Treating element #" << el.elemNum;

    // Fetch reference element and set correct order
    FeHi * myFe = GetFeHi( region, el.type);

    // important: do not adjust the entity order here, as the
    // edge / face order is not yet initialized
    SetElemOrder( elems[i], myFe, regionOrder_[region], false );

    // a) loop over all edges
    // -----------------------
    UInt numEdges = el.edges.GetSize();
    LOG_DBG3(feSpaceHi) << "Checking " << numEdges << " edges ";
    const StdVector<UInt>& edgeOrders = myFe->GetEdgeOrder();
    for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ){
      UInt edgeNum = std::abs(el.edges[iEdge]);
      UInt edgeOrder = edgeOrders[iEdge];

      // check if edge got already mapped
      if( orderEdges_.find(edgeNum) == orderEdges_.end() ) {
        orderEdges_[edgeNum] = edgeOrder;
        LOG_DBG3(feSpaceHi) << "\tedge #" << edgeNum << ", order: " 
            << edgeOrder;
      } else {
        UInt oldOrder = orderEdges_[edgeNum];
        LOG_DBG3(feSpaceHi) << "\tedge #" << edgeNum 
            << " -> Already set to " << oldOrder; 
        // check if previous set order is different from current
        if( edgeOrder != oldOrder) {
          orderEdges_[edgeNum] = std::max( edgeOrder, oldOrder );
          LOG_DBG3(feSpaceHi) << "\t\t=> Re-Set to " 
              << orderEdges_[edgeNum]<< " (max. rule)";
          adjustedEdges_.insert( edgeNum );
        }
      }
    } // loop over edges


    // b) loop over all faces (only in 3D case )
    // -----------------------------------------
    if( ptGrid_->GetDim() == 3 ) {
      UInt numFaces = el.faces.GetSize();
      LOG_DBG3(feSpaceHi) << "Checking " << numFaces << " faces ";
      const StdVector<boost::array<UInt,2> >& faceOrders = 
          myFe->GetFaceOrder();
      for( UInt iFace = 0; iFace < numFaces; ++iFace ){
        UInt faceNum = el.faces[iFace];
        boost::array<UInt,2> faceOrder = faceOrders[iFace];

        // check if face got already mapped
        if( orderFaces_.find(faceNum) == orderFaces_.end() ) {
          orderFaces_[faceNum] = faceOrder;
          LOG_DBG3(feSpaceHi) << "\tface #" << faceNum << ", order: " 
              << faceOrder[0] << ", " << faceOrder[1];
        } else {
          boost::array<UInt,2> oldOrder = orderFaces_[faceNum];
          LOG_DBG3(feSpaceHi) << "\tface #" << faceNum 
              << " -> Already set to " << oldOrder[0] << ", " 
              << oldOrder[1];
          // check if previous set order is different from current
          if( faceOrder != oldOrder) {
            orderFaces_[faceNum][0] = std::max( faceOrder[0], oldOrder[0] );
            orderFaces_[faceNum][1] = std::max( faceOrder[1], oldOrder[1] );
            LOG_DBG3(feSpaceHi) << "\t\t=> Re-Set to " 
                << orderFaces_[faceNum][0] << ", "
                << orderFaces_[faceNum][1] << " (max. rule)";
            adjustedFaces_.insert( faceNum );
          }
        }
      }  // loop over faces
    } // check for 3D case 

  } // loop over elements
}

void FeSpaceHi::SetElemOrder( const Elem* ptEl, FeHi* ptFe, 
                              const ApproxOrder& order ,
                              bool applyMaxRule ) {
  
  LOG_DBG(feSpaceHi) << "In SetElemOrder for elem " << ptEl->elemNum 
                     << ", maxRule: " << applyMaxRule;
  
  // Stage 1: Set order as given from the region template.
  // This order template specified either an isotropic polynomial degree
  // or a direction dependent order (for each unknown).

  // check for isotropy -> This has to be changed to a more
  // general representation soon
  if( order.IsIsotropic() ) {
    ptFe->SetIsoOrder( order.GetIsoOrder() );
  } else {

    // check, if element is a surface element
    if( Elem::shapes[ptEl->type].dim < ptGrid_->GetDim() ) {
      LOG_DBG3(feSpaceHi) << "\t-> Mapping volume dirs to surface";
      
      // map element to volume element
      const Elem* ptVolEl = (dynamic_cast<const SurfElem*>(ptEl))->ptVolElems[0];
      LOG_DBG3(feSpaceHi) << "\tVol element: " << ptVolEl->elemNum;
      assert(ptVolEl);
      // we have surface element 
      StdVector<UInt> surfLocDir;
      shared_ptr<ElemShapeMap> esm = ptGrid_->GetElemShapeMap(ptVolEl, false);

      esm->MapSurfLocDirs( ptEl, surfLocDir );
      LOG_DBG3(feSpaceHi) << "\tSurface directions w.r.t. volume: " 
                          << surfLocDir.Serialize(); 
     

      const StdVector<UInt> & volMaxLocDir  = order.GetMaxOrderLocDir();
      StdVector<UInt> surfOrder(surfLocDir.GetSize());
      for (UInt i = 0; i < surfLocDir.GetSize(); ++i ) {
        surfOrder[i] = volMaxLocDir[surfLocDir[i]]; 
      }
      LOG_DBG3(feSpaceHi) << "\tAnisotropic order w.r.t. to vol dir: " << 
          volMaxLocDir.Serialize();
      LOG_DBG3(feSpaceHi) << "\tAnisotropic order w.r.t. to surface dir: " << 
          surfOrder.Serialize();
      ptFe->SetAnisoOrder(surfOrder, ptEl );
    } else {

      // we have a volume element
      // set maximum order for all dofs per local direction
      ptFe->SetAnisoOrder(order.GetMaxOrderLocDir(), ptEl );
    }
  }

  // Stage 2: Adjust edge / face order according to minimum /
  // maximum rule, in case the element has neighboring elements with
  // different order. We will use the maximum rule
  if( applyMaxRule && orderEdges_.size() ) {
    // loop over all edges
    UInt numEdges = ptEl->edges.GetSize();
    for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ) {
      UInt edgeNum = std::abs( ptEl->edges[iEdge] );
      // check if edge got adjusted
      if( orderEdges_.find(edgeNum) != orderEdges_.end() ) {
        LOG_DBG3(feSpaceHi) << "Setting edge " << edgeNum
            << " to order " << orderEdges_[edgeNum];
        ptFe->SetEdgeOrder( iEdge, orderEdges_[edgeNum] );
      }
    }

  }

  if( applyMaxRule && orderFaces_.size() 
      && ptGrid_->GetDim() == 3  ) {
    // loop over all faces
    UInt numFaces = ptEl->faces.GetSize();
    for( UInt iFace = 0; iFace < numFaces; ++iFace ) {
      UInt faceNum = ptEl->faces[iFace];
      // check if face got adjusted
      if( orderFaces_.find(faceNum) != orderFaces_.end() ) {
        LOG_DBG3(feSpaceHi) << "Setting face " << faceNum
            << " to order " << orderFaces_[faceNum][0]
                                                    << ", " << orderFaces_[faceNum][1];
        ptFe->SetFaceOrder( iFace, orderFaces_[faceNum] );
      }
    }

  }
}

void FeSpaceHi::AdjustEntityOrder() {
  LOG_TRACE(feSpaceHi) << "Adjusting entity order";

  // This method is called in the end.

  // The idea is as follows:
  // Up to now the polynomial orders of all edges and faces ar stored
  // explicitly. However, we need to store only those, which
  // have been adjusted explicitly due to the maximum rule.
  // Thus we store only explicitly the edge and face order of those
  // edges / faces adjusted.

  // Create temporary map for edges
  boost::unordered_map<UInt, UInt> edgeOrders;

  // Loop over all adjusted edges
  LOG_DBG(feSpaceHi) << "Adjusting edge order";
  boost::unordered_set<UInt>::const_iterator edgeIt = adjustedEdges_.begin();
  for( ; edgeIt != adjustedEdges_.end(); ++edgeIt ) {
    LOG_DBG3(feSpaceHi) << "\tOrder of edge #" << *edgeIt 
        << " is " <<  orderEdges_[*edgeIt];
    edgeOrders[*edgeIt] = orderEdges_[*edgeIt];
  } // loop over edges

  // store only adjusted edges back
  orderEdges_ = edgeOrders;

  // treat faces only in 3D
  if( ptGrid_->GetDim() < 3) return;

  // Create temporary map for faces
  boost::unordered_map<UInt, boost::array<UInt,2> > faceOrders;

  // Loop over all adjusted faces
  LOG_DBG(feSpaceHi) << "Adjusting face order";
  boost::unordered_set<UInt>::const_iterator faceIt = adjustedFaces_.begin();
  for( ; faceIt != adjustedFaces_.end(); ++faceIt ) {
    LOG_DBG3(feSpaceHi) << "\tOrder of face #" << *faceIt 
        << " is " <<  orderFaces_[*faceIt][0]
                                           << ", " << orderFaces_[*faceIt][1];
    faceOrders[*faceIt] = orderFaces_[*faceIt];
  } // loop over faces

  // store only adjusted faces back
  orderFaces_ = faceOrders;
}

void FeSpaceHi::FixHigherOrderAnisoDofs() {
  //    std::cerr 
  //    << "===================================================\n"
  //    << " F I X I I N G   A N I S O T R O P I C    D O F S  \n" 
  //    << "===================================================\n";


  
  
  StdVector< shared_ptr<EntityList> > fctEntList = 
      feFunction_->GetEntityList();
  EntityList::ListType actListType = EntityList::NO_LIST;
  UInt numDofs = GetNumDofs();

  // Leave, if we only have a scalar valued problem
  if( numDofs == 1 ) 
    return;
    
  // loop over all entitylists (i.e. regions)
  for(UInt actList = 0;actList <  fctEntList.GetSize(); actList++){
    actListType = fctEntList[actList]->GetType();

    // only handle elements
    if ( ! (actListType == EntityList::ELEM_LIST) &&
        ! (actListType == EntityList::NC_ELEM_LIST))  {
      LOG_DBG(feSpaceHi) << "\tLEAVING";
      continue;
    }
    // create entity list for given region
    EntityList* actElemList = fctEntList[actList].get();
    EntityIterator entIt = actElemList->GetIterator();

    // check if region is anisotropic and get order array
    std::map<RegionIdType,ApproxOrder>::iterator approxIt =  
        regionOrder_.find( entIt.GetElem()->regionId );
    ApproxOrder & ord = approxIt->second;

    // loop over all elements
    for(entIt.Begin(); !entIt.IsEnd(); entIt++ ){
      const Elem * ptElem = entIt.GetElem();

      // Fetch reference element and set correct order
      FeHi * ptFe = GetFeHi( ptElem->regionId, ptElem->type);
      SetElemOrder( ptElem, ptFe, ord, true );

      // get vector of virtual nodes for current element
      StdVector<UInt> nodes;
      GetNodesOfElement( nodes, ptElem, BaseFE::ALL);

      // loop over all dofs
      for( UInt iDof = 0; iDof < numDofs; ++iDof ) {

        // slice out vector of order for given dof
        StdVector<UInt> dofOrder = ord.GetDofOrder( iDof );

        // get indices, which exceed the global dof order
        std::set<UInt> indices; 
        ptFe->GetNodesExceedingOrder( indices, dofOrder, ptElem );

        // loop over all "exceeding" virtual nodes and fix their
        // value by a homogeneous BC
        std::set<UInt>::const_iterator indexIt = indices.begin(); 
        for( ; indexIt != indices.end(); ++indexIt ) {
          UInt vNode = nodes[*indexIt];
          if( nodeMap_.BcKeys.find(vNode) == nodeMap_.BcKeys.end()){
            nodeMap_.BcKeys[vNode] = StdVector<BcType>(numDofs);
            nodeMap_.BcKeys[vNode].Init(FeSpace::NOBC);
          }
          nodeMap_.BcKeys[vNode][iDof] = HDBC;
          bcCounter_[HDBC]++;

        } // loop over exceeding virtual nodes
      } // loop over dofs
    } // loop over all elements
  } // loop over regions
}

} // end of namespace
