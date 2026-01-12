#include "FeSpaceHi.hh"

#include "FeHi.hh"



#include "DataInOut/Logging/LogConfigurator.hh"

#include "OLAS/algsys/AlgebraicSys.hh"
#include "Driver/Assemble.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "PDE/SinglePDE.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

#include "def_use_pardiso.hh"
#include "def_use_suitesparse.hh"
#include "def_use_superlu.hh"

namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(feSpaceHi, "feSpaceHi")

  // -----------------------------------
  FeSpaceHi::MapContext::MapContext() {
    algSys = NULL;
    assemble = NULL;
    sol = NULL;
  }
  
  FeSpaceHi::MapContext::~MapContext() {
    delete algSys;
    algSys = NULL;
    
    delete assemble;
    assemble = NULL;
    
    delete sol;
    sol = NULL;
    
    entityEqns.Clear();
    
  }

FeSpaceHi ::FeSpaceHi (PtrParamNode paramNode, PtrParamNode infoNode, 
                       Grid* ptGrid )
: FeSpace( paramNode, infoNode, ptGrid ) {
  isHierarchical_ = true;
}

FeSpaceHi ::~FeSpaceHi() {
    // delete all mapping contexts
    std::map<std::string, MapContext*>::iterator it = entityCtx_.begin(); 
    for( ; it != entityCtx_.end(); ++it ) {
      delete it->second;
    }
    entityCtx_.clear();
}

//! \copydoc FeSpace::MapCoefFctToSpace
void FeSpaceHi ::MapCoefFctToSpace( StdVector<shared_ptr<EntityList> > support,
                        shared_ptr<CoefFunction> coefFct,
                        shared_ptr<BaseFeFunction> feFct,
                        std::map<Integer, Double>& vals, bool cache,
                        const std::set<UInt>& comp,
                        bool updatedGeo) {
  MapCoefFctToSpacePriv<Double>( support, coefFct, feFct, vals, cache, comp );
}

//! \copydoc FeSpace::MapCoefFctToSpace
void FeSpaceHi ::MapCoefFctToSpace( StdVector<shared_ptr<EntityList> > support,
                        shared_ptr<CoefFunction> coefFct,
                        shared_ptr<BaseFeFunction> feFct,
                        std::map<Integer, Complex>& vals, bool cache,
                        const std::set<UInt>& comp,
                        bool updatedGeo) {
  MapCoefFctToSpacePriv<Complex>( support, coefFct, feFct, vals, cache, comp );
}



void FeSpaceHi::MapNodalBCs() {
  LOG_TRACE(feSpaceHi) << "Mapping Nodal BCs";
  StdVector<UInt> actNodes;

  // check if feFunction was defined
  shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
  if( !feFct ) {
    EXCEPTION("No FeFunction set at FeSpace");
  }

  // check if feFunction has a result assigned
  if( !feFct->GetResultInfo()) {
    EXCEPTION("No resultinfo associated with feFunction");
  }

  //Get Grip of HdBC List for the fefunction
  const HdBcList hdbcs = feFct->GetHomDirichletBCs();
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
      
      // loop over all dofs
      std::set<UInt>::const_iterator dofIt = (*actHBC)->dofs.begin();
      for( ; dofIt != (*actHBC)->dofs.end(); ++dofIt) { 
        LOG_DBG(feSpaceHi) << "\tHDBC for vNode " << vNode << ", dof " <<  *dofIt;
        nodeMap_.BcKeys[vNode][*dofIt] = HDBC;
        bcCounter_[HDBC]++;
      } // dofs
    }
  }
  //Get Grip of IdBC List for the fefunction
  const IdBcList idbcs = feFct->GetInHomDirichletBCs();
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
      // check first, if this node was already processed
      // loop over all dofs
      std::set<UInt>::const_iterator dofIt = (*actIBC)->dofs.begin();
      for( ; dofIt != (*actIBC)->dofs.end(); ++dofIt) { 
        // check first, if this node was already processed
        if( nodeMap_.BcKeys[vNode][*dofIt] != IDBC) {
          nodeMap_.BcKeys[vNode][*dofIt] = IDBC;
          bcCounter_[IDBC]++;
        }
      } // dofs
    } // loop over virtual nodes
  } // loop over idbcs

  //Get Grip of constraint List for the fefunction
  const ConstraintList constraints = feFct->GetConstraints();
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
    UInt numEdges = el.extended->edges.GetSize();
    LOG_DBG3(feSpaceHi) << "Checking " << numEdges << " edges ";
    const StdVector<UInt>& edgeOrders = myFe->GetEdgeOrder();
    for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ){
      UInt edgeNum = std::abs(el.extended->edges[iEdge]);
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
      UInt numFaces = el.extended->faces.GetSize();
      LOG_DBG3(feSpaceHi) << "Checking " << numFaces << " faces ";
      const StdVector<boost::array<UInt,2> >& faceOrders = 
          myFe->GetFaceOrder();
      for( UInt iFace = 0; iFace < numFaces; ++iFace ){
        UInt faceNum = el.extended->faces[iFace];
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
    UInt numEdges = ptEl->extended->edges.GetSize();
    for( UInt iEdge = 0; iEdge < numEdges; ++iEdge ) {
      UInt edgeNum = std::abs( ptEl->extended->edges[iEdge] );
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
    UInt numFaces = ptEl->extended->faces.GetSize();
    for( UInt iFace = 0; iFace < numFaces; ++iFace ) {
      UInt faceNum = ptEl->extended->faces[iFace];
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
  // Up to now the polynomial orders of all edges and faces are stored
  // explicitly. However, we need to store only those, which
  // have been adjusted explicitly due to the maximum rule.
  // Thus we store only explicitly the edge and face order of those
  // edges / faces adjusted.

  // Create temporary map for edges
  std::unordered_map<UInt, UInt> edgeOrders;

  // Loop over all adjusted edges
  LOG_DBG(feSpaceHi) << "Adjusting edge order";
  std::unordered_set<UInt>::const_iterator edgeIt = adjustedEdges_.begin();
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
  std::unordered_map<UInt, boost::array<UInt,2> > faceOrders;

  // Loop over all adjusted faces
  LOG_DBG(feSpaceHi) << "Adjusting face order";
  std::unordered_set<UInt>::const_iterator faceIt = adjustedFaces_.begin();
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


  
  shared_ptr<BaseFeFunction> feFct = feFunction_.lock(); // request a strong pointer
  assert(feFct);
  StdVector< shared_ptr<EntityList> > fctEntList = 
      feFct->GetEntityList();
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


    // loop over all elements
    for(entIt.Begin(); !entIt.IsEnd(); entIt++ ){
      const Elem * ptElem = entIt.GetElem();

      // Fetch reference element and set correct order
      RegionIdType eRegion = GetVolElem(ptElem)->regionId;
      FeHi * ptFe = GetFeHi( eRegion, ptElem->type);

      // Get order approximation
      std::map<RegionIdType,ApproxOrder>::iterator approxIt =  
          regionOrder_.find( eRegion );
      if( approxIt == regionOrder_.end() ) {
        EXCEPTION( "Could not find approximation for element #"
            << ptElem->elemNum << " of region " 
            << ptGrid_->GetRegion().ToString(ptElem->regionId) );
      }
      
      ApproxOrder & ord = approxIt->second;
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

template<typename T>
void FeSpaceHi::MapCoefFctToSpacePriv(StdVector<shared_ptr<EntityList> > entityLists,
                                          shared_ptr<CoefFunction> coefFct,
                                          shared_ptr<BaseFeFunction> feFct,
                                          std::map <Integer, T>& vals,
                                          bool cache,
                                          const std::set<UInt>& comp ) {


  // Perform some checks at first:
  // a) if coefficient function is constant -> perform "easy mapping", without
  //    any caching. Here we simply take the constant value and
  // b) if coefficient function has a general (spatial) dependency,
  //    we create solve FE problem on the subset of the FeSpace, as defined
  //    by the boundary condition, i.e. a mass-bilinearform is created and an
  //    according RHS integrator. The auxilliary data needed can be cached for
  //    repeated access / changing values (e.g. time / frequency dependend boundary+
  //    conditions.

  // assemble complete name
  std::string entityNames;
  for( UInt i = 0; i < entityLists.GetSize()-1; ++i ) {
    entityNames += entityLists[i]->GetName() + ", ";
  }
  entityNames += entityLists.Last()->GetName();
  
  if (IS_LOG_ENABLED(feSpaceHi, trace)) {
    StdVector<UInt> compVec;
    std::set<UInt>::const_iterator it = comp.begin();
    for (; it != comp.end(); ++it ) {
      compVec.Push_back(*it);
    }
    
    LOG_TRACE(feSpaceHi) << "Mapping coeffct " << coefFct->ToString()
                            << " on " << entityNames << " for dofs "
                            << compVec.ToString() << " for FeFunction "
                            << SolutionTypeEnum.ToString(feFct->GetResultInfo()->resultType);
  }

  if( (coefFct->GetDependency() == CoefFunction::CONSTANT ||
      coefFct->GetDependency() == CoefFunction::TIMEFREQ) && 
      type_ != HCURL) {
    // --------------------------
    //  SIMPLE MAPPING MECHANISM
    // --------------------------
    LOG_DBG(feSpaceHi) << "=> Applying simple mapping algorithm for constant values";
    LocPointMapped lpm;
    StdVector<Integer> eqns, vEqns;
    Vector<T> dummyVec;
    if( coefFct->GetDimType() == CoefFunction::SCALAR ) {
      dummyVec.Resize(1);
      coefFct->GetScalar(dummyVec[0], lpm);
    } else {
      coefFct->GetVector( dummyVec, lpm);
    }

    // loop over all lists 
    for( UInt i = 0; i < entityLists.GetSize(); ++i ) {
      shared_ptr<EntityList> actList = entityLists[i];

      // loop over all dofs
      std::set<UInt>::const_iterator it = comp.begin();
      for( ; it != comp.end(); ++it) {

        T val = dummyVec[*it];
        this->GetEntityListEqns( eqns, actList, *it );
        UInt numEqns = eqns.GetSize();

        // switch depending on space type:
        if( !this->IsHierarchical() ) {
          // --- non-hierachical case ---
          for( UInt i = 0; i < numEqns; ++i ) {
            vals[eqns[i]] = val;
          }
        } else {
          // -- hierachical case ---
          for( UInt i = 0; i < numEqns; ++i ) {
            vals[eqns[i]] = 0.0;
          }

          // get nodal vertex nodes
          this->GetEntityListEqns(vEqns, actList, *it, BaseFE::VERTEX );
          UInt numVEqns = vEqns.GetSize();
          for( UInt i = 0; i < numVEqns; ++i ) {
            vals[vEqns[i]] = val;
          }
        } // loop: dofs
      } // loop: lists
    } // loop: dofs
  } else
  {
    // ---------------------------
    //  GENERAL MAPPING MECHANISM
    // ---------------------------
    LOG_DBG(feSpaceHi) << "=> Applying general mapping algorithm for general values";
    std::map<Integer,Integer> eqnMap;
    std::set<Integer> registeredEqns;
    std::string name = entityNames;
    MapContext* ctx = NULL;

    // check, if caching is activated and if entry can be found in cache
    if( cache == true
        && entityCtx_.find(name) != entityCtx_.end() ) {
      ctx = entityCtx_[name];
    } else {
      // generate new context
      ctx = new MapContext();

      bool isComplex = std::is_same<T,Complex>::value;

      // generate new algebraic system
      ctx->olasNode.reset(new ParamNode());
      ctx->olasNode->SetName(std::string("IDBC-") + name);

      // Explicitly define optimized solver, if available.
      PtrParamNode solverListNode(new ParamNode());
      solverListNode->SetName("solverList");

      PtrParamNode solverNode(new ParamNode());
#if defined(USE_PARDISO)
      solverNode->SetName("pardiso");
      PtrParamNode statsNode(new ParamNode());
      statsNode->SetName("stats");
      statsNode->SetValue("no");
      solverNode->AddChildNode(statsNode);      
      PtrParamNode loggingNode(new ParamNode());
      loggingNode->SetName("logging");
      loggingNode->SetValue("no");
      solverNode->AddChildNode(loggingNode);      
#elif defined(USE_SUITESPARSE)
      solverNode->SetName("cholmod");
#elif defined(USE_SUPERLU)
      solverNode->SetName("superlu");
#else
      solverNode->SetName("directLDL");
#endif

      solverListNode->AddChildNode(solverNode);

      PtrParamNode solverIdNode(new ParamNode());
      solverIdNode->SetName("id");
      solverIdNode->SetValue("default");
      solverNode->AddChildNode(solverIdNode);

      ctx->olasNode->AddChildNode(solverListNode);

      ctx->infoNode.reset(new ParamNode(ParamNode::INSERT));
      ctx->algSys = new AlgebraicSys( ctx->olasNode, ctx->infoNode, isComplex );

      // generate new assemble class with dummy info points
      PtrParamNode infoNode = 
          PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));
      BasePDE::AnalysisType aType =
          isComplex ? BasePDE::HARMONIC : BasePDE::STATIC;
             
      MathParser * mp = feFct->GetPDE()->GetDomain()->GetMathParser();
      ctx->assemble = new Assemble( ctx->algSys, aType, mp, infoNode );

      // --------------------------------------------------------------------
      // generate (dimensional and space-dependent) interpolation (bi)linear
      // forms for the auxiliary problem
      // --------------------------------------------------------------------
      for( UInt i = 0; i < entityLists.GetSize(); ++i ) {
        shared_ptr<EntityList> actList = entityLists[i];
        // Only proceed, if list type is elem list or surf elem list
        if( ! ( actList->GetType() == EntityList::ELEM_LIST ||
            actList->GetType() == EntityList::SURF_ELEM_LIST ) ) {
          continue; 
        }

        const Elem* el = actList->GetIterator().GetElem();
        UInt dim = Elem::shapes[el->type].dim;
        UInt dofDim = this->GetNumDofs();
        BiLinearForm *massInt = feFct->GenerateInterpolBilinForm(dim, dofDim, true);
        LinearForm * rhsInt = feFct->GenerateInterpolLinForm(dim, dofDim, coefFct, true);

        BiLinFormContext * massCtx = new BiLinFormContext( massInt, STIFFNESS);
        massInt->SetName("Interpolator");
        massCtx->SetEntities( actList, actList );
        massCtx->SetFeFunctions(feFct, feFct);
        ctx->assemble->AddBiLinearForm(massCtx);

        LinearFormContext * rhsCtx = new LinearFormContext( rhsInt );
        rhsInt->SetName("RhsInterpolator");
        rhsCtx->SetEntities( actList );
        rhsCtx->SetFeFunction(feFct);
        ctx->assemble->AddLinearForm(rhsCtx);

        // generate a mapping global eqnNr -> entityList local one for all dofs
        // loop over all dofs
        if( comp.size() == 0 ) {
          // scalar problem
          this->GetEntityListEqns( ctx->entityEqns, actList );
        } else {
          StdVector<Integer> tmp;
          std::set<UInt>::const_iterator it = comp.begin();
          for( ; it != comp.end(); it++ ) {
            this->GetEntityListEqns( tmp, actList, *it );
            for( UInt k = 0; k < tmp.GetSize(); ++k ) {
              ctx->entityEqns.Push_back( tmp[k] ) ;
            }
          }
        }
        // fill map and pass it to assemble class
        UInt numEqns = ctx->entityEqns.GetSize();
        for( UInt i = 0; i < numEqns; ++i ) {
          eqnMap[ctx->entityEqns[i]] = i+1;
          registeredEqns.insert(i+1);
        }
      } // loop: entityLists
      // ------------------------
      //  setup algebraic system
      // ------------------------
      // 1) setup matrix graph
      std::map<FeFctIdType, FeFctIdType> feFctIdMap;
      ctx->algSys->GraphSetupInit( 1,  false );
      FeFctIdType newFctId = ctx->algSys->ObtainFctId("interpolation");
      ctx->algSys->RegisterFct(newFctId, eqnMap.size(), eqnMap.size() );
      feFctIdMap[newFctId] = feFct->GetFctId();
      ctx->assemble->SetEqnCustomMap(eqnMap, feFctIdMap);

      // define matrix graph and SBM blocks
      AlgebraicSys::SBMBlockDef sbmBlock;
      sbmBlock[newFctId] = registeredEqns;
      ctx->algSys->DefineSBMMatrixBlock( sbmBlock, false );
      ctx->algSys->FinishRegistration();
      ctx->assemble->SetupMatrixGraph(newFctId, newFctId);
      ctx->algSys->GraphSetupDone();
      ctx->algSys->CreateLinSys();
      ctx->algSys->InitMatrix();

      // create solver and preconditioner
      ctx->algSys->CreateSolver();
      ctx->algSys->CreatePrecond();

      // now reset AlgebraicSystem
      ctx->algSys->InitRHS();
      ctx->algSys->InitSol();

      // assemble mapping matrix
      ctx->assemble->AssembleMatrices();

      // setup the preconditioner and solver
      ctx->algSys->SetupPrecond(); // analysis_id ctx->infoNode);
      ctx->algSys->SetupSolver(); // analysis_id ctx->infoNode);

      // initialize solution SBM vector
      ctx->sol = new SBM_Vector();
      ctx->sol->Resize(1);
      Vector<T> * tmp = new Vector<T>(eqnMap.size());
      ctx->sol->SetSubVector(tmp,0);
      ctx->sol->SetOwnership(true);


      // Store the context in the map
      if( cache ) {
        entityCtx_[name] = ctx;
      }
    } // end of first time setup of matrix

    // update the RHS due to the new coefficient vector
    ctx->algSys->InitRHS();
    ctx->assemble->AssembleLinRHS();


    // solve system and aquire solution
    ctx->algSys->Solve(); // analysis_id ctx->infoNode);
    ctx->algSys->GetSolutionVal(*(ctx->sol));

    Vector <T> & sol = dynamic_cast<Vector<T> &>(*(ctx->sol->GetPointer(0)));

    // store result values to vals-map
    UInt numEqns = ctx->entityEqns.GetSize();

    for( UInt i = 0; i < numEqns; ++i ) {
      if( ctx->entityEqns[i] != 0 ) {
        vals[ctx->entityEqns[i]] = sol[i];
      }
    }

    // Print out information about the solution of the system
    // ctx->infoNode->ToXML(std::cerr);

    if (!cache){ // Delete ctx before it leaves scope if it is not stored
      delete ctx;
    }

  } // end of general mapping section
}


} // end of namespace

template void FeSpaceHi::MapCoefFctToSpacePriv<Double>( StdVector<shared_ptr<EntityList> > ,
    shared_ptr<CoefFunction>,
    shared_ptr<BaseFeFunction> feFct,
    std::map <Integer, Double>&,
    bool,const std::set<UInt>&);
template void FeSpaceHi::MapCoefFctToSpacePriv<Complex>( StdVector<shared_ptr<EntityList> > ,
    shared_ptr<CoefFunction>,
    shared_ptr<BaseFeFunction> feFct,
    std::map <Integer, Complex>&,
    bool,const std::set<UInt>&);
