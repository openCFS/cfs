// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "eqnMap.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh" 
#include "Utils/coordSystem.hh"
#include "DataInOut/Logging/cfslog.hh"


namespace CoupledField {

  // declare class specific logging stream
  DECLARE_LOG(eqnMap)
  DEFINE_LOG(eqnMap, "eqnMap")


    EqnMap::EqnMap(Grid* grid, PdeIdType pdeId, bool sortEqns ) {
    ENTER_FCN( "EqnMap::EqnMap" );

    ptGrid_ = grid;

    sortEqns_ = sortEqns;
    pdeId_ = pdeId;
    isFinalized_ = false;


    numEqns_ = 0;
    numRealEqns_ = 0;
    numIdBcs_ = 0;
    numCs_ = 0;
    
    numLocNodes_ = 0;
    numLocElems_ = 0;
    numLocEdges_ = 0;
    numLocFaces_ = 0;

  }

  EqnMap::~EqnMap() {
    ENTER_FCN( "EqnMap::~EqnMap" );

  }

  //! ======================================================================
  //! SET AND INITIALIZATION METHODS
  //! ======================================================================
  void EqnMap::AddResult( ResultInfo& result, 
		          shared_ptr<EntityList> list ) {
    ENTER_FCN( "EqnMap::AddResult" );
    resEntMap_[result].Push_back( list );
  }

  void EqnMap::SetHomoDirichletBCs(HdBcList& hdBcs) {
    ENTER_FCN( "EqnMap::SetHomoDirichletBCs" );

    for( UInt i = 0; i< hdBcs.GetSize(); i++ ) {
      ResultInfo actResult = *(hdBcs[i]->result);
      hdBcs_[actResult].Push_back(hdBcs[i]);
    }
  }
  
  void EqnMap::SetInhomDirichletBCs( IdBcList& idBcs ) {
    ENTER_FCN( "EqnMap::SetInhomDirichletBCs" );

    for( UInt i = 0; i < idBcs.GetSize(); i++ ) {
      ResultInfo actResult = *(idBcs[i]->result);
      idBcs_[actResult].Push_back( idBcs[i] );
    }
  }
  
  void EqnMap::SetConstraints( ConstraintList& constraints ) {
    ENTER_FCN( "EqnMap::SetConstraints" );
    
    for( UInt i = 0; i < constraints.GetSize(); i++ ) {
      ResultInfo actResult = *(constraints[i]->result);
      constraints_[actResult].Push_back( constraints[i] );
    }
  }
  
  void EqnMap::Finalize() {
    ENTER_FCN( "EqnMap::Finalize" );

    LOG_TRACE(eqnMap) << "Starting Initialization";
    
    // 1) Get all entity lists, which consist of elements. 
    //    They are used lateron for global/local element/nodal number mapping    
    ResultEntityMap::iterator it;
    shared_ptr<ElemList> elemList;
    // iterate over all results
    for (it=resEntMap_.begin(); it!=resEntMap_.end(); it++ ) {
      StdVector<shared_ptr<EntityList> > & lists = (*it).second;
    
      // iterate over entity-lists of this result type
      for (UInt iList=0; iList<lists.GetSize(); iList++) {
	if (lists[iList]->GetType() == EntityList::ELEM_LIST ) {

	  // Add element list to list of mapped regions
	  elemList = dynamic_pointer_cast<ElemList, EntityList>(lists[iList]);
	  locElems_.Push_back( elemList );
	}
      }
    }

    // 2.) Get all entity lists, for which edges/surfaces have to be mapped
    // Get information about which region need mapping of
    // edges / surface 

    // ... to be implemented ...

    // 3.) Perform local<->global mapping of elements
    LOG_TRACE(eqnMap) << "Performing NodeElem Mapping";
    CalcNodeElemMapping();

    // 4.) Perform global->local mapping of edges / surfaces

    // ... to be implemented ...
    
    // ==== PHASE 1: Number all free equations ====
    //

    // Now iterate over all types of maps and assign equation
    // numbers; 

    //   a) node <-> eqnNr
    //   -----------------

    // iterate over all resultDofs
    for ( it=resEntMap_.begin(); it!=resEntMap_.end(); it++ ) {

      // check if resultDof is mapped onto nodes
      if( it->first.definedOn == ResultInfo::NODE ||
	  it->first.definedOn == ResultInfo::PFEM ) {
	for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
	  nodeMappedList_[it->first].Push_back( it->second[iList] );
	}
      }
    }
    
    // assign equation numbers to nodes
    LOG_TRACE(eqnMap) << "Calculating Nodal equations";
    CalcNodalEquations( 1 );

    //   b) edge <-> eqnNr
    //   -----------------

    // iterate over all resultDofs
    for ( it=resEntMap_.begin(); it!=resEntMap_.end(); it++ ) {
      // check if resultDof is mapped onto nodes
      if( it->first.definedOn == ResultInfo::EDGE ||
	  it->first.definedOn == ResultInfo::PFEM ) {
	for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
	  edgeMappedList_[it->first].Push_back( it->second[iList] );
	}
      }
    }

    // Only calc global->local mapping of edges, if at
    // least on entry is present in edgeMappedList_
    if( edgeMappedList_.size() > 0 ) {
      
      // Trigger calculation of edges at grid class
      ptGrid_->MapEdges();

      // calc local<->global mapping of edges 
      CalcEdgeMapping();

      // calc edge / surface equations
      CalcEdgeEquations(1);
    }

    //   c) face <-> eqnNr (only in 3D)
    //   -------------------------------
    
    if( ptGrid_->GetDim() == 3 ) {
      // iterate over all resultDofs
      for ( it=resEntMap_.begin(); it!=resEntMap_.end(); it++ ) {
        // check if resultDof is mapped onto nodes
        if( it->first.definedOn == ResultInfo::FACE ||
            it->first.definedOn == ResultInfo::PFEM ) {
          for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
            faceMappedList_[it->first].Push_back( it->second[iList] );
          }
        }
      }
    }

    // Only calc global->local mapping of faces, if at
    // least on entry is present in faceMappedList_
    if( faceMappedList_.size() > 0 ) {
    
      // Trigger calculation of faces at grid class
      ptGrid_->MapFaces();

      // calc local<->global mapping of faces
      CalcFaceMapping();

      // calc face equations
      CalcFaceEquations(1);
    }



    // d) elem <-> eqnNr ( 'bubble functions' )
    // ----------------------------------------
    
    // iterate over all resultDofs
    for( it = resEntMap_.begin(); it != resEntMap_.end(); it++ ) {
      // check if resultDof is apprxomiated using PFEM
      if( it->first.definedOn == ResultInfo::PFEM ) {
        for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
          elemIntMappedList_[it->first].Push_back( it->second[iList] );
        }
      }
    }
    
    CalcElemInteriorEquations(1);

    //   e) elem <-> eqnNr (only for constants)
    //   -----------------

    // iterate over all resultDofs
    for ( it=resEntMap_.begin(); it!=resEntMap_.end(); it++ ) {

      // check if resultDof is mapped onto nodes
      if( it->first.definedOn == ResultInfo::ELEMENT ) {
	for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
	  elemConstMappedList_[it->first].Push_back( it->second[iList] );
	}
      }
    }
    
    // assign equation numbers to nodes
    CalcElemConstEquations( 1 );

    // ==== PHASE 2: Renumber fixed-equations ====

    // Afterwards re-iterate and give dirichletDOFs highest
    // numbers
    CalcNodalEquations( 2 );
    
    CalcElemConstEquations( 2 );

    // Calc number of 'real equations'
    numRealEqns_ = numEqns_ - numIdBcs_ - numCs_;
    LOG_DBG(eqnMap) << "#equations: " << numEqns_;
    LOG_DBG(eqnMap) << "#dirichletBcs: " << numIdBcs_;
    LOG_DBG(eqnMap) << "#constraints: " << numCs_;
    LOG_DBG(eqnMap) << "#realEquations: " << numRealEqns_;

    
    // Now class is finalized
    isFinalized_ = true;

    LOG_TRACE(eqnMap) << "Finished Initialization\n";
  }

  void EqnMap::ReorderMapping( Integer **order ) {
    ENTER_FCN( "EqnMap::ReorderMapping" );

    
    // Check if any reordering array was given
    if( (*order) == NULL ) {
      //std::cerr << "Performing no reordering";
      return;
    }

    // Iterate over all nodal maps and map the new ordering
    EqnMapType::iterator it1;

    // iterate over all nodal resuls
    for ( it1=nodeEqns_.begin(); it1!=nodeEqns_.end(); it1++ ) {
      Matrix<Integer> & actMap = it1->second;

      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSizeRow(); i++ ) {
	for ( UInt j = 0; j < actMap.GetSizeCol(); j++ ) {
	  if ( actMap[i][j] > 0 ) {
	    actMap[i][j] = (*order)[actMap[i][j]-1];
	  }
	  else if(actMap[i][j] < 0 ) {
	    //due to constraints
	    actMap[i][j] = -(*order)[-actMap[i][j]-1];   
	  }
	}
      }
    }
    
    VecEqnMapType::iterator it2;

    // iterate over all face
    for ( it2=faceEqns_.begin(); it2!=faceEqns_.end(); it2++ ) {
      StdVector<Vector <Integer> >& actMap = it2->second;
      
      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSize(); i++ ) {
	for ( UInt j = 0; j < actMap[i].GetSize(); j++ ) {
	  if ( actMap[i][j] > 0 ) {
	    actMap[i][j] = (*order)[actMap[i][j]-1];
	  }
	  else if(actMap[i][j] < 0 ) {
	    //due to constraints
	    actMap[i][j] = -(*order)[-actMap[i][j]-1];   
	  }
	}
      }
      
    }

    // iterate over all edge resuls
    for ( it2=edgeEqns_.begin(); it2!=edgeEqns_.end(); it2++ ) {
      StdVector<Vector <Integer> >& actMap = it2->second;
      
      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSize(); i++ ) {
	for ( UInt j = 0; j < actMap[i].GetSize(); j++ ) {
	  if ( actMap[i][j] > 0 ) {
	    actMap[i][j] = (*order)[actMap[i][j]-1];
	  }
	  else if(actMap[i][j] < 0 ) {
	    //due to constraints
	    actMap[i][j] = -(*order)[-actMap[i][j]-1];   
	  }
	}
      }
      
    }
    
    
    // Iterate over all element maps and map the new ordering
    for ( it2=elemEqns_.begin(); it2!=elemEqns_.end(); it2++ ) {
      StdVector<Vector < Integer> > & actMap = it2->second;
      
      // iterate over all entries in the current map and renumber them
      for ( UInt i = 0; i < actMap.GetSize(); i++ ) {
	for ( UInt j = 0; j < actMap[i].GetSize(); j++ ) {
	  if ( actMap[i][j] > 0 ) {
	    actMap[i][j] = (*order)[actMap[i][j]-1];
	  }
	  else if(actMap[i][j] < 0 ) {
	    //due to constraints
	    actMap[i][j] = -(*order)[-actMap[i][j]-1];   
	  }
	}
      }
    }

    // Delete ordering array at the end
    delete [] (*order);
    (*order) = NULL;
  }
  
  
  //! ======================================================================
  //! EQUATION MAPPING
  //! ======================================================================
  void EqnMap::GetEqns( StdVector<Integer> &eqns,
		        const ResultInfo& result, 
		        const EntityIterator& it ) const{
    ENTER_IFCN( "EqnMap::GetEqns" );
    
    UInt numDofs = result.dofNames.GetSize();          

    // first of all, delete eqns-array
    eqns.Clear();
    
    // temporary
    

    // ============
    //  NODAL PART 
    // ============

    if ( result.definedOn == ResultInfo::NODE ||
	 result.definedOn == ResultInfo::PFEM ) {
      
      // get related nodal equaiton map
      Matrix<Integer> const & map = (nodeEqns_.find( result ) )->second;

      // Distinguish the type the of the list
      if ( it.GetType() == EntityList::ELEM_LIST
	   || it.GetType() == EntityList::SURF_ELEM_LIST ) {
	StdVector<UInt> const & nodes = it.GetElem()->connect;
	
	eqns.Resize( nodes.GetSize() * numDofs );
	eqns.Init();
	
	for (UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
	  Integer localNode = mesh2PdeNode_[ nodes[iNode]-1 ];
	  
	  for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
	    
	    if (localNode < 1 ) {
	      eqns[iNode*numDofs + iDof] = 0;
	    } else {
	      eqns[iNode*numDofs + iDof] = map[localNode-1][iDof];
	    }
	    
	  }
	}
      } else if ( it.GetType() == EntityList::NODE_LIST ) {
	UInt node = it.GetNode();
	UInt numDofs = result.dofNames.GetSize();
	eqns.Resize( numDofs );
	eqns.Init();
	
	Integer localNode = mesh2PdeNode_[ node-1 ];
	
	for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
	  
	  if (localNode < 1 ) {
	    eqns[iDof] = 0;
	  } else {
	    eqns[iDof] = map[localNode-1][iDof];
	  }
	    
	}
      } else {
	EXCEPTION( "This type of entity list is not defined for " 
                   << "equation mapping!" );
      }
    }

    // ===============
    //   EDGE PART 
    // ===============
    
    if ( result.definedOn == ResultInfo::EDGE ||
	 result.definedOn == ResultInfo::PFEM ) {
      
      // Check if entity is of type NODE:
      // In this case, we have to leave, as currently only nodal
      // boundary conditions can be assigned.
      if( it.GetType() == EntityList::NODE_LIST ) {
	return;
      }
      
      // get edge map of current result
      StdVector<Vector<Integer> > const & map = 
	(edgeEqns_.find( result ) )->second;
      
      // get edges
      StdVector<Integer> const & edges = it.GetElem()->edges;
      
      // iterate over all edges
      for( UInt iEdge = 0; iEdge < edges.GetSize(); iEdge++ ) {
	
	// get local edge number
	Integer locEdge = mesh2PdeEdge_[ std::abs(edges[iEdge]) -1 ];
	
	// iterate over all dofs of this edge
	for( UInt iDof = 0; iDof < map[locEdge-1].GetSize(); iDof++ ) {
	  eqns.Push_back( map[locEdge-1][iDof] );
	}
      }
    }

    // ===============
    //   FACE PART 
    // ===============

    if ( result.definedOn == ResultInfo::FACE ||
	 result.definedOn == ResultInfo::PFEM ) {
      
      // Check if entity is of type NODE:
      // In this case, we have to leave, as currently only nodal
      // boundary conditions can be assigned.
      if( it.GetType() == EntityList::NODE_LIST ) {
	return;
      }
      
      // get face map of current result
      StdVector<Vector<Integer> > const & map = 
	(faceEqns_.find( result ) )->second;
      
      // get faces
      StdVector<Integer> const & faces = it.GetElem()->faces;
      
      // iterate over all faces
      for( UInt iFace = 0; iFace < faces.GetSize(); iFace++ ) {
	
	// get local face number
	Integer locFace = mesh2PdeFace_[ faces[iFace] -1 ];
	
	// iterate over all dofs of this face
	for( UInt iDof = 0; iDof < map[locFace-1].GetSize(); iDof++ ) {
	  eqns.Push_back( map[locFace-1][iDof] );
	}
      }
    }
      


    // =============
    //  ELEM PART
    // =============
    if( result.definedOn == ResultInfo::ELEMENT  ||
        result.definedOn == ResultInfo::PFEM ) {
      StdVector< Vector<Integer> >const & elemMap = (elemEqns_.find( result ) )->second;
      UInt localElem = mesh2PdeElem_[(it.GetElem()->elemNum)-1];
      for (UInt iDof = 0; iDof < elemMap[localElem-1].GetSize(); iDof++ ) {
	if (localElem < 1 ) {
	  eqns.Push_back(0);
	} else {
	  eqns.Push_back( elemMap[localElem-1][iDof] );
          LOG_DBG3(eqnMap) << "Pushin back eqn " <<  elemMap[localElem-1][iDof]
                           << " for interior of element #" 
                           << it.GetElem()->elemNum << std::endl;
	}
	
      }
    }

    LOG_DBG3(eqnMap) << "Equations are: " << eqns.Serialize();
    LOG_DBG3(eqnMap) << "Number of equations: " << eqns.GetSize() << std::endl;
  }
    
  Integer EqnMap::GetEqn( const ResultInfo& result, 
		          const EntityIterator& it,
		          UInt dof ) const{
    StdVector<Integer> eqns;
    GetEqns( eqns, result, it );
    return eqns[dof-1];
  }

    


  Integer EqnMap::GetNodeEqn( const ResultInfo& result,
		              UInt nodeNr, UInt dof ) { 
    ENTER_FCN( "EqnMap::GetNodeEqn" );
    

    if ( result.definedOn == ResultInfo::NODE
	 || result.definedOn == ResultInfo::PFEM ) {
      
      
      Matrix<Integer> const & map = (nodeEqns_.find( result ) )->second;
      
      
      Integer localNode = mesh2PdeNode_[ nodeNr-1 ];
      
      if (localNode < 1 ) {
	return 0; 
      } else {
	return map[localNode-1][dof-1];
      }
      
    }    
    return 0;
  }

  Integer EqnMap::GetNodeEqn( UInt nodeNr, UInt dof ) {

    // First check, if more than one type of results are defined
    if ( nodeEqns_.size() != 1 ) {
      EXCEPTION( "GetNodeEqn() can only be used if exactly "
                 << "one nodal result is present!" );
    }

    Matrix<Integer> const & map = (nodeEqns_.begin() )->second;
    Integer localNode = mesh2PdeNode_[ nodeNr-1 ];
    
    if (localNode < 1 ) {
      return 0; 
    } else {
      return map[localNode-1][dof-1];
    }
  }

  
  void EqnMap::GetNodeEqn( const StdVector<UInt>& nodeNrs, 
		           StdVector<Integer>& eqnNrs ) {
    // First check, if more than one type of results are defined
    if ( nodeEqns_.size() != 1 ) {
      EXCEPTION( "GetNodeEqn() can only be used if exactly "
                 << "one nodal result is present!" );
    }
    Matrix<Integer> const & map = (nodeEqns_.begin() )->second;
    UInt dofsPerNode =  (nodeEqns_.begin() )->first.dofNames.GetSize();

    eqnNrs.Resize( nodeNrs.GetSize() * dofsPerNode );
    eqnNrs.Init();
    
    for( UInt iNode=0; iNode<nodeNrs.GetSize(); iNode++) {
      Integer localNode = mesh2PdeNode_[ nodeNrs[iNode]-1 ];
      
      for (UInt iDof = 0; iDof < dofsPerNode; iDof++ ) {
	
	if (localNode < 1 ) {
	  eqnNrs[iNode*dofsPerNode + iDof] = 0;
	} else {
	  eqnNrs[iNode*dofsPerNode + iDof] = map[localNode-1][iDof];
	}
      }
    }
      
  }
	 

    
  
  
  //! ======================================================================
  //! LOCAL/GLOBAL MAPPING OF MESH ENTITIES
  //! ======================================================================
  
  void EqnMap::Mesh2PdeNode(StdVector<UInt> & PdeNodes,
		            const StdVector<UInt> & MeshNodes) const {
    ENTER_FCN( "EqnMap::Mesh2PdeNode" );

    PdeNodes.Resize(MeshNodes.GetSize());
    PdeNodes.Init();
   
    for (UInt i=0; i<MeshNodes.GetSize(); i++) 
      PdeNodes[i] = mesh2PdeNode_[MeshNodes[i]-1];
  }
 
  UInt EqnMap::Mesh2PdeNode(const UInt meshNode) const {
    ENTER_FCN( "EqnMap::Mesh2PdeNode" );

    if ( mesh2PdeNode_[meshNode-1] < 0 ) {
      EXCEPTION(  "MeshNode Nr. " << meshNode 
                  << " has no local node number!" );
    }

    return abs(mesh2PdeNode_[meshNode-1]);
  }
    
  void EqnMap::Pde2MeshNode(StdVector<UInt> & meshNodes,
		            const StdVector<UInt> & pdeNodes) const {
    ENTER_FCN( "EqnMap::Pde2MeshNode" );
    meshNodes.Resize(pdeNodes.GetSize());
    meshNodes.Init();
   
    for (UInt i=0; i<pdeNodes.GetSize(); i++) 
      meshNodes[i] = pde2MeshNode_[pdeNodes[i]-1];
  }
  
  UInt EqnMap::Pde2MeshNode(const UInt pdeNode) const {
    ENTER_FCN( "EqnMap::Pde2MeshNode" );
    return pde2MeshNode_[pdeNode-1];
  }
  
  UInt EqnMap::Mesh2PdeElem(const UInt elemNumGlob) const {
    ENTER_FCN( "EqnMap::Mesh2PdeElem" );

    if ( mesh2PdeElem_[elemNumGlob-1] < 0 ) {
      EXCEPTION( "MeshElem Nr. " << elemNumGlob 
                 << " has no local elem number!" );
    }
    
    return abs(mesh2PdeElem_[elemNumGlob-1]);
  }
  
  UInt EqnMap::Pde2MeshElem(const UInt elemNumLoc) const {
    ENTER_FCN( "EqnMap::Pde2MeshElem" );
    return pde2MeshElem_[elemNumLoc-1];
  }
  
  //! ======================================================================
  //! MISCELLANEOUS
  //! ======================================================================
  
  void EqnMap::Print(std::ostream & out) const {
    ENTER_FCN( "EqnMap::Print" );
    out << "======================================\n"
	<< "  Equation numbering - Information\n"
	<< "======================================\n\n";

    // --- Local<->Global mapping (nodes) ---
    out << "1.) local <-> global mapping (nodes)\n"
	<< "------------------------------------\n"
	<< "#global nodes = " << ptGrid_->GetNumNodes() << "\n"
	<< "#local nodes = " << numLocNodes_ << "\n\n";
      
    out << std::setw(12) << "loc. nodeNr" << " | ";
    out << std::setw(12) << "glob. nodeNr" << "\n";
    out << std::setfill('-') << std::setw(27) << "-" << std::endl;
    out << std::setfill(' ');

    for( UInt iNode = 0; iNode < pde2MeshNode_.GetSize(); iNode++ ) {
      out << std::setw(12) << iNode+1 << "|";
      out << std::setw(12) << pde2MeshNode_[iNode];
      out << "\n";
    }
    out << "\n\n";

    // --- Local<->Global mapping (elements) ---
    out << "2.) local <-> global mapping (elements)\n"
	<< "------------------------------------\n"
	<< "#global nodes = " << ptGrid_->GetNumVolElems() << "\n"
	<< "#local nodes = " << numLocElems_ << "\n\n";
      
    out << std::setw(12) << "loc. elemNr" << " | ";
    out << std::setw(12) << "glob. elemNr" << "\n";
    out << std::setfill('-') << std::setw(27) << "-" << std::endl;
    out << std::setfill(' ');

    for( UInt iElem = 0; iElem < pde2MeshElem_.GetSize(); iElem++ ) {
      out << std::setw(12) << iElem+1 << "|";
      out << std::setw(12) << pde2MeshElem_[iElem];
      out << "\n";
    }
    out << "\n\n";

    // --- Local<->Global mapping (faces) ---
    out << "3.) local <-> global mapping (faces)\n"
	<< "------------------------------------\n\n";

    out << std::setw(12) << "glob. faceNr" << " | ";
    out << std::setw(12) << "loc. faceNr" << "\n";
    out << std::setfill('-') << std::setw(27) << "-" << std::endl;
    out << std::setfill(' ');

    for( UInt iFace = 0; iFace < mesh2PdeFace_.GetSize(); iFace++ ) {
      out << std::setw(12) << iFace+1;
      out << std::setw(12) << mesh2PdeFace_[iFace] << "|";
      out << "\n";
    }
    out << "\n\n";

    // --- Local<->Global mapping (edges) ---
    out << "4.) local <-> global mapping (edges)\n"
	<< "------------------------------------\n\n";

    out << std::setw(12) << "glob. edgeNr" << " | ";
    out << std::setw(12) << "loc. edgeNr" << "\n";
    out << std::setfill('-') << std::setw(27) << "-" << std::endl;
    out << std::setfill(' ');

    for( UInt iEdge = 0; iEdge < mesh2PdeEdge_.GetSize(); iEdge++ ) {
      out << std::setw(12) << iEdge+1;
      out << std::setw(12) << mesh2PdeEdge_[iEdge] << "|";
      out << "\n";
    }
    out << "\n\n";

    // --- Nodal equation mapping ---
    out << "5.) Nodal equation mapping\n"
	<< "------------------------------------\n";
    

    // Loop over all nodal mapped results
    ResultEntityMap::const_iterator listIt;
    for( listIt = nodeMappedList_.begin(); 
	 listIt != nodeMappedList_.end(); 
	 listIt++ ) {

      // get dofspernode and associated eqnMap

      Matrix<Integer> const & eqnMap = (nodeEqns_.find(listIt->first))->second;
      UInt dofsPerNode = listIt->first.dofNames.GetSize();
      
      // Print head for which result this is shown
      std::string type;
      Enum2String( listIt->first.resultType, type );
      out << "Result:" << type << "\n";
      out << "Dofs: " << dofsPerNode << "\n\n";

      out << std::setw(12) << "loc. nodeNr" << " | ";
      out << std::setw(12) << "glob. nodeNr" << " | ";
      out << std::setw(6) << "dof" << " | ";
      out << std::setw(12) << "eqnNr" << "\n";
      out << std::setfill('-') << std::setw(50) << "-" << std::endl;
      out << std::setfill(' ');

      for ( UInt iNode = 0; iNode < pde2MeshNode_.GetSize(); iNode++ ) {
	
	// Print out first dof
	out << std::setw(12) << iNode+1  << " | ";
	out << std::setw(12) << pde2MeshNode_[iNode] << " | ";
	out << std::setw(6) << 1 << " | ";
	
	out << std::setw(12) << eqnMap[iNode][0] << std::endl;
	
	// Print out further dofs if needed
	for ( UInt iDof = 1; iDof < dofsPerNode; iDof++ ) {
	  out << std::setw(12) << " " << " | ";
	  out << std::setw(12) << " " << " | ";
	  out << std::setw(6) << iDof+1 << " | ";
	  out << std::setw(12) << eqnMap[iNode][iDof] << std::endl;
	}
      }
    }

    out << "\n\n";

    // --- Edge equation mapping ---
    out << "6.) Edge equation mapping\n"
	<< "------------------------------------\n";
    

    // Loop over all edge mapped results
    ResultEntityMap::const_iterator edgeListIt;
    for( edgeListIt = edgeMappedList_.begin(); 
	 edgeListIt != edgeMappedList_.end(); 
	 edgeListIt++ ) {

      StdVector<Vector<Integer> >  const & eqnMap = 
	(edgeEqns_.find(edgeListIt->first))->second;
      
      // Print head for which result this is shown
      std::string type;
      Enum2String( edgeListIt->first.resultType, type );
      out << "Result:" << type << "\n\n";

      out << std::setw(12) << "glob. edgeNr " << " | ";
      out << std::setw(12) << "loc. edgeNr" << " | ";
      out << std::setw(6) << "dof" << " | ";
      out << std::setw(12) << "eqnNr" << "\n";
      out << std::setfill('-') << std::setw(50) << "-" << std::endl;
      out << std::setfill(' ');

      UInt numEdges = ptGrid_->GetNumEdges();
      for ( UInt iEdge = 0; iEdge < numEdges; iEdge++ ) {

	// local edge numger
	Integer locEdge = mesh2PdeEdge_[iEdge];
	
	// check if edge has any equations at all
	if ( locEdge <= 0 ) {
	  break;
	}
	
	// Print out first dof
	out << std::setw(12) << iEdge+1  << " | ";
	out << std::setw(12) << mesh2PdeEdge_[iEdge] << " | ";
	out << std::setw(6) << 1 << " | ";
	

	out << std::setw(12);
	if( locEdge > 0 && eqnMap[locEdge-1].GetSize() > 0) {
	  out << eqnMap[locEdge-1][0] << std::endl;

	  // Print out further dofs if needed
	  for ( UInt iDof = 1; iDof < eqnMap[locEdge-1].GetSize(); iDof++ ) {
	    out << std::setw(12) << " " << " | ";
	    out << std::setw(12) << " " << " | ";
	    out << std::setw(6) << iDof+1 << " | ";
	    out << std::setw(12) << eqnMap[locEdge-1][iDof] << std::endl;
	  }
	} else {
	  out << "-" << std::endl;
	}
      }
    }
    out << "\n\n";


    
    // --- Face equation mapping ---
    out << "7.) Face equation mapping\n"
	<< "------------------------------------\n";
    

    // Loop over all face mapped results
    ResultEntityMap::const_iterator faceListIt;
    for( faceListIt = faceMappedList_.begin(); 
	 faceListIt != faceMappedList_.end(); 
	 faceListIt++ ) {

      StdVector<Vector<Integer> >  const & eqnMap = 
	(faceEqns_.find(faceListIt->first))->second;
      
      // Print head for which result this is shown
      std::string type;
      Enum2String( faceListIt->first.resultType, type );
      out << "Result:" << type << "\n\n";

      out << std::setw(12) << "glob. faceNr " << " | ";
      out << std::setw(12) << "loc. faceNr" << " | ";
      out << std::setw(6) << "dof" << " | ";
      out << std::setw(12) << "eqnNr" << "\n";
      out << std::setfill('-') << std::setw(50) << "-" << std::endl;
      out << std::setfill(' ');

      UInt numFaces = ptGrid_->GetNumFaces();
      for ( UInt iFace = 0; iFace < numFaces; iFace++ ) {

	// local face numer
	Integer locFace = mesh2PdeFace_[iFace];
	
	// check if face has any equations at all
	if ( locFace <= 0 ) {
	  break;
	}
	
	// Print out first dof
	out << std::setw(12) << iFace+1  << " | ";
	out << std::setw(12) << mesh2PdeFace_[iFace] << " | ";
	out << std::setw(6) << 1 << " | ";
	

	out << std::setw(12);
	if( locFace > 0 && eqnMap[locFace-1].GetSize() > 0) {
	  out << eqnMap[locFace-1][0] << std::endl;

	  // Print out further dofs if needed
	  for ( UInt iDof = 1; iDof < eqnMap[iFace].GetSize(); iDof++ ) {
	    out << std::setw(12) << " " << " | ";
	    out << std::setw(12) << " " << " | ";
	    out << std::setw(6) << iDof+1 << " | ";
	    out << std::setw(12) << eqnMap[locFace-1][iDof] << std::endl;
	  }
	} else {
	  out << "-" << std::endl;
	}
      }
    }

    // --- Interior equation mapping ---
    out << "8.) Interior equation mapping\n"
	<< "------------------------------------\n";
    

    // Loop over all element interior mapped results
    ResultEntityMap::const_iterator elemIntListIt;
    for( elemIntListIt = elemIntMappedList_.begin(); 
	 elemIntListIt != elemIntMappedList_.end(); 
	 elemIntListIt++ ) {

      StdVector<Vector<Integer> >  const & eqnMap = 
	(elemEqns_.find(elemIntListIt->first))->second;
      
      // Print head for which result this is shown
      std::string type;
      Enum2String( elemIntListIt->first.resultType, type );
      out << "Result: " << type << "\n\n";

      out << std::setw(12) << "loc. elemNr " << " | ";
      out << std::setw(12) << "glob. elemNr" << " | ";
      out << std::setw(6) << "dof" << " | ";
      out << std::setw(12) << "eqnNr" << "\n";
      out << std::setfill('-') << std::setw(50) << "-" << std::endl;
      out << std::setfill(' ');

      for ( UInt iElem = 0; iElem < numLocElems_; iElem++ ) {

	
	// Print out first dof
	out << std::setw(12) << iElem+1  << " | ";
	out << std::setw(12) << pde2MeshElem_[iElem] << " | ";
	out << std::setw(6) << 1 << " | ";
	

	out << std::setw(12);
	if(  eqnMap[iElem].GetSize() > 0) {
	  out << eqnMap[iElem][0] << std::endl;

	  // Print out further dofs if needed
	  for ( UInt iDof = 1; iDof < eqnMap[iElem].GetSize(); iDof++ ) {
	    out << std::setw(12) << " " << " | ";
	    out << std::setw(12) << " " << " | ";
	    out << std::setw(6) << iDof+1 << " | ";
	    out << std::setw(12) << eqnMap[iElem][iDof] << std::endl;
	  }
	} else {
	  out << "0" << std::endl;
	}
      }
    }

    out << "\n\n";

  
    
  }


  //! ======================================================================
  //! PRIVATE HELPER METHODS
  //! ======================================================================

  
  void EqnMap::CalcNodeElemMapping() {
    ENTER_FCN( "EqnMap::CalcNodeElemMapping" );
    
    mesh2PdeNode_.Resize( ptGrid_->GetNumNodes() );
    mesh2PdeNode_.Init( -1 );
    pde2MeshNode_.Clear( );

    mesh2PdeElem_.Resize( ptGrid_->GetNumElems() );
    mesh2PdeElem_.Init( -1 );
    // Note: here we could iterate over all element lists and 
    //       count the number of entries, so we would now from beginning,
    //       how many elements we have in this pde
    pde2MeshElem_.Clear( );
 

    UInt nodeCounter = 0;
    UInt elemCounter = 1;
    StdVector<Elem*> subdom;
 
    // Note: The check, if the elemList of this region corresponds to all 
    //       elements of the pde can not yet be performed
    // First, check, if Pde is defined on all regions of domain
    //  StdVector<RegionIdType> allRegions;
    bool allFound = false;
    //     Integer index = -1;

    //     ptGrid_->GetVolRegionIds( allRegions );
    //     for ( UInt i = 0; i < allRegions.GetSize(); i++ ) {
    //       index = subdoms_.Find(allRegions[i]);
    //       if ( index == -1 ) {
    //         allFound = false;
    //         break;
    //       }
    //     }
    
    // Case 1: This pde is defined on all volume regions of the grid
    //         --> Do a simple 1-to-1 mapping of global node number
    //             and local one
    if ( allFound ) {
      //    pde2MeshNode_.Resize(mesh2PdeNode_.GetSize());
      //       for (UInt i = 0; i<mesh2PdeNode_.GetSize(); i++) {
      //         mesh2PdeNode_[i] = i+1;
      //         pde2MeshNode_[i] = i+1;
      //     }
      
      //       pde2MeshElem_.Resize(ptGrid_->GetNumElems());
      //       for (UInt i = 0; i<mesh2PdeElem_.GetSize(); i++) {
      //         pde2MeshElem_[i] = i+1;
      //         mesh2PdeElem_[i] = i+1;
      //       }
      
      
      //       nodeCounter=mesh2PdeNode_.GetSize();
    } else {
      // Case 2: This pde is defined on a subset of all regions
      //         --> Perform normal node renumbering
      
      //  // iterate over all element lists
      for ( UInt iList = 0; iList < locElems_.GetSize(); iList++ ) {

	// Get iterator of current element list
	EntityIterator it = locElems_[iList]->GetIterator();
	
	// iterate over all elements in element list
	for ( it.Begin(); !it.IsEnd(); it++ ) {

	  // Store current element
	  const Elem* actEl = it.GetElem();
	  
	  // *** Mapping of Elements ***
	  mesh2PdeElem_[actEl->elemNum - 1 ] = elemCounter;
	  pde2MeshElem_.Push_back( actEl->elemNum );
	  elemCounter++;
	  
	  
	  // *** Mapping of Nodes ***
	  // iterate over all nodes in elem
	  for ( UInt iNode=0; iNode < actEl->connect.GetSize(); iNode++)

	    // Check if node was already assigned
	    if (mesh2PdeNode_[actEl->connect[iNode]-1] == -1) {
	      mesh2PdeNode_[actEl->connect[iNode]-1] = ++nodeCounter;
	      pde2MeshNode_.Push_back(actEl->connect[iNode]);
	    }
	}
      }
    }
    
    // remember number of local nodes and elements
    numLocNodes_ = pde2MeshNode_.GetSize();
    numLocElems_ = pde2MeshElem_.GetSize();
  }
  
  
  void EqnMap::CalcEdgeMapping()  {
    ENTER_FCN( "EqnMap::CalcEdgeMapping" );

    mesh2PdeEdge_.Resize( ptGrid_->GetNumEdges() );
    mesh2PdeEdge_.Init( -1 );

    UInt edgeCounter = 0;
    
    // iterate over all element lists
    for ( UInt iList = 0; iList < locElems_.GetSize(); iList++ ) {

      // Get iterator of current element list
      EntityIterator it = locElems_[iList]->GetIterator();
      
      // iterate over all elements in element list
      for ( it.Begin(); !it.IsEnd(); it++ ) {
	
	// Store current element
	const Elem* actEl = it.GetElem();
	
	// iterate over all edges
	for ( UInt iEdge=0; iEdge < actEl->edges.GetSize(); iEdge++) {
	  
	  // Check if edge was already assigned
	  if( mesh2PdeEdge_[std::abs(actEl->edges[iEdge])-1] 
	      == -1 ) {
	    mesh2PdeEdge_[std::abs(actEl->edges[iEdge])-1] 
	      = ++edgeCounter;
	  }
	}
      }
    }
    
    // store number of local edges
    numLocEdges_ = edgeCounter;
    
  }
  

  void EqnMap::CalcFaceMapping()  {
    ENTER_FCN( "EqnMap::CalcFaceMapping" );

    LOG_TRACE(eqnMap) << "Starting local<->global face mapping\n";
    
    mesh2PdeFace_.Resize( ptGrid_->GetNumFaces() );
    mesh2PdeFace_.Init( -1 );
    
    UInt faceCounter = 0;
    
    // iterate over all element lists
    for ( UInt iList = 0; iList < locElems_.GetSize(); iList++ ) {
      
      // Get iterator of current element list
      EntityIterator it = locElems_[iList]->GetIterator();
      
      // iterate over all elements in element list
      for ( it.Begin(); !it.IsEnd(); it++ ) {
	
	// Store current element
	const Elem* actEl = it.GetElem();
	
	// iterate over element faces
	for ( UInt iFace=0; iFace < actEl->faces.GetSize(); iFace++) {
          
	  // Check if face was already assigned
	  if( mesh2PdeFace_[actEl->faces[iFace]-1] 
	      == -1 ) {
	    mesh2PdeFace_[actEl->faces[iFace]-1] 
              = ++faceCounter;
	  }
	}
      }
    }
    
    // store number of local faces
    numLocFaces_ = faceCounter;

    LOG_DBG(eqnMap) << "There are " << numLocFaces_ << " local faces";
    LOG_TRACE(eqnMap) << "Finished local<->global face mapping\n";
    
  }
  
  void EqnMap::CalcNodalEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcNodalEquations" );


    // Big outer loop over all nodal mapped element lists
    ResultEntityMap::iterator listIt;

    for( listIt = nodeMappedList_.begin(); 
	 listIt != nodeMappedList_.end(); 
	 listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;


      // Get grip of homogeneous and in-homogeneous boundary conditions
      // for this tpye of result
      ResultHdBcMap::iterator hdBcIt = hdBcs_.find( actRes );
      ResultIdBcMap::iterator idBcIt = idBcs_.find( actRes );
      ResultConstraintMap::iterator csIt = constraints_.find( actRes );

      
      std::string type;
      Enum2String( actRes.resultType, type );
      Matrix<Integer> & actMap = nodeEqns_[actRes];

      // Get number of dofs
      UInt dofsPerNode = actRes.dofNames.GetSize();
    
      // Idea of the algorithm:
      //
      // -- PHASE 1 --
      // Step 1:  Initialize actMap with -1
      // Step 2:  For each entry in homoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 2b: For each entry in inhomoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 3:  For each entry in constraintSlaveNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 4:  Loop over all nodes of given entity lists
      //          and assign each non-zero entry an equation number
      // Step 5:  Loop over the whole map and set all entries with -1 to 0

      // -- PHASE 2 --
      // Step 6:  Afterwards loop again over all nodes in constraintSlaveNodes_
      //          and set the corresponding entry in pdeNode2EQN_ to the
      //          negative of the value of constraintMasterNode
      // Step 6b: Loop over all entries in inhomoDirichletNodes_ and assign that
      //          dof an equation number after the hightest equation number of
      //          the free dofs
      //
      // Note:    Steps 2b and 5b are only performed, if sortEQNs = true
      
      if( phase == 1 ) {
	
	// ------
	// STEP 1
	// ------
	//UInt multipleBCs = 0;
	
	actMap.Resize( numLocNodes_, dofsPerNode );
	actMap.Init( -1 );
	
	
	// ------
	// STEP 2
	// ------
	
	
	Matrix<UInt> countNodes;
	countNodes.Resize( numLocNodes_, dofsPerNode );
	countNodes.Init( 0 );
	
	if( hdBcIt != hdBcs_.end() ) {
          HdBcList const & actHdBcList = hdBcIt->second;
	  for ( UInt i = 0; i < actHdBcList.GetSize(); i++ ) {
	    StdVector<UInt> nodes;
	    GetNodesOfEntities( nodes, actHdBcList[i]->entities );
	    UInt actDof = actHdBcList[i]->dof;
	    
	    for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
	      // Check if homDirichletNode belongs to one of my subdomains
	      if ( mesh2PdeNode_[nodes[iNode]-1]-1 < 0 ) {
		(*warning) << "EqnMap::CalcNodalEquations: Homogen. Dirichlet node "
		           << "nr. " << nodes[iNode]
		           << " is not contained in any of the regions for this PDE";
		Warning( __FILE__, __LINE__ );
	      }
	      else if ( countNodes[mesh2PdeNode_[nodes[iNode]-1]-1]
		        [actDof-1] != 0 ) {
// 		(*warning) << "EqnMap::CalcNodalEquations: HomDirichletNode # "
// 		           << nodes[i]
// 		           << "\nappeared already at least once in the list of "
// 		           << "boundary nodes for this PDE!\n Please check, if this "
// 		           << "node is defined in more than one level of boundary "
// 		           << "nodes!";
// 		Warning( __FILE__, __LINE__ );
	      }
	      else {
                actMap[mesh2PdeNode_[nodes[iNode]-1]-1][actDof-1] = 0;
                countNodes[mesh2PdeNode_[nodes[iNode]-1]-1][actDof-1]++;
	      }
	    }
	  }
	}
	
	
	// -------
	// STEP 2b
	// -------

	countNodes.Init(0);

	// Check if any inhom. boundary condition is defined for the current
	// result
	if( idBcIt != idBcs_.end() ) {
	  IdBcList const & actIdBcList = idBcIt->second;

	  for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
	    StdVector<UInt> nodes;
	    GetNodesOfEntities( nodes, actIdBcList[i]->entities );
	    UInt actDof = actIdBcList[i]->dof;

	    for( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
	      
	      if ( mesh2PdeNode_[ nodes[iNode] - 1 ] - 1 < 0 ) {
		(*warning) << "EqnMap::CalcNodalEquations: Inhom. Dirichlet "
		           << "node #" << nodes[iNode]
		           << " is not contained in any of the regions for "
		           << "this Pde";
		Warning( __FILE__, __LINE__ );
	      }
	      else if ( countNodes[mesh2PdeNode_[nodes[iNode]-1]-1]
		        [actDof-1] != 0 ) {
	// 	(*warning) << "EqnMap::CalcNodalEquations: Inhom. Dirichlet "
// 		           << "node #" << nodes[iNode]
// 		           << "\nappeared already at least once in the list of "
// 		           << "boundary nodes for this Pde!\n Please check, if "
// 		           << "this node is defined in more than one level of "
// 		           << "boundary nodes!";
// 		Warning( __FILE__, __LINE__ );
	      }
	      else {

                // Only set equation number to zero, if we sort equations

                if( sortEqns_ ) {
                  actMap[mesh2PdeNode_[nodes[iNode]-1]-1] [actDof-1] = 0;
                  countNodes[mesh2PdeNode_[nodes[iNode]-1]-1][actDof-1]++;
                }
		// In any case we have to increment the number of idBC-conditions
		numIdBcs_++;
	      }
	    }
	  }
	}
	
	
	// ------
	// STEP 3
	// ------
        
	// Check if any inhom. constraint is defined for the current
	// result
	if( csIt != constraints_.end() ) {
	  ConstraintList const & actCsList = csIt->second;
	  
	  for ( UInt i = 0; i < actCsList.GetSize(); i++ ) {
	    StdVector<UInt> slaveNodes;
	    GetNodesOfEntities( slaveNodes, actCsList[i]->slaveEntities );
	    UInt slaveDof = actCsList[i]->slaveDof;
	    
	    //for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
            for ( UInt iNode = 0; iNode < slaveNodes.GetSize(); iNode++ ) {
	      actMap[mesh2PdeNode_[slaveNodes[iNode]-1]-1] [slaveDof-1] = 0;
	    }
	  }
	}
	
	 
        // ------
	// STEP 4
	// ------
	// Initialize countNodes to zero. It will be used to count if
	// a node got already an equation number
	countNodes.Init(0);
	
	StdVector<UInt> nodes;
	// Iterate over all element list belonging to this result
	// and assign each of them a node number
	for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {
	  
	  // Get nodes of current entityList
	  GetNodesOfEntities( nodes, actLists[iList] );
	  
	  Integer locNode = 0;
	  for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
	    locNode = mesh2PdeNode_[nodes[iNode]-1];
	    for ( UInt iDof = 0; iDof < dofsPerNode; iDof++ ) {
	      if ( actMap[locNode-1][iDof] != 0 &&
		   countNodes[locNode-1][iDof] == 0) {
		numEqns_++;
                LOG_DBG3(eqnMap) << "Adding equation " << numEqns_ 
                                 << " to local node " << locNode << std::endl;
		actMap[locNode-1][iDof] = numEqns_;
		countNodes[locNode-1][iDof] = 1;
	      }
	    }
	  }
	}

        // ------
	// STEP 5
	// ------
        // Re-iterate over the whole equation map and set all entries 
        // with eqn-number of -1 to 0
        for( UInt i = 0; i < actMap.GetSizeRow(); i++ ) {
          for( UInt j = 0; j < actMap.GetSizeCol(); j++ ) {
            if( actMap[i][j] == -1 )
              actMap[i][j] = 0;
          }
        }
	
        LOG_DBG2(eqnMap) << "Final equation map looks like: \n" 
                         << actMap << std::endl;
        
      } else if( phase == 2 ) {
	
	// ------
	// STEP 6
	// ------
	// Check if any inhom. constraint is defined for the current
	// result
	if( csIt != constraints_.end() ) {
	  ConstraintList const & actCsList = csIt->second;
	  
	  for ( UInt i = 0; i < actCsList.GetSize(); i++ ) {
	    StdVector<UInt> slaveNodes;
	    GetNodesOfEntities( slaveNodes, actCsList[i]->slaveEntities );
	    UInt slaveDof = actCsList[i]->slaveDof;
	    UInt masterDof = actCsList[i]->masterDof;
	    UInt masterNode = slaveNodes[0];

            // assign the master node/dof a equation number
            actMap[mesh2PdeNode_[slaveNodes[0]-1]-1]
              [slaveDof-1] = ++numEqns_;

	    for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
	      actMap[mesh2PdeNode_[slaveNodes[iNode]-1]-1] [slaveDof-1] =
		-actMap[mesh2PdeNode_[masterNode-1]-1] [masterDof-1];
              numCs_++;
	    }
	  }
	}
	

	// -------
	// STEP 6b
	// -------
	if( idBcIt != idBcs_.end() ) {
	  IdBcList const & actIdBcList = idBcIt->second;
	  for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
	    StdVector<UInt> nodes;
	    GetNodesOfEntities( nodes, actIdBcList[i]->entities );
	    UInt actDof = actIdBcList[i]->dof;
	    
	    for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
              // only assign an equation number, if the map contains
              // a 0. Otherwise, we have already labeled this node
              Integer locNode = mesh2PdeNode_[nodes[iNode]-1];
              if( locNode > 0 ) { 
                if(  actMap[locNode-1] [actDof-1] 
                     == 0 ) {
                  numEqns_++;
                  actMap[locNode-1] [actDof-1] = numEqns_;
                }
              }
	    }
	  }
	}
	
      } else {
	EXCEPTION( "Phase '" << phase << "' does not exist!" );
      }
    }
  }
   
  void EqnMap::CalcElemInteriorEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcElemInteriorEquations" );
    
    // Big outer loop over all element interior mapped lists
    ResultEntityMap::iterator listIt;
    
    for( listIt = elemIntMappedList_.begin(); 
         listIt != elemIntMappedList_.end(); 
	 listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      StdVector<Vector<Integer> > & actMap = elemEqns_[actRes];
      
      UInt dofsPerElem = actRes.dofNames.GetSize();
      actMap.Resize( numLocElems_  );

      // Fetch also nodal map
      Matrix<Integer> & actNodeMap = nodeEqns_[actRes];

      Integer locElem = 0;

      // Iterate over all entitylists
      for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

	EntityIterator it = actLists[iList]->GetIterator();
	
	// Iterate over all elements within this list
	for( it.Begin(); !it.IsEnd(); it++ ) {

	  // Get grip of element
	  const Elem & actEl = *(it.GetElem());

	  // Get number of unknowns for element
          StdVector<Vector<UInt> > numFcns; 
          numFcns.Resize( dofsPerElem );

          // iterate over all dofs of this result
          for( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
            actEl.ptElem->GetNumFncs( numFcns[iDof], actRes.fctType, 
                                      AnsatzFct::INTERIOR, iDof );
          }

          // determine local element number
          locElem = mesh2PdeElem_[actEl.elemNum - 1];
          
          // Check if this elem was already mapped
          if ( actMap[locElem-1].GetSize() == 0 && locElem > 0 ) {
            UInt sum = 0;
            UInt max = 0;
            for( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
              sum +=  numFcns[iDof][0];
              if( numFcns[iDof][0] > max )
                max = numFcns[iDof][0];
            }
            //LOG_DBG3(eqnMap) << "1-dof: " << numFcns[0].Serialize() << std::endl;
            //LOG_DBG3(eqnMap) << "2-dof: " << numFcns[1].Serialize() << std::endl;
            //LOG_DBG3(eqnMap) << "3-dof: " << numFcns[2].Serialize() << std::endl;
            //LOG_DBG3(eqnMap) << "maximum for this elem"  << max << std::endl;
            //LOG_DBG3(eqnMap) << "Elem got " << sum << " equation numbers!\n";
             

            // iterate over all element functions
            UInt pos = 0;
            for( UInt iFcn = 0; iFcn < max; iFcn++ ) {
              if( actMap[locElem-1].GetSize() == 0 ) {
                actMap[locElem-1].Resize( dofsPerElem * max );
                actMap[locElem-1].Init(0);
              }
              
              // iterate over all element dofs
              for( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
                

                // Check if the related nodes have a 0 equations number
                // (= are (in-)hom. Dirichlet nodes)
                bool allZero = true;
                StdVector<UInt> const & elemNodes = actEl.connect;
                for( UInt iNode = 0; iNode < elemNodes.GetSize(); iNode++ ) {
                  if( actNodeMap[mesh2PdeNode_[elemNodes[iNode]-1]-1][iDof] != 0) {
                    allZero = false;
                    break;
                  }
                }
                if(  !allZero && iFcn < numFcns[iDof][0] ) {
                  actMap[locElem-1][pos++] = ++numEqns_;
                } else {
                  actMap[locElem-1][pos++] = 0;
                } // check for number of functions
              } // loop over dofs
            } // loop over functions
          } // check if element was already mapped
        } // loop over elements
      } // loop over entitylists
    } // loop over results
  }
  

  void EqnMap::CalcElemConstEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcElemEquations" );

    // Big outer loop over all nodal mapped element lists
    ResultEntityMap::iterator listIt;
    
    for( listIt = elemConstMappedList_.begin(); 
	 listIt != elemConstMappedList_.end(); 
	 listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      
      
      // Get grip of homogeneous and in-homogeneous boundary conditions
      // for this tpye of result
      ResultHdBcMap::iterator hdBcIt = hdBcs_.find( actRes );
      ResultIdBcMap::iterator idBcIt = idBcs_.find( actRes );
      ResultConstraintMap::iterator csIt = constraints_.find( actRes );
      
      
      std::string type;
      Enum2String( actRes.resultType, type );
      StdVector<Vector<Integer> > & actMap = elemEqns_[actRes];

      // Get number of dofs
      UInt dofsPerElem = actRes.dofNames.GetSize();
    
      // Idea of the algorithm:
      //
      // -- PHASE 1 --
      // Step 1:  Initialize actMap with 1
      // Step 2:  For each entry in homoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 2b: For each entry in inhomoDirichletNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 3:  For each entry in constraintSlaveNodes_ set the corresponding
      //          entry in pdeNode2eqn_ to 0
      // Step 4:  Loop over all entries in pde2Meshnode
      //          and assign each non-zero entry an equation number

      // -- PHASE 2 --
      // Step 5:  Afterwards loop again over all nodes in constraintSlaveNodes_
      //          and set the corresponding entry in pdeNode2EQN_ to the
      //          negative of the value of constraintMasterNode
      // Step 5b: Loop over all entries in inhomoDirichletNodes_ and assign that
      //          dof an equation number after the hightest equation number of
      //          the free dofs
      //
      // Note:    Steps 2b and 5b are only performed, if sortEQNs = true
      
      if( phase == 1 ) {
	
	// ------
	// STEP 1
	// ------
	//UInt multipleBCs = 0;
	Matrix<UInt> countElems;
	countElems.Resize( numLocElems_, dofsPerElem );
	countElems.Init( 0 );
	
	actMap.Resize( numLocElems_);
	for( UInt i = 0; i < actMap.GetSize(); i++) {
	  actMap[i].Resize( dofsPerElem );
	  actMap[i].Init( 1 );
	}

	// -------
	// STEP 2b
	// -------

	countElems.Init(0);

	// Check if any inhom. boundary condition is defined for the current
	// result
	if( idBcIt != idBcs_.end() ) {
	  IdBcList const & actIdBcList = idBcIt->second;

	  for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
	    StdVector<UInt> nodes;
	    EntityIterator elemIt = actIdBcList[i]->entities->GetIterator();
	    
	    UInt actDof = actIdBcList[i]->dof;
	    
	    for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
	      UInt actElem = elemIt.GetElem()->elemNum;
	      if ( mesh2PdeElem_[ actElem - 1 ] - 1 < 0 ) {
		(*warning) << "EqnMap::CalcElemEquations: Inhom. Dirichlet "
		           << "elem #" << actElem
		           << " is not contained in any of the regions for "
		           << "this Pde";
		Warning( __FILE__, __LINE__ );
	      }
	      else if ( countElems[mesh2PdeElem_[actElem-1]-1]
		        [actDof-1] != 0 ) {
		(*warning) << "EqnMap::CalcElemEquations: Inhom. Dirichlet "
		           << "elem #" << actElem
		           << "\nappeared already at least once in the list of "
		           << "boundary nodes for this Pde!\n Please check, if "
		           << "this node is defined in more than one level of "
		           << "boundary nodes!";
		Warning( __FILE__, __LINE__ );
	      }
	      else {
		if ( sortEqns_ == true ) {
		  actMap[mesh2PdeElem_[actElem-1]-1] [actDof-1] = 0;
		  countElems[mesh2PdeElem_[actElem-1]-1][actDof-1]++;
		}
		// In any case we have to increment the number of idBC-conditions
		numIdBcs_++;
	      }
	    }
	  }
	}
	
	// ------
	// STEP 4
	// ------
	// Initialize countNodes to zero. It will be used to count if
	// a node got already an equation number
	countElems.Init(0);
	
	StdVector<Elem*> elems;
	// Iterate over all element list belonging to this result
	// and assign each of them a node number
	for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

	  EntityIterator it = actLists[iList]->GetIterator();
	  
	  Integer locElem = 0;
	  for ( it.Begin(); !it.IsEnd(); it++ ) {
	    locElem = mesh2PdeElem_[(it.GetElem()->elemNum)-1];

	    // Check if this element was already mapped
	    if ( actMap[locElem-1].GetSize() == 0 ) {
	      actMap[locElem-1].Resize( dofsPerElem );
	      actMap[locElem-1].Init(0);
	    }
	    
	    for ( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
	      if ( actMap[locElem-1][iDof] != 0 &&
		   countElems[locElem-1][iDof] == 0) {
		numEqns_++;
		actMap[locElem-1][iDof] = numEqns_;
		countElems[locElem-1][iDof] = 1;
	      }
	    }
	  }
	}
	
      } else if( phase == 2 ) {
	
	
	
      } else {
	EXCEPTION( "Phase '" << phase << "' does not exist!" );
      }
    }

  }
  
  void EqnMap::CalcEdgeEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcEdgeEquations" );

    // Big outer loop over all edge mapped element lists
    ResultEntityMap::iterator listIt;

    for( listIt = edgeMappedList_.begin(); 
	 listIt != edgeMappedList_.end(); 
	 listIt++ ) {

      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      StdVector<Vector<Integer> > & actMap = edgeEqns_[actRes];
      
      UInt dofsPerEdge = actRes.dofNames.GetSize();
      actMap.Resize( numLocEdges_  );

      // Fetch also nodal map
      Matrix<Integer> & actNodeMap = nodeEqns_[actRes];

      Integer locEdge = 0;

      // Iterate over all entitylists
      for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {
	EntityIterator it = actLists[iList]->GetIterator();
	
	// Iterate over all elements within this list
	for( it.Begin(); !it.IsEnd(); it++ ) {

	  // Get grip of element
	  const Elem & actEl = *(it.GetElem());
                      
          // Get number of unknowns for each edge
          StdVector<Vector<UInt> > numFcns; 
          numFcns.Resize( dofsPerEdge );

          // iterate over all dofs of this result
          for( UInt iDof = 0; iDof < dofsPerEdge; iDof++ ) {
            actEl.ptElem->GetNumFncs( numFcns[iDof], actRes.fctType, 
                                      AnsatzFct::EDGE, iDof );
            //LOG_DBG3(eqnMap) << "numFncs, dof " << iDof+1 << ": " 
            //          << numFcns[iDof].Serialize() << std::endl;
            // assert that we have as many entries as we have edges
            assert( numFcns[iDof].GetSize() == actEl.edges.GetSize() );
          }
          
          // Iterate over all edges of this element
          for( UInt iEdge = 0; iEdge < actEl.edges.GetSize(); iEdge++ ) {
            
            // determine local edge
            locEdge = mesh2PdeEdge_[std::abs(actEl.edges[iEdge]) - 1];
            Edge const & edge = ptGrid_->GetEdge( std::abs(actEl.edges[iEdge] ));
            
            // Check if this edge was already mapped
            if ( actMap[locEdge-1].GetSize() == 0 && locEdge > 0 ) {
              
              // sum up unknowns of this dof
              UInt sum = 0;
              UInt max = 0;
              for( UInt iDof = 0; iDof < dofsPerEdge; iDof++ ) {
                sum +=  numFcns[iDof][iEdge];
                if( numFcns[iDof][iEdge] > max )
                  max = numFcns[iDof][iEdge];
              }
              //LOG_DBG3(eqnMap) << "1-dof: " << numFcns[0].Serialize();
              //LOG_DBG3(eqnMap) << "2-dof: " << numFcns[1].Serialize();
              //LOG_DBG3(eqnMap) << "3-dof: " << numFcns[2].Serialize();
              LOG_DBG3(eqnMap) << "maximum for edge " << iEdge 
                               << ": " << max;
              LOG_DBG3(eqnMap) << "Edge got " << sum 
                               << " equation numbers!" << std::endl;
              
              
              // iterate over all edge functions
              UInt pos = 0;
              for( UInt iFcn = 0; iFcn < max; iFcn++ ) {
                if( actMap[locEdge-1].GetSize() == 0 ) {
                  actMap[locEdge-1].Resize( dofsPerEdge * max );
                  actMap[locEdge-1].Init(0);
                }
                
                // iterate over all dofs
                for( UInt iDof = 0; iDof < dofsPerEdge; iDof++ ) {
                  
                  // Check if the related nodes have a 0 equations number
                  // (= are (in-)hom. Dirichlet nodes)
                  // of if the edge with this dof has a smaller order
                  // than the maximum for this edge -> assign 0 equation number
                  if(  (actNodeMap[mesh2PdeNode_[edge.nodes[0]-1]-1][iDof] != 0
                        || actNodeMap[mesh2PdeNode_[edge.nodes[1]-1]-1][iDof] != 0)
                       && (iFcn < numFcns[iDof][iEdge]) ) {
                    actMap[locEdge-1][pos++] = ++numEqns_;
                  } else {
                    actMap[locEdge-1][pos++] = 0;
                  }// check if edge has one non-zero node
		} // loop over dofs
              } // loop over functions
	    } // check if edge was already mapped
	  } // loop over all edges
	} // looop over all elements
      } // loop over all entitylists
    } // loop over all results
    
  } 

  void EqnMap::CalcFaceEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcFaceEquations" );

    LOG_TRACE(eqnMap) << "Starting to map face equations\n";
    
    // Big outer loop over all face mapped element lists
    ResultEntityMap::iterator listIt;
    
    for( listIt = edgeMappedList_.begin(); 
	 listIt != edgeMappedList_.end(); 
	 listIt++ ) {
      
      // Remeber current result and list of elementLists
      const ResultInfo & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      StdVector<Vector<Integer> > & actMap = faceEqns_[actRes];
      
      UInt dofsPerFace = actRes.dofNames.GetSize();
      actMap.Resize( numLocFaces_  );

      // Fetch also nodal map
      Matrix<Integer> & actNodeMap = nodeEqns_[actRes];

      Integer locFace = 0;

      // Iterate over all entitylists
      for( UInt iList = 0; iList < actLists.GetSize(); iList++ ) {

	EntityIterator it = actLists[iList]->GetIterator();
	
	// Iterate over all elements within this list
	for( it.Begin(); !it.IsEnd(); it++ ) {

	  // Get grip of element
	  const Elem & actEl = *(it.GetElem());
          
          // Get number of unknowns for each face
          StdVector<Vector<UInt> > numFcns; 
          numFcns.Resize( dofsPerFace );

          // iterate over all dofs of this result
          for( UInt iDof = 0; iDof < dofsPerFace; iDof++ ) {
            actEl.ptElem->GetNumFncs( numFcns[iDof], actRes.fctType, 
                                      AnsatzFct::FACE, iDof );
          }
          
          // Iterate over all faces of this element
          for( UInt iFace = 0; iFace < actEl.faces.GetSize(); iFace++ ) {
            
            // Check if there are any faces at all present
            //LOG_DBG3(eqnMap) << "Element has " << numFcns[iFace]
            //                 << "face functions!";
            //if( numFcns[iFace].GetSize() == 0 ) {
            //  continue;
            //}
            
            // determine local face number
            locFace = mesh2PdeFace_[actEl.faces[iFace] - 1];
            Face const & face = ptGrid_->GetFace( actEl.faces[iFace] );
            LOG_DBG3(eqnMap) << "Actual face: " << actEl.faces[iFace] << " / "
                             << locFace << " (global / local )";
            
            // Check if this face was already mapped
            if ( actMap[locFace-1].GetSize() == 0 && locFace > 0 ) {
              
              // sum up unknowns of this dof
              UInt sum = 0;
              UInt max = 0;
              for( UInt iDof = 0; iDof < dofsPerFace; iDof++ ) {
                sum +=  numFcns[iDof][iFace];
                if( numFcns[iDof][iFace] > max )
                  max = numFcns[iDof][iFace];
              }
              //LOG_DBG3(eqnMap) << "1-dof: " << numFcns[0].Serialize();
              //LOG_DBG3(eqnMap) << "2-dof: " << numFcns[1].Serialize();
              //LOG_DBG3(eqnMap) << "3-dof: " << numFcns[2].Serialize();
              LOG_DBG3(eqnMap) << "maximum for face " << iFace 
                               << ": " << max;
              LOG_DBG3(eqnMap) << "Face got " << sum << " equation numbers!\n";
             
              
              // iterate over all functions of this face
              UInt pos = 0;
              for( UInt iFcn = 0; iFcn < max; iFcn++ ) {
                
                if( actMap[locFace-1].GetSize() == 0 ) {
                  actMap[locFace-1].Resize( dofsPerFace * max);
                  actMap[locFace-1].Init(0);
                }
              
                // iterate over all dofs
                for( UInt iDof = 0; iDof < dofsPerFace; iDof++ ) {
                  
                  // Check if the related nodes have a 0 equations number
                  // (= are (in-)hom. Dirichlet nodes)
                  bool allZero = true;
                  for( UInt iNode = 0; iNode < face.nodes.GetSize(); iNode++ ) {
                    if( actNodeMap[mesh2PdeNode_[face.nodes[iNode]-1]-1][iDof] != 0) {
                      allZero = false;
                      break;
                    }
                  }
                  if( !allZero && (iFcn < numFcns[iDof][iFace]) ) {
                    actMap[locFace-1][pos++] = ++numEqns_;
                    LOG_DBG3(eqnMap) << "Face #" << actEl.faces[iFace]
                                     << " got equation number " << numEqns_-1;
                  } else {
                    actMap[locFace-1][pos++] = 0;
                    LOG_DBG3(eqnMap) << "-> Nodal equations 0 for dof" << iDof+1;
                  } // check for zero nodes
                } // loop over dofs
              }  // loop over functions
            } else {
              LOG_DBG3(eqnMap) << "--> Already mapped\n ";
            }// check if already mapped
          } // loop over faces
        } // loop over elements
      } // loop over results
    }
    
    LOG_TRACE(eqnMap) << "Finished mapping face equations\n";
    
  }

  void EqnMap::GetNodesOfEntities( StdVector<UInt>& nodes,
		                   shared_ptr<EntityList> ent ) {
    ENTER_FCN( "EqnMap::GetNodesOfEntities" );
    
    // Get type of entries of the particular entity list
    EntityList::ListType type = ent->GetType();
    
    shared_ptr<ElemList> elemList;
    shared_ptr<SurfElemList> sElemList;
    shared_ptr<NodeList> nodeList;
    shared_ptr<RegionList> regionList;
    StdVector<UInt> helpNodes;
    EntityIterator it;
    switch( type ) {
      
    case EntityList::ELEM_LIST:
      elemList= 
	dynamic_pointer_cast<ElemList, EntityList>(ent);
      ptGrid_->GetNodesByRegion( nodes, elemList->GetRegion() );
      break;
      
    case EntityList::SURF_ELEM_LIST:
      sElemList = 
	dynamic_pointer_cast<SurfElemList, EntityList>(ent);
      ptGrid_->GetNodesByRegion( nodes, sElemList->GetRegion() );
      break;

    case EntityList::REGION_LIST:
      regionList= 
	dynamic_pointer_cast<RegionList, EntityList>(ent);
      it = regionList->GetIterator();
      for( ; !it.IsEnd(); it++ ) {
        helpNodes.Clear();
        ptGrid_->GetNodesByRegion( helpNodes, it.GetRegion() );
        for( UInt i = 0; i < helpNodes.GetSize(); i++ ) {
          nodes.Push_back( helpNodes[i] );
        }
      }
      break;

    case EntityList::NODE_LIST:
      nodeList = 
	dynamic_pointer_cast<NodeList, EntityList>(ent);        
      nodes = nodeList->GetNodes();
      break;
      
    default :
      std::string listString;
      EntityList::Enum2String( type, listString );
      EXCEPTION( "'" << listString
                 << "' is no EntityList with nodal information." );
      
    }
  }

}
