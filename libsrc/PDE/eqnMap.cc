#include "eqnMap.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "Domain/domain.hh" 
#include "Utils/coordSystem.hh"

namespace CoupledField {


  EqnMap::EqnMap(Grid* grid, PdeIdType pdeId, bool sortEqns ) {
    ENTER_FCN( "EqnMap::EqnMap" );

    ptGrid_ = grid;
    sortEqns_ = sortEqns;
    pdeId_ = pdeId;
    isFinalized_ = false;


    numEqns_ = 0;
    numRealEqns_ = 0;
    numIdBcs_ = 0;
    
    numLocNodes_ = 0;
    numLocElems_ = 0;
    numLocEdges_ = 0;
    numLocSurfaces_ = 0;

  }

  EqnMap::~EqnMap() {
    ENTER_FCN( "EqnMap::~EqnMap" );

  }

  //! ======================================================================
  //! SET AND INITIALIZATION METHODS
  //! ======================================================================
  void EqnMap::AddResult( ResultDof& result, 
                          shared_ptr<EntityList> list ) {
    ENTER_FCN( "EqnMap::AddResult" );
    resEntMap_[result].Push_back( list );
  }

  void EqnMap::SetHomoDirichletBCs(HdBcList& hdBcs) {
    ENTER_FCN( "EqnMap::SetHomoDirichletBCs" );

    for( UInt i = 0; i< hdBcs.GetSize(); i++ ) {
      ResultDof actResult = *(hdBcs[i]->result);
      hdBcs_[actResult].Push_back(hdBcs[i]);
    }
  }
  
  void EqnMap::SetInhomDirichletBCs( IdBcList& idBcs ) {
    ENTER_FCN( "EqnMap::SetInhomDirichletBCs" );

    for( UInt i = 0; i < idBcs.GetSize(); i++ ) {
      ResultDof actResult = *(idBcs[i]->result);
      idBcs_[actResult].Push_back( idBcs[i] );
    }
  }
  
  void EqnMap::SetConstraints( ConstraintList& constraints ) {
    ENTER_FCN( "EqnMap::SetConstraints" );
    
    for( UInt i = 0; i < constraints.GetSize(); i++ ) {
      ResultDof actResult = *(constraints[i]->result);
      constraints_[actResult].Push_back( constraints[i] );
    }
  }
  
  void EqnMap::Finalize() {
    ENTER_FCN( "EqnMap::Finalize" );
    
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
          nodalLocElems_.Push_back( elemList );
        }
      }
    }

    // 2.) Get all entity lists, for which edges/surfaces have to be mapped
    // Get information about which region need mapping of
    // edges / surface 

    // ... to be implemented ...

    // 3.) Perform local<->global mapping of elements
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
      if( it->first.definedOn == ResultDof::NODE ) {
        for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
          nodeMappedList_[it->first].Push_back( it->second[iList] );
        }
      }
    }
    
    // assign equation numbers to nodes
    CalcNodalEquations( 1 );

    //   b) elem <-> eqnNr
    //   -----------------

    // iterate over all resultDofs
    for ( it=resEntMap_.begin(); it!=resEntMap_.end(); it++ ) {

      // check if resultDof is mapped onto nodes
      if( it->first.definedOn == ResultDof::ELEMENT ) {
        for( UInt iList = 0; iList < it->second.GetSize(); iList++ ) {
          elemMappedList_[it->first].Push_back( it->second[iList] );
        }
      }
    }
    
    // assign equation numbers to nodes
    CalcElemEquations( 1 );

    // ==== PHASE 2: Renumber fixed-equations ====

    // Afterwards re-iterate and give dirichletDOFs highest
    // numbers
    CalcNodalEquations( 2 );

    CalcElemEquations( 2 );

    // Now class is finalized
    isFinalized_ = true;
    
  }

  void EqnMap::ReorderMapping( Integer **order ) {
    ENTER_FCN( "EqnMap::ReorderMapping" );

    // Check if any reordering array was given
    if( (*order) == NULL ) {
      return;
    }

    // Iterate over all nodal maps and map the new ordering
    EqnMapType::iterator it;

    // iterate over all resuls
    for ( it=nodeEqns_.begin(); it!=nodeEqns_.end(); it++ ) {
      Matrix<Integer> & actMap = it->second;

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
    // Iterate over all element maps and map the new ordering
    for ( it=elemEqns_.begin(); it!=elemEqns_.end(); it++ ) {
      Matrix<Integer> & actMap = it->second;
      
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
    delete [] (*order);
    (*order) = NULL;
  }
  
  
  //! ======================================================================
  //! EQUATION MAPPING
  //! ======================================================================
  void EqnMap::GetEqns( StdVector<Integer> &eqns,
                        const ResultDof& result, 
                        const EntityIterator& it ) const{
    ENTER_IFCN( "EqnMap::GetEqns" );
    
    UInt numDofs = result.dofNames.GetSize();      
    Matrix<Integer> const & elemMap = (elemEqns_.find( result ) )->second;
    Matrix<Integer> const & map = (nodeEqns_.find( result ) )->second;
    

    if ( result.definedOn == ResultDof::NODE ) {
      
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
        return;
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
        (*error) << "This type of entity list is not defined for " 
                 << "equation mapping!";
        Error( __FILE__, __LINE__ );
      }

      return;
    } else if( result.definedOn == ResultDof::NODELIST ) {
      Error( "Not yet implemented", __FILE__, __LINE__ );
    }else if( result.definedOn == ResultDof::EDGE ) {

    } else if( result.definedOn == ResultDof::SURFACE ) {

    } else if( result.definedOn == ResultDof::REGION ) {
    
    } else if( result.definedOn == ResultDof::ELEMENT ) {
      eqns.Resize(numDofs);
      UInt localElem = mesh2PdeElem_[(it.GetElem()->elemNum)-1];
      for (UInt iDof = 0; iDof < numDofs; iDof++ ) {
        
        if (localElem < 1 ) {
          eqns[iDof] = 0;
        } else {
          eqns[iDof] = elemMap[localElem-1][iDof];
        }
        
      }
    } else {
      *error << "The entity type '" << result.definedOn
             << "' is not known for result types!";
    }
    
    
    // 1) NODE / NODELIST
    // -> Get entries from nodeEqns_ array
    // -> If it is higher order, fetch also entries from nodeEdgeEQNs_
    
    // 2) ELEM:
    // -> Get entries from elemEQNs_

    // 3) FREE:
    // -> Get entries from freeEQNs_
    
  }
  
  Integer EqnMap::GetEqn( const ResultDof& result, 
                          const EntityIterator& it,
                          UInt dof ) const{
    StdVector<Integer> eqns;
    GetEqns( eqns, result, it );
    return eqns[dof-1];
  }

    


  Integer EqnMap::GetNodeEqn( const ResultDof& result,
                              UInt nodeNr, UInt dof ) { 
    ENTER_FCN( "EqnMap::GetNodeEqn" );
    

    if ( result.definedOn == ResultDof::NODE ) {
      
      
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
      (*error) << "GetNodeEqn() can only be used if exactly "
               << "one nodal result is present!";
      Error( __FILE__, __LINE__ );
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
      (*error) << "GetNodeEqn() can only be used if exactly "
               << "one nodal result is present!";
      Error( __FILE__, __LINE__ );
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
      (*error) << "MeshNode Nr. " << meshNode 
               << " has no local node number!";
      Error( __FILE__, __LINE__);
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
      (*error) << "MeshElem Nr. " << elemNumGlob 
               << " has no local elem number!";
      Error( __FILE__, __LINE__);
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

    // --- Local<->Global mapping (edges) ---
    out << "3.) local <-> global mapping (edges)\n"
        << "------------------------------------\n";

    // --- Nodal equation mapping ---
    out << "4.) Nodal equation mapping\n"
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
      for ( UInt iList = 0; iList < nodalLocElems_.GetSize(); iList++ ) {

        // Get iterator of current element list
        EntityIterator it = nodalLocElems_[iList]->GetIterator();
        
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
  
  
  void EqnMap::CalcEdgeSurfaceMapping()  {
    ENTER_FCN( "EqnMap::CalcEdgeSurfaceMapping" );
  }
  
  
  void EqnMap::CalcNodalEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcNodalEquations" );


    // Big outer loop over all nodal mapped element lists
    ResultEntityMap::iterator listIt;

    for( listIt = nodeMappedList_.begin(); 
         listIt != nodeMappedList_.end(); 
         listIt++ ) {

      UInt eqnCounter = numRealEqns_;
      
      // Remeber current result and list of elementLists
      const ResultDof & actRes = listIt->first;
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
        
        actMap.Resize( numLocNodes_, dofsPerNode );
        actMap.Init( 1 );
        
        
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
                (*warning) << "EqnMap::CalcNodalEquations: HomDirichletNode # "
                           << nodes[i]
                           << "\nappeared already at least once in the list of "
                           << "boundary nodes for this PDE!\n Please check, if this "
                           << "node is defined in more than one level of boundary "
                           << "nodes!";
                Warning( __FILE__, __LINE__ );
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
                (*warning) << "EqnMap::CalcNodalEquations: Inhom. Dirichlet "
                           << "node #" << nodes[iNode]
                           << "\nappeared already at least once in the list of "
                           << "boundary nodes for this Pde!\n Please check, if "
                           << "this node is defined in more than one level of "
                           << "boundary nodes!";
                Warning( __FILE__, __LINE__ );
              }
              else {
                if ( sortEqns_ == true ) {
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
            
            for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
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
                eqnCounter++;
                actMap[locNode-1][iDof] = eqnCounter;
                countNodes[locNode-1][iDof] = 1;
              }
            }
          }
        }
        
        // now we know the number of 'real' dofs
        numRealEqns_ = eqnCounter;

      } else if( phase == 2 ) {
        
        // ------
        // STEP 5
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
            for ( UInt iNode = 1; iNode < slaveNodes.GetSize(); iNode++ ) {
              actMap[mesh2PdeNode_[slaveNodes[iNode]-1]-1] [slaveDof-1] =
                -actMap[mesh2PdeNode_[masterNode-1]-1] [masterDof-1];
            }
          }
        }
        
        // -------
        // STEP 5b
        // -------
        if ( sortEqns_ == true ) {
          if( idBcIt != idBcs_.end() ) {
            IdBcList const & actIdBcList = idBcIt->second;
            for ( UInt i = 0; i < actIdBcList.GetSize(); i++ ) {
              StdVector<UInt> nodes;
              GetNodesOfEntities( nodes, actIdBcList[i]->entities );
              UInt actDof = actIdBcList[i]->dof;
              
              for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
                eqnCounter++;
                actMap[mesh2PdeNode_[nodes[iNode]-1]-1] [actDof-1] = eqnCounter;
              }
            }
          }
        }
        
        
        numEqns_ = eqnCounter;        
        
        //numDroppedDofs_ = numLocNodes_ * dofsPerNode - numEqns_ + multipleBCs;
 
      } else {
        *error << "Phase '" << phase << "' does not exist!";
        Error( __FILE__, __LINE__ );
      }
    }
  }
   
  void EqnMap::CalcElemEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcElemEquations" );

     // Big outer loop over all nodal mapped element lists
    ResultEntityMap::iterator listIt;
    
    for( listIt = elemMappedList_.begin(); 
         listIt != elemMappedList_.end(); 
         listIt++ ) {

      UInt eqnCounter = numRealEqns_;
      
      // Remeber current result and list of elementLists
      const ResultDof & actRes = listIt->first;
      StdVector<shared_ptr<EntityList> > & actLists = listIt->second;
      
      
      // Get grip of homogeneous and in-homogeneous boundary conditions
      // for this tpye of result
      ResultHdBcMap::iterator hdBcIt = hdBcs_.find( actRes );
      ResultIdBcMap::iterator idBcIt = idBcs_.find( actRes );
      ResultConstraintMap::iterator csIt = constraints_.find( actRes );
      
      
      std::string type;
      Enum2String( actRes.resultType, type );
      Matrix<Integer> & actMap = elemEqns_[actRes];

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
        
        actMap.Resize( numLocElems_, dofsPerElem );
        actMap.Init( 1 );

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
            for ( UInt iDof = 0; iDof < dofsPerElem; iDof++ ) {
              if ( actMap[locElem-1][iDof] != 0 &&
                   countElems[locElem-1][iDof] == 0) {
                eqnCounter++;
                actMap[locElem-1][iDof] = eqnCounter;
                countElems[locElem-1][iDof] = 1;
              }
            }
          }
        }
        
        // now we know the number of 'real' dofs
        numRealEqns_ = eqnCounter;

      } else if( phase == 2 ) {
        
        // ------
        
        
        
        numEqns_ = eqnCounter;        
        
        //numDroppedDofs_ = numLocNodes_ * dofsPerNode - numEqns_ + multipleBCs;
        
      } else {
        *error << "Phase '" << phase << "' does not exist!";
        Error( __FILE__, __LINE__ );
      }
    }

  }
  
  void EqnMap::CalcEdgeSurfEquations( UInt phase ) {
    ENTER_FCN( "EqnMap::CalcEdgeSurfEquations" );
  }

  void EqnMap::GetNodesOfEntities( StdVector<UInt>& nodes,
                                   shared_ptr<EntityList> ent ) {
    ENTER_FCN( "EqnMap::GetNodesOfEntities" );
    
    // Get type of entries of the particular entity list
    EntityList::Type type = ent->GetType();
    
    shared_ptr<ElemList> elemList;
    shared_ptr<SurfElemList> sElemList;
    shared_ptr<NodeList> nodeList;
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

    case EntityList::NODE_LIST:
      nodeList = 
        dynamic_pointer_cast<NodeList, EntityList>(ent);        
      ptGrid_->GetNodesByName( nodes, nodeList->GetNodesName() );
      break;
      
    default :
      *error <<"'" << EntityList::TypeToString(type)
             << "' is no EntityList with nodal information.";
      Error( __FILE__, __LINE__ );
      
    }
  }

}
