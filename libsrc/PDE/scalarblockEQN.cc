#include <iomanip>

#include "scalarblockEQN.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"


namespace CoupledField {
  

  // ***************
  //   Constructor
  // ***************
  ScalarBlockEQN::ScalarBlockEQN( Grid *aptGrid, 
                                  StdVector<RegionIdType> &asubdoms, 
                                  UInt dofsPerNode, bool sortEQNs )
    : NodeEQN( aptGrid, asubdoms, dofsPerNode, sortEQNs ) {
    ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
    isBlockMapped_ = FALSE;
    dofsPerEQN_    = 1;
  }


  // **************
  //   Destructor
  // **************
  ScalarBlockEQN::~ScalarBlockEQN() {
    ENTER_FCN( "ScalarBlockEQN::ScalarBlockEQN" );
  }


  // ***************
  //   CalcMapping
  // ***************
  void ScalarBlockEQN::CalcMapping() {

    ENTER_FCN( "ScalarBlockEQN::CalcMapping" );

    // First apply Mapping from global to
    // local node/elem numbers and back
    CalcLocalGlobalMapping(mesh2PDENode_,
                           pde2MeshNode_,
                           mesh2PDEElem_,
                           pde2MeshElem_);

 
    // Idea of the algorithm:
    //
    // Step 1:  Initialize pdeNode2eqn_ with 1
    // Step 2:  For each entry in homoDirichletNodes_ set the corresponding
    //          entry in pdeNode2eqn_ to 0
    // Step 2b: For each entry in inhomoDirichletNodes_ set the corresponding
    //          entry in pdeNode2eqn_ to 0
    // Step 3:  For each entry in constraintSlaveNodes_ set the corresponding
    //          entry in pdeNode2eqn_ to 0
    // Step 4:  Loop over all entries in pde2Meshnode
    //          and assign each non-zero entry an equation number
    // Step 5:  Afterwards loop again over all nodes in constraintSlaveNodes_
    //          and set the corresponding entry in pdeNode2EQN_ to the
    //          negative of the value of constraintMasterNode
    // Step 5b: Loop over all entries in inhomoDirichletNodes_ and assign that
    //          dof an equation number after the hightest equation number of
    //          the free dofs
    //
    // Note:    Steps 2b and 5b are only performed, if sortEQNs = TRUE


    // ------
    // STEP 1
    // ------
    UInt multipleBCs = 0;
    pdeNode2EQN_.Resize( numPDENodes_, dofsPerNode_ );
    pdeNode2EQN_.Init( 1 );


    // ------
    // STEP 2
    // ------
    Matrix<UInt> countNodes;
    countNodes.Resize( numPDENodes_, dofsPerNode_ );
    countNodes.Init( 0 );
    for ( UInt i = 0; i < homoDirichletNodes_.GetSize(); i++ ) {

      // Check if homDirichletNode belongs to one of my subdomains
      if ( mesh2PDENode_[homoDirichletNodes_[i]-1]-1 < 0 ) {
        (*warning) << "ScalarBlockEQN::CalcMapping: Homogen. Dirichlet node "
                   << "nr. " << homoDirichletNodes_[i]
                   << " is not contained in any of the regions for this PDE";
        Warning( __FILE__, __LINE__ );
        multipleBCs++;
      }
      else if ( countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
                [homoDirichletDofs_[i]-1] != 0 ) {
        (*warning) << "ScalarBlockEQN::CalcMapping: HomDirichletNode # "
                 << homoDirichletNodes_[i]
                 << "\nappeared already at least once in the list of "
                 << "boundary nodes for this PDE!\n Please check, if this "
                 << "node is defined in more than one level of boundary "
                 << "nodes!";
        Warning( __FILE__, __LINE__ );
        multipleBCs++;
      }
      else {
        pdeNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
          [homoDirichletDofs_[i]-1] = 0;
        countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]
          [homoDirichletDofs_[i]-1]++;
      }
    }


    // -------
    // STEP 2b
    // -------
    if ( sortEQNs_ == TRUE ) {
      countNodes.Init(0);
      for ( UInt i = 0; i < inhomDirichletNodes_.GetSize(); i++ ) {
        if ( mesh2PDENode_[ inhomDirichletNodes_[i] - 1 ] - 1 < 0 ) {
          (*warning) << "ScalarBlockEQN::CalcMapping: Inhom. Dirichlet "
                     << "node #"
                     << inhomDirichletNodes_[i]
                     << " is not contained in any of the regions for "
                     << "this PDE";
          Warning( __FILE__, __LINE__ );
          multipleBCs++;
        }
        else if ( countNodes[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1]
                  [inhomDirichletDofs_[i]-1] != 0 ) {
          (*warning) << "ScalarBlockEQN::CalcMapping: Inhom. Dirichlet "
                     << "node #"
                     << inhomDirichletNodes_[i]
                     << "\nappeared already at least once in the list of "
                     << "boundary nodes for this PDE!\n Please check, if "
                     << "this node is defined in more than one level of "
                     << "boundary nodes!";
          Warning( __FILE__, __LINE__ );
          multipleBCs++;
        }
        else {
          pdeNode2EQN_[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1]
            [inhomDirichletDofs_[i]-1] = 0;
          countNodes[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1]
            [inhomDirichletDofs_[i]-1]++;
        }
      }
    }


    // ------
    // STEP 3
    // ------
    for ( UInt i = 0; i < constraintSlaveNodes_.GetSize(); i++ ) {
      pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]
        [constraintDofs_[i]-1] = 0;
    }

  
    // ------
    // STEP 4
    // ------
    UInt eqnCounter = 0;
    for ( UInt iNode = 0; iNode < pde2MeshNode_.GetSize(); iNode++ ) {
      for ( UInt iDof = 0; iDof < dofsPerNode_; iDof++ ) {
        if ( pdeNode2EQN_[iNode][iDof] != 0 ) {
          eqnCounter++;
          pdeNode2EQN_[iNode][iDof] = eqnCounter;
        }
      }
    }


    // now we know the number of 'real' dofs
    numRealEqns_ = eqnCounter;


    // ------
    // STEP 5
    // ------
    for ( UInt i = 0; i < constraintSlaveNodes_.GetSize(); i++ ) {
      pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1]
        [constraintDofs_[i]-1] =
        pdeNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1]
        [constraintDofs_[i]-1];
    }


    // -------
    // STEP 5b
    // -------
    if ( sortEQNs_ == TRUE ) {
      for ( UInt i = 0; i < inhomDirichletNodes_.GetSize(); i++ ) {
        eqnCounter++;
        pdeNode2EQN_[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1]
          [inhomDirichletDofs_[i]-1] = eqnCounter;
      }
    }


    // Now object is initialized
    isInitialized_ = TRUE;
    numEqns_ = eqnCounter;
    numDroppedDofs_ = numPDENodes_ * dofsPerNode_ - numEqns_ + multipleBCs;


    // Print information to LAS file
    if ( commandLine->GetShowEqnMap() == true ) {
      Print( (*cla) );
    }

  }


  // *********
  //   Print
  // *********
  void ScalarBlockEQN::Print( std::ostream &out ) const {

    ENTER_FCN( "ScalarBlockEQN::Print" );

    out << "Equation numbering - Information\n"
        << "================================\n"
        << "DOFs per Node: " << dofsPerNode_ << '\n'
        << "Using SCALAR-BLOCK numbering of equations\n\n"
        << "highest equation numbers             = " << numEqns_ << '\n'
        << "number of free degrees of freedom    = " << numRealEqns_ << '\n'
        << "number of dropped degrees of freedom = " << numDroppedDofs_
        << '\n' << std::endl;

    // Print pde2MeshNode_ and pdeNode2EQN_
    out << std::setw(10) << "PDE NodeNr" << " | ";
    out << std::setw(13) << "Global NodeNr" << " | ";
    out << std::setw(10) << "EQ number";
    out << std::endl;
    out << std::setfill('-') << std::setw(39) << "-" << std::endl;
    out << std::setfill(' ');

    for ( UInt iNode = 0; iNode < pde2MeshNode_.GetSize(); iNode++ ) {

      // Print out first dof
      out << std::setw(10) << iNode+1  << " | ";
      out << std::setw(13) << pde2MeshNode_[iNode] << " | ";
      out << std::setw(10) << pdeNode2EQN_[iNode][0] << std::endl;
      
      // Print out further dofs if needed
      for ( UInt iDof = 1; iDof < dofsPerNode_; iDof++ ) {
        out << std::setw(10) << " " << " | ";
        out << std::setw(13) << " " << " | ";
        out << std::setw(10) << pdeNode2EQN_[iNode][iDof] << std::endl;
      }
    }
  }


  // ************************
  //   Node2Eqn (Version 1)
  // ************************
  void ScalarBlockEQN::Node2EQN( const UInt nodeNr, 
                                 const UInt dof,
                                 Integer &eqnNr,
                                 UInt &eqnDof) const {

    ENTER_FCN( "ScalarBlockEQN::Node2EQN" );

#ifdef CHECK_INDEX
    if ( nodeNr > mesh2PDENode_.GetSize() ) {
      Error( "ScalarNodeEQN::Node2EQN: Index out of bounds",
             __FILE__, __LINE__ );
    }
#endif

    Integer localNode = mesh2PDENode_[nodeNr-1];
    if (localNode < 1 ) {
      eqnNr = 0;
    } else {
      eqnNr = pdeNode2EQN_[localNode-1][dof-1];
    }
    eqnDof = 1;
  }


  // ************************
  //   Node2Eqn (Version 2)
  // ************************
  void ScalarBlockEQN::Node2EQN( const UInt nodeNr,
                                 StdVector<Integer> &eqns ) const {

    ENTER_FCN( "ScalarBlockEQN::Node2EQN" );

#ifdef CHECK_INDEX
    if ( nodeNr > mesh2PDENode_.GetSize() ) {
      Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
            __FILE__, __LINE__);
    }
#endif

    eqns.Resize(dofsPerNode_);
    Integer localNode = 0;

    for ( UInt i = 0; i < dofsPerNode_; i++ ) {
      localNode = mesh2PDENode_[nodeNr-1];
      if (localNode < 1 ) {
        eqns[i] = 0;
      } else {
        eqns[i] = pdeNode2EQN_[localNode-1][i];
      }
    }
  }



  // ************************
  //   Node2Eqn (Version 3)
  // ************************
  void ScalarBlockEQN::Node2EQN( const StdVector<UInt> &nodeNr,
                                 StdVector<Integer> &eqnNr ) const {

    ENTER_FCN( "ScalarBlockEQN::Node2EQN" );
    Integer localNode = 0;
    eqnNr.Resize(nodeNr.GetSize()*dofsPerNode_);

    for( UInt iNode = 0; iNode < nodeNr.GetSize(); iNode++ ) {
      
#ifdef CHECK_INDEX
      if ( nodeNr[iNode] > mesh2PDENode_.GetSize() ) {
        Error("ScalarNodeEQN::Node2EQN: Index out of bounds", 
              __FILE__, __LINE__);
      }
#endif
      for (UInt iDof = 0; iDof < dofsPerNode_; iDof++ ) {
        localNode = mesh2PDENode_[nodeNr[iNode]-1];
        if (localNode < 1 ) {
          eqnNr[iNode*dofsPerNode_ + iDof] = 0;
        } else {
          eqnNr[iNode*dofsPerNode_ + iDof] = pdeNode2EQN_[localNode-1][iDof];
        }
      }
    }
  }
    

  // ******************
  //   ReorderMapping
  // ******************
  
  // Equation re-mapping according to re-ordering
  void ScalarBlockEQN::ReorderMapping( Integer **order ) {
    
    ENTER_FCN( "ScalarBlockEQN::ReorderMapping" );
    
    if ( (*order) != NULL ) {
      for ( UInt i = 0; i < pdeNode2EQN_.GetSizeRow(); i++ ) {
        for ( UInt j = 0; j < pdeNode2EQN_.GetSizeCol(); j++ ) {
          if ( pdeNode2EQN_[i][j] > 0 ) {
            pdeNode2EQN_[i][j] = (*order)[pdeNode2EQN_[i][j]-1];
          }
          else if(pdeNode2EQN_[i][j] < 0 ) {
            //due to constraints
            pdeNode2EQN_[i][j] = -(*order)[-pdeNode2EQN_[i][j]-1];   
          }
        }
      }
      delete [] (*order);
      (*order) = NULL;
    }
  }
  
} // end of namespace
