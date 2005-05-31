#include <iomanip>

#include "scalarnodeEQN.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"


namespace CoupledField {

  
  // ===============
  //   Constructor
  // ===============
  ScalarNodeEQN::ScalarNodeEQN(Grid * aptGrid, 
                               StdVector<RegionIdType> & asubdoms, 
                               UInt dofsPerNode)
    : NodeEQN(aptGrid, asubdoms, dofsPerNode) {

    ENTER_FCN( "NodeEQN::NodeEQN" );
  
    isBlockMapped_ = FALSE;
    dofsPerEQN_ = 1;
 
  }


  // ==============
  //   Destructor
  // ==============
  ScalarNodeEQN::~ScalarNodeEQN() {
    ENTER_FCN( "ScalarNodeEQN::~ScalarNodeEQN" );
    if ( commandLine->GetShowEqnMap() == true ) {
      Print( (*cla) );
    }
  }


  // ===============
  //   CalcMapping
  // ===============
  void ScalarNodeEQN::CalcMapping() {

    ENTER_FCN( "ScalarNodeEQN::CalcMapping" );
 
    // First apply Mapping from global to
    // local node numbers and back
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
    //          and set the corresponding entry in pdeNode2EQN_ to the negative
    //          of the value of constraintMasterNode

    // Each local PDE-node number is assigned a 
    // corresponding EQN number:
    //
    // 0      : homog. Dirichlet BC
    // +number: normal Equation
    // -number: slave node (constraint)
    //
    // Therefore both pdeNode2EQN has the same size
    // as number of PDE nodes.
    // Furthermore, eqn2Pos has the length of
    // total number of equations, which is the total
    // number of pdenodes minus the number of hom.
    // Dirichlet nodes - the number of constraintNodes


    // **********
    //   STEP 1
    // **********
    UInt notIncludedBCs = 0;
    UInt multipleBCs = 0;
    pdeNode2EQN_.Clear();
    pdeNode2EQN_.Resize(numPDENodes_);
    pdeNode2EQN_.Init(1);
 
    // **********
    //   STEP 2
    // **********
    StdVector<UInt> countNodes;
    countNodes.Resize(numPDENodes_);
    for ( UInt i = 0; i < homoDirichletNodes_.GetSize(); i++ ) {
      if (mesh2PDENode_[homoDirichletNodes_[i]-1]-1 < 0) {
        (*warning) << "ScalarNodeEQN::CalcMapping: Homogen. Dirichlet node #"
                   << homoDirichletNodes_[i]
                   << " is not contained in any of the regions for this PDE";
        Warning( __FILE__, __LINE__ );
        notIncludedBCs++;
      }
      else if (countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1] != 0) {
        (*warning) << "ScalarNodeEQN::CalcMapping: HomDirichletNode # "
                 << homoDirichletNodes_[i]
                 << "\nappeared already at least once in the list of "
                 << "boundary nodes for this PDE!\n Please check, if this "
                 << "node is defined in more than one level of boundary "
                 << "nodes!";
        Warning( __FILE__, __LINE__ );
        multipleBCs++;
      }
      else {
        pdeNode2EQN_[mesh2PDENode_[homoDirichletNodes_[i]-1]-1] = 0;
        countNodes[mesh2PDENode_[homoDirichletNodes_[i]-1]-1]++;
      }
    }
 
    // ***********
    //   STEP 2b
    // ***********
    if ( sortEQNs_ == TRUE ) {
      countNodes.Init(0);
      for ( UInt i = 0; i < inhomDirichletNodes_.GetSize(); i++ ) {
        if (mesh2PDENode_[inhomDirichletNodes_[i]-1]-1 < 0) {
          (*warning) << "ScalarNodeEQN::CalcMapping: Inhom. Dirichlet node #"
                     << inhomDirichletNodes_[i]
                     << " is not contained in any of the regions for this PDE";
          Warning( __FILE__, __LINE__ );
          notIncludedBCs++;
        }
        else if (countNodes[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1] != 0) {
          (*error) << "ScalarNodeEQN::CalcMapping: InhomDirichletNode # "
                   << inhomDirichletNodes_[i]
                   << "\nappeared already at least once in the list of "
                   << "boundary nodes for this PDE!\n Please check, if this "
                   << "node is defined in more than one level of boundary "
                   << "nodes!";
          Warning( __FILE__, __LINE__ );
          multipleBCs++;
        }
        else {
          pdeNode2EQN_[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1] = 0;
          countNodes[mesh2PDENode_[inhomDirichletNodes_[i]-1]-1]++;
        }
      }
    }
  
    // **********
    //   STEP 3
    // **********
    for ( UInt i = 0; i < constraintSlaveNodes_.GetSize(); i++ ) {
      pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] = 0;
    }

    // **********
    //   STEP 4
    // **********
    UInt eqnCounter = 0;
    for ( UInt i = 0; i < pde2MeshNode_.GetSize(); i++ ) {
      if ( pdeNode2EQN_[i] != 0 ) {
        eqnCounter++;
        pdeNode2EQN_[i] = eqnCounter;
      }
    }


    // Now we know how many "real" equations we have
    if ( sortEQNs_ == TRUE ) {
      numRealEqns_ = eqnCounter;
    }
    else {
      numRealEqns_ = 0;
    }


    // **********
    //   STEP 5
    // **********
    for ( UInt i = 0; i < constraintSlaveNodes_.GetSize(); i++ ) {
      pdeNode2EQN_[mesh2PDENode_[constraintSlaveNodes_[i]-1]-1] =
        - pdeNode2EQN_[mesh2PDENode_[constraintMasterNodes_[i]-1]-1];
    }

    // ***********
    //   STEP 5b
    // ***********
    if ( sortEQNs_ == TRUE ) {
      for ( UInt i = 0; i < inhomDirichletNodes_.GetSize(); i++ ) {
        eqnCounter++;
        pdeNode2EQN_[i] = eqnCounter;
      }
    }

    // Now the object is initialized
    isInitialized_ = TRUE;
    numEqns_ = eqnCounter;

    numBuildInDirichletEQNs_ = numPDENodes_ - numEqns_ + multipleBCs;
  }


  // ====================
  //   CalcMpcciMapping
  // ====================
  void ScalarNodeEQN::CalcMpcciMapping() {

    ENTER_FCN( "ScalarNodeEQN::CalcMpcciMapping" );

    // Apply Mapping from global to
    // local node numbers and back
    CalcLocalGlobalMapping( mesh2PDENode_, pde2MeshNode_, mesh2PDEElem_,
                            pde2MeshElem_ );
  }


  // =========
  //   Print
  // =========
  void ScalarNodeEQN::Print(std::ostream & out) const {

    ENTER_FCN( "ScalarNodeEQN::Print" );
  
    out << "Equation numbering - Information" << std::endl;
    out << "================================" << std::endl;
    out << "DOFs per Node: 1 " << std::endl;
    out << "Using SCALAR numbering of equations" << std::endl;
    out << std::endl;

    // Print pde2MeshNode_ and pdeNode2EQN_
    out << std::setw(10) << "PDE NodeNr" << " | ";
    out << std::setw(13) << "Global NodeNr" << " | ";
    out << std::setw(10) << "EQ number";
    out << std::endl;
    out << std::setfill('-') << std::setw(39) << "-" << std::endl;
    out << std::setfill(' ');
  
    for ( UInt i = 0; i < pde2MeshNode_.GetSize(); i++ ) {
      out << std::setw(10) << i+1  << " | ";
      out << std::setw(13) << pde2MeshNode_[i] << " | ";
      out << std::setw(10) << pdeNode2EQN_[i] << std::endl;
    }
  }


  // ============
  //   Node2EQN
  // ============
  void ScalarNodeEQN::Node2EQN( const UInt nodeNr, 
                                const UInt dof,
                                Integer &eqnNr,
                                UInt &eqnDof ) const { 

    ENTER_FCN( "ScalarNodeEQN::Node2EQN" );

#ifdef CHECK_INDEX
    if ( nodeNr > mesh2PDENode_.GetSize() ) {
      (*error) << "ScalarNodeEQN::Node2EQN: Index out of bounds"
               << nodeNr << " > " << mesh2PDENode_.GetSize();
      Error( __FILE__, __LINE__ );
    }
#endif

    eqnNr = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
    eqnDof = 1;
  }


  // ============
  //   Node2EQN
  // ============
  void ScalarNodeEQN::Node2EQN( const UInt nodeNr,
                                StdVector<Integer> &eqnNr ) const {

    ENTER_FCN( "ScalarNodeEQN::Node2EQN" );

#ifdef CHECK_INDEX
    if ( nodeNr > mesh2PDENode_.GetSize() ) {
      (*error) << "ScalarNodeEQN::Node2EQN: Index out of bounds"
               << nodeNr << " > " << mesh2PDENode_.GetSize();
      Error( __FILE__, __LINE__ );
    }
#endif 

    eqnNr.Resize(dofsPerNode_);

    eqnNr[0] = pdeNode2EQN_[mesh2PDENode_[nodeNr-1]-1];
  }


  // ============
  //   Node2EQN
  // ============
  void ScalarNodeEQN::Node2EQN( const StdVector<UInt> &nodeNr,
                                StdVector<Integer> &eqnNr ) const {

    ENTER_FCN( "ScalarNodeEQN::Node2EQN" );
  
    eqnNr.Resize(nodeNr.GetSize());

    for ( UInt i = 0; i < nodeNr.GetSize(); i++ ) {

#ifdef CHECK_INDEX
      if ( nodeNr[i] > mesh2PDENode_.GetSize() ) {
        (*error) << "ScalarNodeEQN::Node2EQN: Index " << i << " out of bounds"
                 << nodeNr[i] << " > " << mesh2PDENode_.GetSize();
        Error( __FILE__, __LINE__ );
      }
#endif

      eqnNr[i] =  pdeNode2EQN_[mesh2PDENode_[nodeNr[i]-1]-1];
    }
  }


  // ================================================
  //   Equation re-mapping according to re-ordering
  // ================================================
  void ScalarNodeEQN::ReorderMapping(Integer *order) {

    ENTER_FCN( "ScalarNodeEQN::ReorderMapping" );

    if ( order != NULL ) {
      for ( UInt i = 0; i < pdeNode2EQN_.GetSize(); i++ ) {
        if ( pdeNode2EQN_[i] > 0 ) {
          pdeNode2EQN_[i] = order[pdeNode2EQN_[i]-1];
        }
        else if(pdeNode2EQN_[i] < 0 ) {
          //due to constraints
          pdeNode2EQN_[i] = -order[-pdeNode2EQN_[i]-1];   
        }
      }
    }
  }


} // end of namespace
