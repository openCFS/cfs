#include "baseEQN.hh"
#include <list>
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  BaseEQN::BaseEQN( Grid * aptGrid, 
		    BCs * aptBCs,
		    StdVector<std::string>& asubdoms, 
		    Integer actlevel, 
		    Integer dofsPerNode ) {

    ENTER_FCN( "BaseEQN::BaseEQN" );

    ptGrid_   = aptGrid;
    ptBCs_    = aptBCs;
    subdoms_  = asubdoms;
    actlevel_ = actlevel;
    dofsPerNode_ = dofsPerNode;

    isInitialized_ = FALSE;
    sortEQNs_ = FALSE;
    numRealEqns_ = -1;
  }


  // **************
  //   Destructor
  // **************
  BaseEQN::~BaseEQN() {
    ENTER_FCN( "~BaseEQN::~BaseEQN" );
  }


  // ***********************
  //   SetHomDirichletBCs
  // ***********************
  void BaseEQN::SetHomoDirichletBCs( const StdVector<std::string> &nodeLevel,
				     const StdVector<std::string> &dofs ) {

    ENTER_FCN( "BaseEQN::SetHomoDirichletBCs" );

    std::list<Integer> tempNodeList;
    homoDirichletNodes_.Clear();

    for ( Integer i = 0; i < nodeLevel.GetSize(); i++ ) {
      
      tempNodeList = ptBCs_->GetNodesLevel(nodeLevel[i]);
      
      std::list<Integer>::iterator it;
      
      for ( it = tempNodeList.begin(); it != tempNodeList.end(); it++ ) {
	homoDirichletNodes_.Push_back(*it);
	if (dofsPerNode_ > 1) {
	  homoDirichletDofs_.Push_back(GetBCDof(dofs[i]));
	}
      } 
    }
  }


  // ************************
  //   SetInhomDirichletBCs
  // ************************
  void BaseEQN::SetInhomDirichletBCs( const StdVector<std::string> &nodeLevel,
				      const StdVector<std::string> &dofs ) {

    ENTER_FCN( "BaseEQN::SetInhomDirichletBCs" );

    // Only do this, if we have to sort the equation numbers
    // with the inhom. Dirichlet part on the top
    if ( sortEQNs_ == TRUE ) {
      std::list<Integer> tempNodeList;
      inhomDirichletNodes_.Clear();

      for ( Integer i = 0; i < nodeLevel.GetSize(); i++ ) {
      
	tempNodeList = ptBCs_->GetNodesLevel(nodeLevel[i]);
	std::list<Integer>::iterator it;
      
	for ( it = tempNodeList.begin(); it != tempNodeList.end(); it++ ) {
	  inhomDirichletNodes_.Push_back(*it);
	  if (dofsPerNode_ > 1) {
	    inhomDirichletDofs_.Push_back(GetBCDof(dofs[i]));
	  }
	} 
      }
    }
  }


  // ******************
  //   SetConstraints
  // ******************
  void BaseEQN::SetConstraints( const StdVector<Integer> &slaveNodeNrs,
				const StdVector<Integer> &masterNodeNrs,
				const StdVector<std::string> &dofs ) {

    ENTER_FCN( "BaseEQN::SetConstraints" );

    constraintSlaveNodes_ = slaveNodeNrs;
    constraintMasterNodes_ = masterNodeNrs;

    if ( dofsPerNode_ > 1 ) {
      constraintDofs_.Resize(dofs.GetSize());

      for ( Integer i = 0; i < dofs.GetSize(); i++ ) {
	constraintDofs_[i] = GetBCDof(dofs[i]);
      }
    }
  }


  // *************
  //   GetBCDofs
  // *************
  Integer BaseEQN::GetBCDof( const std::string dofString ) const {

    ENTER_FCN( "BaseEQN::GetBCDof" );
  
    Integer retVal = 0;
  
    if ( dofString == "ux" ) {
      retVal = 1;
    }
    else if ( dofString == "uy" ) {
      retVal = 2;
    }
    else if ( dofString == "uz" ) {
      retVal = 3;
    }
  
    // hard-coded for Piezo PDE
    else if ( dofString == "ep" ) {
      retVal = dofsPerNode_;
    }

    else {
      // According to the Schema definition of the parameter file this cannot
      // happen. Did the parser not perform validation?
      (*error) << "Direction should be one of ux, uy, uz and not "
	       << dofString;
      Error( __FILE__, __LINE__ );
    }

    return retVal;
  }

} // end of namespace
