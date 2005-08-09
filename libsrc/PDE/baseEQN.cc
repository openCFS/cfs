#include "baseEQN.hh"
#include <list>
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  BaseEQN::BaseEQN( Grid *aptGrid, StdVector<RegionIdType>& asubdoms, 
                    UInt dofsPerNode, bool sortEQNs ) {

    ENTER_FCN( "BaseEQN::BaseEQN" );

    ptGrid_        = aptGrid;
    subdoms_       = asubdoms;
    dofsPerNode_   = dofsPerNode;
    isInitialized_ = FALSE;
    numRealEqns_   = 0;
    numEqns_       = 0;
    sortEQNs_      = sortEQNs;
  }


  // **************
  //   Destructor
  // **************
  BaseEQN::~BaseEQN() {
    ENTER_FCN( "BaseEQN::~BaseEQN" );
  }


  // ***********************
  //   SetHomDirichletBCs
  // ***********************
  void BaseEQN::SetHomoDirichletBCs( const StdVector<std::string> &nodeLevel,
                                     const StdVector<std::string> &dofs ) {

    ENTER_FCN( "BaseEQN::SetHomoDirichletBCs" );

    StdVector<UInt> tempNodeList;
    UInt dof;
    homoDirichletNodes_.Clear();

    for ( UInt i = 0; i < nodeLevel.GetSize(); i++ ) {
      
      ptGrid_->GetNodesByName( tempNodeList, nodeLevel[i]);
      
      for ( UInt iNode = 0; iNode < tempNodeList.GetSize(); iNode++ ) {
        homoDirichletNodes_.Push_back(tempNodeList[iNode]);
        if (dofsPerNode_ > 1) {
          dof = domain->GetCoordSystem()->GetVecComponent(dofs[i]);
          homoDirichletDofs_.Push_back(dof);
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

    UInt dof = 0;

    // Only do this, if we have to sort the equation numbers
    // with the inhom. Dirichlet part on the top
    if ( sortEQNs_ == TRUE ) {
      StdVector<UInt> tempNodeList;
      inhomDirichletNodes_.Clear();

      for ( UInt i = 0; i < nodeLevel.GetSize(); i++ ) {
      
        ptGrid_->GetNodesByName( tempNodeList, nodeLevel[i] );
      
        for ( UInt iNode =0; iNode < tempNodeList.GetSize(); iNode++ ) {
          inhomDirichletNodes_.Push_back(tempNodeList[iNode]);
          if (dofsPerNode_ > 1) {
            dof = domain->GetCoordSystem()->GetVecComponent(dofs[i]);
            inhomDirichletDofs_.Push_back(dof);
          }
        } 
      }
    }
  }


  // ******************
  //   SetConstraints
  // ******************
  void BaseEQN::SetConstraints( const StdVector<UInt> &slaveNodeNrs,
                                const StdVector<UInt> &masterNodeNrs,
                                const StdVector<std::string> &dofs ) {

    ENTER_FCN( "BaseEQN::SetConstraints" );

    constraintSlaveNodes_ = slaveNodeNrs;
    constraintMasterNodes_ = masterNodeNrs;

    if ( dofsPerNode_ > 1 ) {
      constraintDofs_.Resize(dofs.GetSize());

      for ( UInt i = 0; i < dofs.GetSize(); i++ ) {
        constraintDofs_[i] =
          domain->GetCoordSystem()->GetVecComponent(dofs[i]);
      }
    }
  }


} // end of namespace
