// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <string> 

#include "algsys/basesystem.hh"

#include "matvec/matvec.hh"
#include "graph/graph.hh"
#include "precond/precond.hh"
#include "solver/solver.hh"
#include "utils/utils.hh"
#include "algsys/olasparams.hh"
#include "algsys/baseentrymanipulator.hh"
#include "algsys/baseidbchandler.hh"
#include "graph/basegraphmanager.hh"


namespace OLAS {

  // ***********************
  //   Default Constructor
  // ***********************
  BaseSystem::BaseSystem() {

    ENTER_FCN( "BaseSystem::BaseSystem" );

    graphManager_   = NULL;
    assemble_       = NULL;
    idbcHandler_    = NULL;

    solver_         = NULL;
    eigenSolver_    = NULL;

    bcOffsets_      = NULL;
    sizePerPDE_     = NULL;
    numLastFreeDof_ = NULL;

    numPDEs_            = 0;
    numRegPDEs_         = 0;
    size_               = 0;
    nne_                = 0;
    blockSize_          = 0;
    numDirichletValues_ = 0;

    algSysType_     = NOALGSYSTYPE;

#ifdef DEBUG
    if ( debug == NULL ) {
      (*error) << "DEBUG macro is set, but debug stream pointer is NULL!\n"
               << " Please check your settings in OLAS/main/Makefile.option";
      Error( __FILE__, __LINE__ );
    }
#endif

#ifdef MEMTRACE
    if ( memtrace == NULL ) {
      (*error) << "MEMTRACE macro is set, but memtrace stream pointer is "
               << "NULL!\n"
               << " Please check your settings in OLAS/main/Makefile.option";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Default is to always use a system matrix
    matrixTypes_.insert( SYSTEM );

    // Set flag for insertion of penalty terms into matrix
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      assembleDirichletToSysMat_ = true;
    }
    else {
      assembleDirichletToSysMat_ = false;
    }
  }


  // **************
  //   Destructor
  // **************
  BaseSystem::~BaseSystem() {

    ENTER_FCN( "BaseSystem::~BaseSystem" );

    delete solver_;
    solver_ = NULL;

    delete graphManager_;
    graphManager_ = NULL;

    delete assemble_;
    assemble_ = NULL;

    delete idbcHandler_;
    idbcHandler_ = NULL;

    // Delete arrays with PDE specific information
    DeleteArray( numLastFreeDof_ );
    DeleteArray( sizePerPDE_     );
    DeleteArray( bcOffsets_      );

  }


  // ***************
  //   ObtainPDEId
  // ***************
  PdeIdType BaseSystem::ObtainPDEId( const std::string &pdeType ) {

    ENTER_FCN( "BaseSystem::ObtainPDEId" );

#ifdef DEBUG_BASESYSTEM
    (*debug) << "\n Obtain PDE Id for pde " << pdeType;
#endif

    // Check if PDE is already registered
    std::map<PdeIdType,std::string>::iterator it;
    for ( it = pdeNames_.begin(); it != pdeNames_.end(); it++ ) {
      if ( (*it).second == pdeType ) {
        (*error) << "A PDE with name '" << pdeType 
                 << "' was already registerd!";
        Error( __FILE__, __LINE__ );
      }
    }
    
    // Create Id
    PdeIdType id = pdeNames_.size() + 1;
    pdeNames_[id] = pdeType;

#ifdef DEBUG_BASESYSTEM
    (*debug) << "\n --> Id is " << id << std::endl;
#endif

    return id;
  }


  // ***************
  //   RegisterPDE
  // ***************
  void BaseSystem::RegisterPDE( const PdeIdType pdeId,
                                UInt const numEqns,
                                UInt const numLastFreeDof ) {
    
    ENTER_FCN( "BaseSystem::RegisterPDE" );
    
#ifdef DEBUG_BASESYSTEM
    (*debug) << "\n RegisterPDE:"
             << "\n --> pdeId        = '" << pdeId << "'"
             << "\n --> numEqns        = " << numEqns
             << "\n --> numLastFreeDof = " << numLastFreeDof
             << std::endl;
#endif

    // Remember number of equation numbers for each PDE
    sizePerPDE_[pdeId] = numEqns;

    // Remember equation numbers of last free dof each PDE
    numLastFreeDof_[pdeId] = numLastFreeDof;

    // Reordering must be handled differently for different
    // types of algebraic system
    switch ( algSysType_ ) {

    case STANDARD_SYSTEM:

      // Obtain reordering type from OLAS-Params
      ReorderingType reorder;
      myParams_.GetEnumValue( "GRAPH_reordering", reorder );

      // Register PDE with the graphmanager
      graphManager_->RegisterPDE( pdeId, numEqns, numLastFreeDof, reorder );

      break;

    case SBM_SYSTEM:

      // Register PDE with the graphmanager
      graphManager_->RegisterPDE( pdeId, numEqns, numLastFreeDof,
                                  NOREORDERING );
      break;

    default:
      (*error) << "Derived class did not set algSysType_ to appropriate "
               << "value!";
      Error( __FILE__, __LINE__ );
    }
  }


  // **********************
  //   SetNumDirichletBCs
  // **********************
  void BaseSystem::SetNumDirichletBCs( const PdeIdType pdeID,
                                       const UInt numBCs ) {

    ENTER_FCN( "BaseSystem::SetNumDirichletBCs" );

    // check if pdeID is known
    if ( pdeID  > numPDEs_ ) {
      (*error) << "BaseSystem::SetnumDirichletBCs: A PDE with the ID '" 
               << pdeID << "' was not registered!";
      Error( __FILE__, __LINE__ );
    }

    // determine offsets
    if ( pdeID == 1 ) {
      bcOffsets_[1] = 0;
      bcOffsets_[2] = numBCs;
    } else {
      bcOffsets_[pdeID+1] = bcOffsets_[pdeID] + numBCs;
    }

    // store number of inhomogeneous Dirichlet boundary values
    numDirichletValues_ += numBCs;

    // check consistency
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == false ) {
      if ( numBCs != sizePerPDE_[pdeID] - numLastFreeDof_[pdeID] ) {
        (*error) << "BaseSystem::SetNumDirichletBCs:"
                 << " Inconsistency detected!\n numBCs != "
                 << "(numEQNs - numLastFreeDof) for PDE '"
                 << pdeNames_[pdeID] << "'\n"
                 << "\n numBCs ........... " << numBCs
                 << "\n sizePerPDE ....... " << sizePerPDE_[pdeID]
                 << "\n numLastFreeDof ... " << numLastFreeDof_[pdeID];
        Error( __FILE__, __LINE__ );
      }
    }
  }


  // *******************
  //   GetFEMatrixType
  // *******************
  void BaseSystem::GetFEMatrixTypes( std::set<FEMatrixType> &matTypes ) const{
    ENTER_FCN( "BaseSystem::GetFEMatrixTypes" );
    matTypes = matrixTypes_;
  }


  // =======================================================================
  //   Methods for creating and assembling the matrix graph
  // =======================================================================


  // ******************
  //   GraphSetupInit
  // ******************
  void BaseSystem::GraphSetupInit( UInt numPDEs ) {

    ENTER_FCN( "BaseSystem::GraphSetupInit" );
    
    // Store number of PDEs and resize bcOffset vector and
    // sizePerPDE-array
    numPDEs_ = numPDEs;
    NewArray( bcOffsets_     , UInt, numPDEs_ + 1 );
    NewArray( sizePerPDE_    , UInt, numPDEs_     );
    NewArray( numLastFreeDof_, UInt, numPDEs_     );

    // Depending on the type of algebraic system and on the number of PDEs,
    // we generate a fitting graph manager object
    if ( algSysType_ == STANDARD_SYSTEM ) {
      if( numPDEs == 1 ) {
        graphManager_ = new GraphManagerSimple();
      }
      else {
        graphManager_ = new GraphManagerStdMat();
      }
    }
    else {
      graphManager_ = new GraphManagerSBMMat();
    }

    // start setup of GraphManager
    graphManager_->SetupInit( numPDEs_ );
  }


  // ******************
  //   GraphSetupDone
  // ******************
  void BaseSystem::GraphSetupDone() {
    ENTER_FCN( "BaseSystem::GraphSetupDone" );

    // Print information about registered PDEs
    PrintRegistrationInfo( cla );
    
    // Finalise graph manager setup
    graphManager_->SetupDone();
  }


  // ****************
  //   AssembleInit
  // ****************
  void BaseSystem::AssembleInit( const PdeIdType identifierPDE1,
				 const PdeIdType identifierPDE2,
                                 bool assemblingTranspose ) {
    ENTER_FCN( "BaseSystem::AssembleInit" );
    graphManager_->AssembleInit( identifierPDE1, identifierPDE2,
                                 assemblingTranspose );
  }


  // ****************
  //   AssembleDone
  // ****************
  void BaseSystem::AssembleDone( const PdeIdType identifierPDE1,
				 const PdeIdType identifierPDE2,
                                 bool assemblingTranspose ) {
      ENTER_FCN( "BaseSystem::AssembleDone" );
      graphManager_->AssembleDone( identifierPDE1, identifierPDE2,
                                   assemblingTranspose );
  }


  // *****************
  //   SetElementPos
  // *****************
  void BaseSystem::SetElementPos( const PdeIdType identifierPDE1,
				  Integer *connect1,
				  Integer elemSize1,
				  const PdeIdType identifierPDE2,
				  Integer *connect2,
				  Integer elemSize2,
                                  bool setCounterPart ) {

    ENTER_FCN( "BaseSystem::SetElementPos" );

    // Adjust 0-based index to 1-based
    Integer* connectOneBased1 = connect1;
    Integer* connectOneBased2 = connect2;

    connectOneBased1--;
    if ( connect2 != NULL ) {
      connectOneBased2--;
    }

    graphManager_->SetElementPos( identifierPDE1, connectOneBased1, 
				  elemSize1, identifierPDE2,
				  connectOneBased2, elemSize2,
                                  setCounterPart );
  }


  // *****************
  //   GetReordering
  // *****************
  Integer* BaseSystem::GetReordering( const PdeIdType identifier ) {
    ENTER_FCN( "BaseSystem::GetReordering" );
    Integer *newOrder = graphManager_->GetReordering(identifier);
    if ( newOrder != NULL ) {
      newOrder++;
    }
    return newOrder;
  }


  // *************************
  //   PrintRegistrationInfo
  // *************************
  void BaseSystem::PrintRegistrationInfo( std::ostream *log ) const {
    ENTER_FCN( "BaseSystem::PrintRegistrationInfo" );

    // Print name, id, number of unknowns and number of 
    // boundary conditions

    (*log) << " ================================\n"
	   << "   PDE Registration Information  \n"
	   << " ================================\n\n";

    (*log) << "     PDE-name    |  Id  |  Unknowns  | numLastFreeDof \n"
	   << " -----------------------------------------------------\n";

    std::map<PdeIdType,std::string>::const_iterator it; 

    for ( it = pdeNames_.begin(); it != pdeNames_.end(); it++ ) {
      (*log) << " "
             << std::setw(15) << (*it).second             << " | "
	     << std::setw(4)  << (*it).first              << " | "
	     << std::setw(10) << sizePerPDE_[(*it).first] << " | "
	     << std::setw(8)  << numLastFreeDof_[(*it).first] 
	     << std::endl;
    }
    (*log) << std::endl;

  }


}// end of Namespace
