// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <string> 

#include "OLAS/algsys/basesystem.hh"

#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/algsys/baseidbchandler.hh"
#include "OLAS/graph/basegraphmanager.hh"
#include "OLAS/graph/graphmanagersimple.hh"
#include "OLAS/graph/graphmanagerstdmat.hh"
#include "OLAS/graph/graphmanagersbmmat.hh"
#include "OLAS/solver/basesolver.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField {

  //#define DEBUG_BASESYSTEM

  // ***********************
  //   Default Constructor
  // ***********************
  BaseSystem::BaseSystem(PtrParamNode pn) 
  {
    graphManager_   = NULL;
    assemble_       = NULL;
    idbcHandler_    = NULL;

    solver_         = NULL;
    eigenSolver_    = NULL;

    bcOffsets_      = NULL;
    sizePerPDE_     = NULL;
    numLastFreeDof_ = NULL;

    xml_             = pn; 

    numPDEs_            = 0;
    numRegPDEs_         = 0;
    size_               = 0;
    nne_                = 0;
    blockSize_          = 1;
    numDirichletValues_ = 0;

    algSysType_     = NOALGSYSTYPE;
    
    olasInfo_ = info->Get("OLAS");
    systemInfo_ = olasInfo_->Get("system");

#ifdef MEMTRACE
    if ( memtrace == NULL ) {
      EXCEPTION("MEMTRACE macro is set, but memtrace stream pointer is "
               << "NULL!\n"
               << " Please check your settings in OLAS/main/Makefile.option");
    }
#endif

    // Default is to always use a system matrix
    matrixTypes_.insert( SYSTEM );
    
    // Set flag for insertion of penalty terms into matrix
    PtrParamNode setupNode;
    setupNode = xml_->Get("setup", ParamNode::INSERT );

    usingPenalty_ = true;
    std::string aux = "penalty";
    setupNode->GetValue("idbcHandling", aux, ParamNode::INSERT );
    usingPenalty_ = aux == "penalty" ? true : false;

    // Set flag for insertion of penalty terms into matrix
    if ( usingPenalty_ ) {
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


    delete solver_;
    solver_ = NULL;

    delete graphManager_;
    graphManager_ = NULL;

    delete assemble_;
    assemble_ = NULL;

    delete idbcHandler_;
    idbcHandler_ = NULL;

    // Delete arrays with PDE specific information
    delete [] ( numLastFreeDof_ );
    delete [] ( sizePerPDE_     );
    delete [] ( bcOffsets_      );

  }

  // ***************
  //   ObtainPDEId
  // ***************
  PdeIdType BaseSystem::ObtainPDEId( const std::string &pdeType ) {


   // (*debug) << "\n Obtain PDE Id for pde " << pdeType;

    // Check if PDE is already registered
    std::map<PdeIdType,std::string>::iterator it;
    for ( it = pdeNames_.begin(); it != pdeNames_.end(); it++ ) {
      if ( (*it).second == pdeType ) {
        EXCEPTION("A PDE with name '" << pdeType 
                 << "' was already registerd!");
      }
    }
    
    // Create Id
    PdeIdType id = pdeNames_.size();
    pdeNames_[id] = pdeType;

    // (*debug) << "\n --> Id is " << id << std::endl;

    return id;
  }


  // ***************
  //   RegisterPDE
  // ***************
  void BaseSystem::RegisterPDE( const PdeIdType pdeId,
                                UInt const numEqns,
                                UInt const numLastFreeDof ) {
    
     
   /* (*debug) << "\n RegisterPDE:"
             << "\n --> pdeId        = '" << pdeId << "'"
             << "\n --> numEqns        = " << numEqns
             << "\n --> numLastFreeDof = " << numLastFreeDof
             << std::endl;
   */

    // Remember number of equation numbers for each PDE
    sizePerPDE_[pdeId] = numEqns;

    // Remember equation numbers of last free dof each PDE
    numLastFreeDof_[pdeId] = numLastFreeDof;

    // Reordering must be handled differently for different
    // types of algebraic system
    switch ( algSysType_ ) {

    case STANDARD_SYSTEM:
      {
        // Obtain reordering type from OLAS-Params
        BaseOrdering::ReorderingType reorder = BaseOrdering::NOREORDERING;
        std::string reorderStr = "noReordering";
        xml_->Get( "matrix",ParamNode::INSERT)
        ->GetValue( "reordering", reorderStr, ParamNode::INSERT );
        reorder = BaseOrdering::reorderingType.Parse( reorderStr ); 

        // Register PDE with the graphmanager
        graphManager_->RegisterPDE( pdeId, numEqns, numLastFreeDof, reorder );
      }
      break;

    case SBM_SYSTEM:

      // Register PDE with the graphmanager
      graphManager_->RegisterPDE( pdeId, numEqns, numLastFreeDof,
                                  BaseOrdering::NOREORDERING );
      break;

    default:
      EXCEPTION("Derived class did not set algSysType_ to appropriate "
               << "value!");
    }
  }


  // **********************
  //   SetNumDirichletBCs
  // **********************
  void BaseSystem::SetNumDirichletBCs( const PdeIdType pdeID,
                                       const UInt numBCs ) {


    // check if pdeID is known
    if ( pdeID  > static_cast<Integer>(numPDEs_) ) {
      EXCEPTION("BaseSystem::SetnumDirichletBCs: A PDE with the ID '" 
               << pdeID << "' was not registered!");
    }

    // determine offsets
    if ( pdeID == 0 ) {
      bcOffsets_[0] = 0;
      bcOffsets_[1] = numBCs;
    } else {
      bcOffsets_[pdeID+1] = bcOffsets_[pdeID] + numBCs;
    }

    // store number of inhomogeneous Dirichlet boundary values
    numDirichletValues_ += numBCs;

    // check consistency
    if ( !usingPenalty_ ) {
      if ( numBCs != sizePerPDE_[pdeID] - numLastFreeDof_[pdeID] ) {
        EXCEPTION("BaseSystem::SetNumDirichletBCs:"
                 << " Inconsistency detected!\n numBCs != "
                 << "(numEQNs - numLastFreeDof) for PDE '"
                 << pdeNames_[pdeID] << "'\n"
                 << "\n numBCs ........... " << numBCs
                 << "\n sizePerPDE ....... " << sizePerPDE_[pdeID]
                 << "\n numLastFreeDof ... " << numLastFreeDof_[pdeID]);
      }
    }
  }


  // *******************
  //   GetFEMatrixType
  // *******************
  void BaseSystem::GetFEMatrixTypes( std::set<FEMatrixType> &matTypes ) const{
    matTypes = matrixTypes_;
  }


  // =======================================================================
  //   Methods for creating and assembling the matrix graph
  // =======================================================================


  // ******************
  //   GraphSetupInit
  // ******************
  void BaseSystem::GraphSetupInit( UInt numPDEs ) {

    
    // Store number of PDEs and resize bcOffset vector and
    // sizePerPDE-array
    numPDEs_ = numPDEs;
    NEWARRAY( bcOffsets_     , UInt, numPDEs_ + 1 );
    NEWARRAY( sizePerPDE_    , UInt, numPDEs_     );
    NEWARRAY( numLastFreeDof_, UInt, numPDEs_     );

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
    graphManager_->AssembleInit( identifierPDE1, identifierPDE2,
                                 assemblingTranspose );
  }


  // ****************
  //   AssembleDone
  // ****************
  void BaseSystem::AssembleDone( const PdeIdType identifierPDE1,
                                 const PdeIdType identifierPDE2,
                                 bool assemblingTranspose ) {
      graphManager_->AssembleDone( identifierPDE1, identifierPDE2,
                                   assemblingTranspose );
  }


  // *****************
  //   SetElementPos
  // *****************
  void BaseSystem::SetElementPos( const PdeIdType identifierPDE1,
                                  const StdVector<Integer>& eqnNrs1,
                                  const PdeIdType identifierPDE2,
                                  const StdVector<Integer>& eqnNrs2,
                                  bool setCounterPart ) {


    graphManager_->SetElementPos( identifierPDE1, eqnNrs1,
                                  identifierPDE2, eqnNrs2,
                                  setCounterPart );
  }


  // *****************
  //   GetReordering
  // *****************
  void BaseSystem::GetReordering( const PdeIdType identifier,
                                  StdVector<UInt>& newOrder ) {
    graphManager_->GetReordering(identifier, newOrder);
  }


  // *************************
  //   PrintRegistrationInfo
  // *************************
  void BaseSystem::PrintRegistrationInfo( std::ostream *log ) const {

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

  Integer BaseSystem::GetNumIter() {
    return olasInfo_->Get( "NumIter" )->As<Integer>();
  }

}// end of Namespace
