// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "algsys/olascomm.hh"
#include "algsys/sbmsystem.hh"

#include "utils/utils.hh"
#include "graph/graphmanagersbmmat.hh"
#include "matvec/matvec.hh"
#include "solver/solver.hh"
#include "precond/precond.hh"

#include "algsys/generateentrymanipulator.hh"
#include "algsys/baseentrymanipulator.hh"
#include "algsys/generateidbchandler.hh"
#include "algsys/baseidbchandler.hh"


namespace OLAS {


  // ***********************
  //   Default Constructor
  // ***********************
  SBM_System::SBM_System(ParamNode* pn) : BaseSystem(pn) {


    precond_             = NULL;
    rhs_                 = NULL;
    sol_                 = NULL;
    feSubMatrices_       = NULL;

    // We are an SBM_System obviously
    algSysType_          = SBM_SYSTEM;

    // Initialize pointer to finite-element matrices
    NewArray( sysMat_, SBM_Matrix*, MAX_NUM_FE_MATRICES );
    for ( UInt i = 1; i <= MAX_NUM_FE_MATRICES; i++ ) {
      sysMat_[i] = NULL;
    }

    // Assuming non-symmetric matrices by default
    sbmSymm_ = false;

    // Generate array for storing FE-matrix pattern information
    feSubMatrices_ = new(std::nothrow) std::set<SubMatrixID,SortSubMatrixID>
      [MAX_NUM_FE_MATRICES];
    if ( feSubMatrices_ == NULL ) {
      (*error) << "Failed to allocate feSubMatrices_ array";
      Error( __FILE__, __LINE__ );
    }
    feSubMatrices_--;
  }


  // **************
  //   Destructor
  // **************
  SBM_System::~SBM_System() {


    // Delete our own objects
    delete precond_;
    delete sol_;
    delete rhs_;

    // Delete finite-element matrices and sub-matrix info
    for ( UInt i = 1; i <= MAX_NUM_FE_MATRICES; i++ ) {
      delete sysMat_[i];
    }
    DeleteArray( sysMat_ );
    DeleteArray( feSubMatrices_ );
  }


  // ****************
  //   SetupPrecond
  // ****************
  void SBM_System::SetupPrecond() {
    precond_->Setup( *sysMat_[SYSTEM] );
  }


  // ****************
  //   SetBlockSize
  // ****************
  void SBM_System::SetBlockSize( const PdeIdType identifier, const UInt bs ) {


    // Check consistency
    if ( bs > 1 ) {
      (*error) << "SBM_System::SetBlockSize: "
               << "We do not support block-size = " << bs
               << ", but only scalar entries!";
      Error( __FILE__, __LINE__ );
    }
    else {
      blockSize_ = 1;
    }
  }


  Integer SBM_System::GetSolutionVal( Double* &ptSol,
                                      const PdeIdType identifierPDE ) {

    // iteratve over all rhs-vectors and copy entries
    Integer actPos = 1;
    Vector<Double> & bufVec = dynamic_cast<Vector<Double>&>(*solBuffer_);
    bufVec.Resize( size_ );
    bufVec.Init();
    for( Integer i = 1; i <= numPDEs_; i++ ) {
      Vector<Double> & actVec = dynamic_cast<Vector<Double>&>
      ( *(sol_->GetPointer( i )) );

      for( Integer iEqn = 1; iEqn <= numLastFreeDof_[i]; iEqn++, actPos++ ) {
        bufVec[actPos] = actVec[iEqn];
      }
    }
    Double * ptr  = (dynamic_cast<Vector<Double>& >
    (*solBuffer_)).GetPointer();
    ptr++;
    ptSol = ptr;
    Warning( "The ordering of results in the result file will not be correct!");
    return size_;
  }

  Integer SBM_System::GetSolutionVal( Complex* &ptSol,
                                      const PdeIdType identifierPDE ) {

    // iteratve over all rhs-vectors and copy entries
     Integer actPos = 1;
     Vector<Complex> & bufVec = dynamic_cast<Vector<Complex>&>(*solBuffer_);
     bufVec.Resize( size_ );
     bufVec.Init();
     for( Integer i = 1; i <= numPDEs_; i++ ) {
       Vector<Complex> & actVec = dynamic_cast<Vector<Complex>&>
       ( *(sol_->GetPointer( i )) );

       for( Integer iEqn = 1; iEqn <= numLastFreeDof_[i]; iEqn++, actPos++ ) {
         bufVec[actPos] = actVec[iEqn];
       }
     }
     Complex * ptr  = (dynamic_cast<Vector<Complex>& >
     (*solBuffer_)).GetPointer();
     ptr++;
     ptSol = ptr;
     Warning( "The ordering of results in the result file will not be correct!");
     return size_;
  }

  Integer SBM_System::GetRHSVal( Double* &ptRhs,
                                 const PdeIdType identifierPDE ) {
    Double * ptr  = (dynamic_cast<Vector<Double>& >
      (*rhsBuffer_)).GetPointer();
     ptr++;
     ptRhs = ptr;
     Warning( "The ordering of results in the result file will not be correct!");
     return size_;

  }

  Integer SBM_System::GetRHSVal( Complex* &ptRhs,
                                 const PdeIdType identifierPDE ) {
    Complex * ptr  = (dynamic_cast<Vector<Complex>& >
    (*rhsBuffer_)).GetPointer();
    ptr++;
    ptRhs = ptr;
    Warning( "The ordering of results in the result file will not be correct!");
    return size_;
  }


  // *******************
  //   SetFEMatrixType
  // *******************
  void SBM_System::SetFEMatrixType( const FEMatrixType matType,
                                    const PdeIdType pdeID1,
                                    const PdeIdType pdeID2 ) {


    // Determine row and column index of sub-matrix
    UInt rowInd = pdeID1;
    UInt colInd = pdeID2 == NO_PDE_ID ? pdeID1 : pdeID2;

    // Do nothing for nothing :)
    if( matType != NOTYPE ) {

      // Insert FE-matrix type
      matrixTypes_.insert( matType );

      // Construct sub-matrix identifier
      SubMatrixID sID;
      sID.rowInd = rowInd;
      sID.colInd = colInd;

      // Insert sub-matrix identifier into
      // corresponding FE-Matrix set
      feSubMatrices_[matType].insert( sID );

      // -----------------------------------------------------------
      // NOTE: We temporarily insert also the transposed sub-matrix
      if ( rowInd != colInd ) {
        sID.colInd = rowInd;
        sID.rowInd = colInd;

        // Insert sub-matrix identifier into
        // corresponding FE-Matrix set
        feSubMatrices_[matType].insert( sID );
      }
      // -----------------------------------------------------------

    }
  }


  // ****************
  //   CreateLinSys
  // ****************
  void SBM_System::CreateLinSys() {

    std::set<FEMatrixType>::iterator fIt;
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;

    // Obtain symmetry flag
    sbmSymm_ = myParams_.GetBoolValue( "SBM_Symmetry" );

    // Determine the set of sub-matrices that we need for the system matrix
    // There are two different cases:
    //
    // 1) We have other FE matrices, then we generate the set by creating
    //    the union of all other sets
    //
    // 2) There are no other sets, so we generate a sub-matrix for each
    //    sub-graph that the graph manager stores
    if ( matrixTypes_.size() > 1 ) {
      for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
        if ( *fIt != SYSTEM ) {
          feSubMatrices_[SYSTEM].insert( feSubMatrices_[*fIt].begin(),
                                         feSubMatrices_[*fIt].end() );
        }
      }
    }
    else {
      GraphManagerSBMMat *sbmManager =
        dynamic_cast<GraphManagerSBMMat*>(graphManager_);
      SubMatrixID sID;
      for ( UInt i = 1; i <= numPDEs_; i++ ) {
        for ( UInt j = 1; j <= numPDEs_; j++ ) {
          if ( sbmManager->SubGraphExists( i, j ) == true ) {
            sID.rowInd = i;
            sID.colInd = j;
            feSubMatrices_[SYSTEM].insert( sID );
          }
        }
      }
    }


    // determine overall number of unknowns
    size_ = 0;
    for( Integer i = 1; i <= numPDEs_; i++ ) {
      size_ += numLastFreeDof_[i];
    }

    // ------------------------------
    //  Generation of matrix objects
    // ------------------------------

    // Log what we will do
    PrintFeMatrixInfo( cla );

    // Obtain some info from parameter file
    MatrixEntryType   entryType;
    myParams_.GetEnumValue( "MatrixEntryType"  , entryType   );

    for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
      sysMat_[*fIt] = GenerateSBM_Matrix( *fIt, entryType );
    }


    // --------------------------------------------
    //  Treatment of Dirichlet Boundary Conditions
    // --------------------------------------------

    // Check the case we have (elimination vs. penalty) and generate
    // an appropriate object for handling inhomogeneous Dirichlet
    // boundary conditions
    if ( myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      idbcHandler_ = GenerateIDBC_HandlerObject( numDirichletValues_,
                                                 numPDEs_, bcOffsets_,
                                                 entryType );
    }
    else {
      idbcHandler_ = GenerateIDBC_HandlerObject( matrixTypes_,
                                                 graphManager_,
                                                 numPDEs_,
                                                 numDirichletValues_,
                                                 entryType,
                                                 true );
    }


    // ------------------------------
    //  Generation of vector objects
    // ------------------------------

    // Generate empty SBM vectors
    rhs_ = dynamic_cast<SBM_Vector*>
      ( GenerateVectorObject( *(sysMat_[SYSTEM]) ) );

    sol_ = dynamic_cast<SBM_Vector*>
      ( GenerateVectorObject( *(sysMat_[SYSTEM]) ) );

    if ( rhs_ == NULL || sol_ == NULL ) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

    // For the moment we insert a sub-vector for each position.
    // In the case of the right-hand side we might actually be
    // more economic. How do we get the information which sub-vectors
    // are really needed, however?
    StdMatrix *stdMat = NULL;
    BaseVector *bVec = NULL;
    SparseVector *sVec = NULL;
    for ( UInt k = 1; k <= numPDEs_; k++ ) {

      // Get diag matrix for vector generation
      stdMat = sysMat_[SYSTEM]->GetPointer( k, k );

      // Insert sub-vector into solution
      bVec = GenerateVectorObject( *stdMat );
      sVec = dynamic_cast<SparseVector*>( bVec );
      sol_->SetSubVector( sVec, k );

      // Insert sub-vector into right-hand side
      bVec = GenerateVectorObject( *stdMat );
      sVec = dynamic_cast<SparseVector*>( bVec );
      rhs_->SetSubVector( sVec, k );
    }

    // ----------------------------------
    //  Generation of auxilliary objects
    // ----------------------------------

    // Create the assemble object
    assemble_ = GenerateEntryManipulatorObject( entryType, 1 );

    // Generate sparsevectors

    MatrixStorageType sType = stdMat->GetStorageType();
    solBuffer_ = dynamic_cast<SparseVector*>
        ( GenerateSparseVectorObject( sType, entryType, 1, size_ ) );
    rhsBuffer_ = dynamic_cast<SparseVector*>
        ( GenerateSparseVectorObject( sType, entryType, 1, size_ ) );

    // -----------------
    //  Memory clean-up
    // -----------------

    // At this point, hopefully, the graph object is no longer
    // required by anyone, so release pointer and delete manager
    // to free memory
    delete graphManager_;
    graphManager_ = NULL;
  }


  // *********************
  //   PrintFeMatrixInfo
  // *********************
  void SBM_System::PrintFeMatrixInfo( std::ostream *os ) {


    std::set<FEMatrixType>::iterator fIt;
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;

    // Print header
    UInt tw = 60;
    (*os) << "\n " << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
           << std::setfill( ' ' );

    // Print information about matrix types to .las file
    if ( sbmSymm_ == false ) {
      (*os) << " SBM_System will generate the following FE matrices:"
            << "\n (sub-matrix set = 'x', not set = '.' )\n"
            << std::endl;
    }
    else {
      (*os) << " SBM_System will generate the following FE matrices:"
            << "\n (sub-matrix set = 'x', not set = '.', dropped due to "
            << "symmetry = '*' )\n" << std::endl;
    }

    for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {

      // Log FE matrix name
      (*os) << " " << Enum2String(*fIt) << std::endl;

      // Now plot sub-matrix pattern
      SubMatrixID sID;
      for ( UInt i = 1; i <= numPDEs_; i++ ) {
        sID.rowInd = i;
        for ( UInt j = 1; j <= numPDEs_; j++ ) {
          sID.colInd = j;
          sIt = feSubMatrices_[*fIt].find(sID);
          if ( sIt != feSubMatrices_[*fIt].end() ) {
            if ( sbmSymm_ == true && sID.rowInd > sID.colInd ) {
              (*os) << " *";
            }
            else {
              (*os) << " x";
            }
          }
          else {
            (*os) << " .";
          }
        }
        (*os) << std::endl;
      }
    }

    // Print footer
    (*os) << std::setw(tw) << std::setfill( '=' ) << '=' << "\n"
          << std::setfill( ' ' );
  }


  // **********************
  //   GenerateSBM_Matrix
  // **********************
  SBM_Matrix* SBM_System::GenerateSBM_Matrix( FEMatrixType matType,
                                              MatrixEntryType entryType ) {


    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;
    retMat = New SBM_Matrix( numPDEs_, numPDEs_, sbmSymm_ );
    if ( retMat == NULL ) {
      (*error) << "SBM_System::GenerateSBM_Matrix: "
               << "This is the end my friend!\n"
               << "Generation of empty SBM_Matrix failed!";
      Error( __FILE__, __LINE__ );
    }

    // STEP 2: Populate with sub-matrices
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;
    BaseGraph *graph;
    for ( sIt = feSubMatrices_[matType].begin();
          sIt != feSubMatrices_[matType].end(); sIt++ ) {

      // Determine associated PDE identifiers
      PdeIdType pde1 = (*sIt).rowInd;
      PdeIdType pde2 = (*sIt).colInd;

      // Determine number of matrix rows and columns
      // graph = graphManager_->GetGraph( pde1, pde1 );
      // UInt nrows = graph->GetSize();
      // graph = graphManager_->GetGraph( pde2, pde2 );
      // UInt ncols = graph->GetSize();

      // Determine number of matrix rows and columns
      UInt nrows = numLastFreeDof_[pde1];
      UInt ncols = numLastFreeDof_[pde2];

      // Check for necessity of generation
      if ( pde1 <= pde2 || sbmSymm_ == false ) {

        // Trigger generation of sub-matrix
        graph = graphManager_->GetGraph( pde1, pde2 );
        if ( pde1 == pde2 && sbmSymm_ == true ) {
          retMat->SetSubMatrix ( pde1, pde2, entryType, SPARSE_SYM, 1,
                                 nrows, ncols, graph->GetNNE() );
        }
        else {
          retMat->SetSubMatrix ( pde1, pde2, entryType, SPARSE_NONSYM, 1,
                                 nrows, ncols, graph->GetNNE() );
        }

        // Set sparsity pattern of sub-matrix
        (*retMat)( pde1, pde2 ).SetSparsityPattern( *graph );
      }
    }

    return retMat;
  }


  // ********************
  //   SetElementMatrix
  // ********************
  void SBM_System::SetElementMatrix( FEMatrixType matrixID, Double *elemMat,
                                     PdeIdType idPDE1, Integer *eqnNrs1,
                                     Integer numEqn1, PdeIdType idPDE2,
                                     Integer *eqnNrs2, Integer numEqn2,
                                     bool setCounterPart ) {

    // Set flag for setting the symmetric counter-part of
    // the element matrix
    if ( setCounterPart == true && idPDE1 == idPDE2 ) {
      setCounterPart = false;
    }

    // Do not try to set lower-part if symmetric
    bool reverse = false;
    if ( sbmSymm_ == true ) {

      // Matrix in upper part, then do not set counter part
      if ( idPDE1 <= idPDE2 ) {
        setCounterPart = false;
      }

      // Matrix in lower part and we shall not set counter part
      // then we have a major problem, since SBM matrix should
      // not be symmetric
      else if ( setCounterPart == false ) {
        (*error) << "SBM_System::SetElementMatrix:\n"
                 << " idPDE1 ........... " << idPDE1 << '\n'
                 << " idPDE2 ........... " << idPDE2 << '\n'
                 << " setCounterPart ... " << std::boolalpha << setCounterPart
                 << " sbmSymm .......... true\n"
                 << " leads to loss of element matrix!";
        Error( __FILE__, __LINE__ );
      }

      // Matrix in lower part and counter part should be set,
      // so we change to the counterPart
      else {
        setCounterPart = false;
        reverse = true;
      }
    }

    // Extract associated sub-matrix
    StdMatrix *subMat1 = &((*sysMat_[matrixID])( idPDE1, idPDE2 ));
    StdMatrix *subMat2 = &((*sysMat_[matrixID])( idPDE2, idPDE1 ));

    // Delegate remaining work to assemble
    if ( reverse == false ) {
      assemble_->SetElementMatrix( matrixID, idPDE1, idPDE2,
                                   subMat1, subMat2,
                                   idbcHandler_, elemMat,
                                   eqnNrs1, (UInt)numEqn1,
                                   eqnNrs2, (UInt)numEqn2,
                                   numLastFreeDof_[idPDE1],
                                   numLastFreeDof_[idPDE2],
                                   setCounterPart);
    }
    else {
      assemble_->SetCounterPartOnly( matrixID, idPDE1, idPDE2,
                                     subMat1, subMat2,
                                     idbcHandler_, elemMat,
                                     eqnNrs1, (UInt)numEqn1,
                                     eqnNrs2, (UInt)numEqn2,
                                     numLastFreeDof_[idPDE1],
                                     numLastFreeDof_[idPDE2] );
    }
  }


  // **************
  //   InitMatrix
  // **************
  void SBM_System::InitMatrix( FEMatrixType matrixID,
                               const PdeIdType pdeID ) {


    std::set<FEMatrixType>::iterator fIt;
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;

    // If a matrix is specified only initialise that one
    if ( matrixID != NOTYPE ) {
      sysMat_[matrixID]->Init();
      idbcHandler_->InitMatrix(matrixID);
    }

    // otherwise initialise every Finite-Element matrix
    else {
      for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
        sysMat_[*fIt]->Init();
        idbcHandler_->InitMatrix(*fIt);
      }
    }

    // Do we need to re-insert penalty values into matrix?
    if ( ( matrixID == NOTYPE || matrixID == SYSTEM ) &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      assembleDirichletToSysMat_ = true;
    }
  }


  // ***********************
  //   InitRHS (Version 1)
  // ***********************
  void SBM_System::InitRHS( PdeIdType pdeID ) {


    if ( pdeID == NO_PDE_ID ) {
      rhs_->Init();
    }
    else {
      rhs_->GetPointer( pdeID )->Init();
    }
  }


  // ***********
  //   InitSol
  // ***********
  void SBM_System::InitSol( PdeIdType pdeID ) {


    if ( pdeID == NO_PDE_ID ) {
      sol_->Init();
    }
    else {
      sol_->GetPointer( pdeID )->Init();
    }
  }


  // ****************************
  //   ConstructEffectiveMatrix
  // ****************************
  void SBM_System::ConstructEffectiveMatrix( const factorMap &matFactors ) {

    factorMap::const_iterator it;
    SBM_Matrix *sys = sysMat_[SYSTEM];

    // It's okay, if there are no factors, if there is only a system
    // matrix and no other ones
    if ( matFactors.empty() == true ) {
      if ( matrixTypes_.size() == 1 && sys != NULL ) {
        //(*debug) << " ConstructEffectiveMatrix: Nothing to be done, since"
        //         << " there is only one FE = system matrix\n\n";

        // Also assemble the effective auxilliary system matrix for moving
        // IDBCs to the right-hand side
        idbcHandler_->BuiltSystemMatrix( matFactors );

        // Now we are done
        return;
      }
      else {
        (*warning) << "SBM_System::ConstructEffectiveMatrix: "
                   << "Map with factors is empty, but there are "
                   << matrixTypes_.size() << " FE matrices in the game!";
        Warning( __FILE__, __LINE__ );
      }
    }

    /*
    (*debug) << "ConstructEffectiveMatrix with factors: " << std::endl;
    for ( it = matFactors.begin(); it != matFactors.end(); it++ ) {
      (*debug) << Enum2String((*it).first) << ":"
               << (*it).second << std::endl;
    }
*/

    for ( it = matFactors.begin(); it != matFactors.end(); it++ ) {
      if ( sysMat_[(*it).first] != NULL  && (*it).second != 0.0 ) {
        sys->Add( (*it).second, *sysMat_[(*it).first] );
      }
    }

    // Also assemble the effective auxilliary system matrix for moving
    // IDBCs to the right-hand side
    idbcHandler_->BuiltSystemMatrix( matFactors );
  }


  // *****************
  //   SetElementRHS
  // *****************
  void SBM_System::SetElementRHS( Double *elemRHS,
                                  const PdeIdType idPDE,
                                  Integer *connect,
                                  UInt length ) {


    // Delegate work to EntryManipulator
    assemble_->SetElementRHS( rhs_->GetPointer(idPDE), elemRHS, connect,
                              length, numLastFreeDof_[idPDE] );
  }


  // *************
  //   UpdateRHS
  // *************
  void SBM_System::UpdateRHS( FEMatrixType matrix_id, Double *fup ) {
    (*error) << "Method not yet implemented!";
    Error( __FILE__, __LINE__ );
  }


  // *************************
  //   SetDirichlet (Double)
  // *************************
  void SBM_System::SetDirichlet(const PdeIdType pdeID,
                                Integer eqnNum, const Double &val ) {


    // Perform some consistency checks
#ifdef DEBUG_SBMSYSTEM

    // Some situation dependend modifiers
    bool usingPenalty = myParams_.GetBoolValue( "UsingPenaltyFormulation" );
    UInt maxVal = usingPenalty == true ? sizePerPDE_[pdeID] :
      numLastFreeDof_[pdeID];
    UInt minVal = usingPenalty == true ? 1 : numLastFreeDof_[pdeID] + 1;

    if ( eqnNum > maxVal || eqnNum < minVal  ) {
      (*error) << "SBMSystem::SetDirichlet: Inconsistency detected:"
               << "\n pdeID           = " << pdeID
               << "\n eqnNum          = " << eqnNum
               << "\n val             = " << val
               << "\n numBC           = " << numDirichletValues_
               << "\n size            = " << sizePerPDE_[pdeID]
               << "\n numLastFreeDof_ = " << numLastFreeDof_[pdeID]
               << "\n minVal          = " << minVal
               << "\n maxVal          = " << maxVal
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( pdeID, eqnNum,  val );
  }


  // **************************
  //   SetDirichlet (Complex)
  // **************************
  void SBM_System::SetDirichlet( const PdeIdType pdeID,
                                 Integer eqnNum, const Complex &val ) {


    // Perform some consistency checks
#ifdef DEBUG_SBMSYSTEM

    // Some situation dependend modifiers
    bool usingPenalty = myParams_.GetBoolValue( "UsingPenaltyFormulation" );
    UInt maxVal = sizePerPDE_[pdeID];
    UInt minVal = usingPenalty == true ? 1 : numLastFreeDof_[pdeID] + 1;

    if ( eqnNum > maxVal || eqnNum < minVal ) {
      (*error) << "SBMSystem::SetDirichlet: Inconsistency detected:"
               << "\n pdeID           = " << pdeID
               << "\n eqnNum          = " << eqnNum
               << "\n val             = " << val
               << "\n numBC           = " << numDirichletValues_
               << "\n sizePerPDE      = " << sizePerPDE_[pdeID]
               << "\n numLastFreeDof_ = " << numLastFreeDof_[pdeID]
               << "\n minVal          = " << minVal
               << "\n maxVal          = " << maxVal
               << "\n SystemName is '"
               << myParams_.GetStringValue( "SystemName" ) << "'";
      Error( __FILE__, __LINE__ );
    }
#endif

    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( pdeID, eqnNum, val );
  }


  // *******************
  //   BuildInDirichlet
  // *******************
  void SBM_System::BuildInDirichlet() {


    // If necessary modify matrix diagonal for penalty approach
    if ( assembleDirichletToSysMat_ == true ) {
      idbcHandler_->AdaptSystemMatrix( *(sysMat_[SYSTEM]) );
      assembleDirichletToSysMat_ = false;
    }
  }


  // *********
  //   Solve
  // *********
  void SBM_System::Solve(const std::string& comment) {


    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the righ-hand side is
    // "contaminated" with penalty terms
    if ( numDirichletValues_ > 0 &&
      myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
      myParams_.SetValue( "RHSwithPenalty", false );
    }

#ifdef PROFILING
    Double t1 = Profiler::GetRealTime();
#endif

    // Iterative solvers require an initial guess and in the penalty case
    // we should insert the Dirichlet values into it
    if ( dynamic_cast<BaseIterativeSolver*>(solver_) != NULL &&
         myParams_.GetBoolValue( "UsingPenaltyFormulation" ) == true ) {
    idbcHandler_->SetDofsToIDBC( sol_ );
    (*cla) << " Inserted Dirichlet values into initial guess"
           << std::endl;
    }

    // Perform a simple sanity check
    // if ( sysmat_[SYSTEM] == NULL || precond_ == NULL || rhs_ == NULL ||
    //      sol_ == NULL ) {
    //   Error( "Detected NULL pointer where there should be none!", __FILE__,
    //          __LINE__ );
    // }

    // Assume that everything will go well
    myReport_.SetValue( "solutionIsOkay", true );

    // Now modifiy the right-hand side vector
    idbcHandler_->AddIDBCToRHS( rhs_ );

    // Remove the export linear system stuff, it has changed in standardsys.cc an as below
    // the solve part is commentet out, I see no reason to export linsys also here, it would
    // require a generalization anyway. Fabian 16.11.07
    // check if we do export stuff
     ParamNode* els = xml  != NULL && xml->Has("exportLinSys") ? xml->Get("exportLinSys") : NULL;
     std::string file;
     std::string base;

     // need it common even when exclusive solution
     if(els) {
       std::ostringstream os;
       os << els->Get("baseName")->AsString();
       if(comment != "") os << "_" << comment;
       base = os.str();
     }

     // check if we do not only want the solution
     if(els && els->Get("solution")->AsString() != "exclusive")
     {
       // two formats. The harwell-boing format includes the rhs!
       if(els->Get("format")->AsString()  == "harwell-boeing")
       {
         Error( "Harwell-Boeing Format not implemented for SBM-case",
                __FILE__, __LINE__ );
       }
       else // classical (default) matrix-market
       {
         sysMat_[SYSTEM]->Export((base+".mtx").c_str() );

         if(els->Has("damping", true) && sysMat_[DAMPING] != NULL)
           sysMat_[DAMPING]->Export((base+"_damping.mtx").c_str() );

         if(els->Has("auxiliary", true) && sysMat_[AUXILIARY] != NULL)
           sysMat_[AUXILIARY]->Export((base+"_aux.mtx").c_str() );

         // rhs is only in harwell-boing included
         rhs_->Export((base+".vec").c_str() );
       }
     }
     if(els && els->Has("initialGuess", true))
           sol_->Export((base+"_intial_guess.vec").c_str());

    // Trigger solution
    solver_->Solve( *sysMat_[SYSTEM], *precond_, *rhs_, *sol_ );

    // Now de-modifiy the right-hand side vector
    idbcHandler_->RemoveIDBCFromRHS( rhs_ );

    // Check that solution went fine, if not issue a warning
    if ( myReport_.GetBoolValue( "solutionIsOkay" ) == false ) {
      (*warning) << "Solver reports a problem! Consult .las file for "
                 << "further diagnostics!";
      Warning( __FILE__, __LINE__ );
    }

#ifdef PROFILING
    Double t2 = Profiler::GetRealTime();
    (*cla)  << "solution time: " << t2-t1 << " seconds " << std::endl;
    Profiler::WriteReport();
#endif

    // Set invaldiation flag for assembling routine
    assemble_->InvalidateSolBuffer();

  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void SBM_System::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                         const PdeIdType pdeID,
                                         Integer eqnNum,
                                         Double *val ) {

    // Determine sub-matrix
    StdMatrix *stdMat = sysMat_[matrixID]->GetPointer( pdeID, pdeID );

    // Delegate work to implementation in assemble class
    assemble_->AddToDiagMatrixEntry( stdMat, eqnNum,  val );
  }

 // ******************
  //   GetMatrixEntry
  // ******************
  void SBM_System::GetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum,
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum,
                                       Double & val ) {


    Error( "Not implemented", __FILE__, __LINE__ );

  }

  void SBM_System::GetMatrixEntry( FEMatrixType matrixID,
                                       const PdeIdType rowPdeID,
                                       Integer rowEqnNum,
                                       const PdeIdType colPdeID,
                                       Integer colEqnNum2,
                                       Complex & val ) {

    Error( "Not implemented", __FILE__, __LINE__ );
  }

  // ******************
  //   GetMatrixEntry
  // ******************

  void SBM_System::SetMatrixEntry( FEMatrixType matrixID,
                                   const PdeIdType rowPdeID,
                                   Integer rowEqnNum,
                                   const PdeIdType colPdeID,
                                   Integer colEqnNum,
                                   Double val, bool setCounterPart ) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }

  void SBM_System::SetMatrixEntry( FEMatrixType matrixID,
                                   const PdeIdType rowPdeID,
                                   Integer rowEqnNum,
                                   const PdeIdType colPdeID,
                                   Integer colEqnNum,
                                   Complex val, bool setCounterPart ) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }

  // ********************
  //   SetMatrixRowVals
  // ********************
  void SBM_System::SetMatrixRowVals( FEMatrixType matrixID,
                                         const PdeIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Double val ) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }
  void SBM_System::SetMatrixRowVals( FEMatrixType matrixID,
                                         const PdeIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Complex val ) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }

  // ********************
  //   SetMatrixColVals
  // ********************
  void SBM_System::SetMatrixColVals( FEMatrixType matrixID,
                                         const PdeIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Double val ) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }

  void SBM_System::SetMatrixColVals( FEMatrixType matrixID,
                                         const PdeIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Complex val ) {
    Error( "Not implemented", __FILE__, __LINE__ );
  }


  // **********
  //   Export
  // **********
  void SBM_System::Export( FEMatrixType matrixID, Char *filename,
                           Char *comment ) const {
    sysMat_[matrixID]->Export( filename, comment );
  }


  // ****************
  //   CreateSolver
  // ****************
  void SBM_System::CreateSolver(){


    solver_ = GenerateSolverObject( *(sysMat_[SYSTEM]), CG, xml,
                                    &myParams_, &myReport_ );
  }


  // ****************
    //   CreateSolver
    // ****************
    void SBM_System::CreatePrecond(){
      PrecondType precond;
      myParams_.GetEnumValue("Precond", precond);
      precond_ = GenerateSBMPrecondObject( *(sysMat_[SYSTEM]), precond,
                                           &myParams_, &myReport_ );
    }

  // ***************
  //   SetupSolver
  // ***************
  void SBM_System::SetupSolver() {
    solver_->Setup( *sysMat_[SYSTEM] );
  }

}
