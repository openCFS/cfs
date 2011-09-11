// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/algorithm/string/replace.hpp>
#include "MatVec/sbmmatrix.hh"
#include "MatVec/generatematvec.hh"

#include "OLAS/algsys/sbmsystem.hh"

#include "OLAS/graph/graphmanagersbmmat.hh"
#include "OLAS/precond/generateprecond.hh"
#include "OLAS/precond/baseprecond.hh"
#include "OLAS/solver/basesolver.hh"
#include "OLAS/solver/generatesolver.hh"

#include "OLAS/algsys/baseentrymanipulator.hh"
#include "OLAS/algsys/generateidbchandler.hh"
#include "OLAS/algsys/baseidbchandler.hh"

#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  SBM_System::SBM_System(PtrParamNode pn) : BaseSystem(pn) {


    precond_             = NULL;
    rhs_                 = NULL;
    sol_                 = NULL;

    // We are an SBM_System obviously
    algSysType_          = SBM_SYSTEM;

    // Initialize pointer to finite-element matrices
    sysMat_.Resize( MAX_NUM_FE_MATRICES );
    sysMat_.Init( NULL );

    // Assuming non-symmetric matrices by default
    sbmSymm_ = false;

    // Generate array for storing FE-matrix pattern information
    feSubMatrices_.Resize( MAX_NUM_FE_MATRICES );
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
    for ( UInt i = 0; i < MAX_NUM_FE_MATRICES; i++ ) {
      delete sysMat_[i];
    }
    sysMat_.Clear();
    sysMat_.Clear();
  }


  // ****************
  //   SetupPrecond
  // ****************
  void SBM_System::SetupPrecond() {
    precond_->Setup( *sysMat_[SYSTEM]);
  }


  // ****************
  //   SetBlockSize
  // ****************
  void SBM_System::SetBlockSize( const FeFctIdType identifier, const UInt bs ) {


    // Check consistency
    if ( bs > 1 ) {
      EXCEPTION( "SBM_System::SetBlockSize: "
               << "We do not support block-size = " << bs
               << ", but only scalar entries!");
    }
    else {
      blockSize_ = 1;
    }
  }


  void SBM_System::GetSolutionVal( SingleVector& ptSol,
                                   const FeFctIdType identifierPDE ) {

    if( ptSol.GetEntryType() == BaseMatrix::COMPLEX ) {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>& >( ptSol );
      retVec.Resize( size_ );
      UInt index = 0;
      for( UInt i = 0; i < numPDEs_; i++ ) {
        const Vector<Complex> & dVec1 = dynamic_cast<Vector<Complex>&>( (*sol_)(i));
        for( UInt j = 0; j < dVec1.GetSize(); j++ ) {
          retVec[index++] = dVec1[j];
        }
      }
    }  else {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>& >( ptSol );
      retVec.Resize( size_ );
      UInt index = 0;
      for( UInt i = 0; i < numPDEs_; i++ ) {
        const Vector<Double> & dVec1 = dynamic_cast<Vector<Double>&>( (*sol_)(i));
        for( UInt j = 0; j < dVec1.GetSize(); j++ ) {
          retVec[index++] = dVec1[j];
        }
      }
      std::cerr << "retVec = " << retVec << std::endl;
    }
  }


  void SBM_System::GetRHSVal( SingleVector &ptRhs,
                              const FeFctIdType identifierPDE ) {

    //Warning( "SBM_System::GetRHSVal: At the moment we are not respecting the splitting"
    //         " of the result vecto and return just the complete vector, regardless of the pdeId ");

    if( ptRhs.GetEntryType() == BaseMatrix::COMPLEX ) {
      Vector<Complex> & retVec = dynamic_cast<Vector<Complex>& >( ptRhs );
      retVec.Resize( size_ );
      UInt index = 0;
      for( UInt i = 0; i < numPDEs_; i++ ) {
        const Vector<Complex> & dVec1 = dynamic_cast<Vector<Complex>&>( (*rhs_)(i));
        for( UInt j = 0; j < dVec1.GetSize(); j++ ) {
          retVec[index++] = dVec1[j];
        }
      }
    }  else {
      Vector<Double> & retVec = dynamic_cast<Vector<Double>& >( ptRhs );
      retVec.Resize( size_ );
      UInt index = 0;
      for( UInt i = 0; i < numPDEs_; i++ ) {
        const Vector<Double> & dVec1 = dynamic_cast<Vector<Double>&>( (*rhs_)(i));
        for( UInt j = 0; j < dVec1.GetSize(); j++ ) {
          retVec[index++] = dVec1[j];
        }
      }

    }
  }

  // *******************
  //   SetFEMatrixType
  // *******************
  void SBM_System::SetFEMatrixType( const FEMatrixType matType,
                                    const FeFctIdType pdeID1,
                                    const FeFctIdType pdeID2 ) {


    // Determine row and column index of sub-matrix
    UInt rowInd = pdeID1;
    UInt colInd = pdeID2 == NO_FCT_ID ? pdeID1 : pdeID2;

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
    sbmSymm_ = false;    
    PtrParamNode matrixNode;
    
    matrixNode = xml_->Get("sbmMatrix", ParamNode::INSERT);
    matrixNode->GetValue("symmetric", sbmSymm_, ParamNode::INSERT);

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
      for ( UInt i = 0; i < numPDEs_; i++ ) {
        for ( UInt j = 0; j < numPDEs_; j++ ) {
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
    for( UInt i = 0; i < numPDEs_; i++ ) {
      size_ += numLastFreeDof_[i];
    }

    // ------------------------------
    //  Generation of matrix objects
    // ------------------------------

    // Log what we will do
    PrintFeMatrixInfo( cla );

    // Obtain some info from parameter file
    BaseMatrix::EntryType entryType;
    std::string entryStr = "double";
    matrixNode->GetValue("entry", entryStr, ParamNode::INSERT);
    entryType = BaseMatrix::entryType.Parse( entryStr );

    for ( fIt = matrixTypes_.begin(); fIt != matrixTypes_.end(); fIt++ ) {
      sysMat_[*fIt] = GenerateSBM_Matrix( *fIt, entryType );
    }


    // --------------------------------------------
    //  Treatment of Dirichlet Boundary Conditions
    // --------------------------------------------

    // Check the case we have (elimination vs. penalty) and generate
    // an appropriate object for handling inhomogeneous Dirichlet
    // boundary conditions
    PtrParamNode setupNode;
    std::string idbcHandlingStr = "penalty";
    bool usePenaltyFormulation = true;
    
    setupNode = xml_->Get("setup", ParamNode::INSERT);
    setupNode->GetValue("idbcHandling", idbcHandlingStr, ParamNode::INSERT);
    if(idbcHandlingStr == "penalty") {
      usePenaltyFormulation = true;
    } else {
      usePenaltyFormulation = false;
    }
    
    if ( usePenaltyFormulation ) {
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
      EXCEPTION( WRONG_CAST_MSG );
    }

    // For the moment we insert a sub-vector for each position.
    // In the case of the right-hand side we might actually be
    // more economic. How do we get the information which sub-vectors
    // are really needed, however?
    StdMatrix *stdMat = NULL;
    BaseVector *bVec = NULL;
    SingleVector *sVec = NULL;
    for ( UInt k =0; k < numPDEs_; k++ ) {

      // Get diag matrix for vector generation
      stdMat = sysMat_[SYSTEM]->GetPointer( k, k );

      // Insert sub-vector into solution
      bVec = GenerateVectorObject( *stdMat );
      sVec = dynamic_cast<SingleVector*>( bVec );
      sol_->SetSubVector( sVec, k );

      // Insert sub-vector into right-hand side
      bVec = GenerateVectorObject( *stdMat );
      sVec = dynamic_cast<SingleVector*>( bVec );
      rhs_->SetSubVector( sVec, k );
    }

    // ----------------------------------
    //  Generation of auxilliary objects
    // ----------------------------------

    // Create the assemble object
    assemble_ = new BaseEntryManipulator();

    // Generate SingleVectors

    BaseMatrix::StorageType sType = stdMat->GetStorageType();
    solBuffer_ = dynamic_cast<SingleVector*>
        ( GenerateSingleVectorObject( sType, entryType, size_ ) );
    rhsBuffer_ = dynamic_cast<SingleVector*>
        ( GenerateSingleVectorObject( sType, entryType, size_ ) );

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
      std::string tmp;
      Enum2String( *fIt, tmp );

      (*os) << " " << tmp << std::endl;

      // Now plot sub-matrix pattern
      SubMatrixID sID;
      for ( UInt i = 0; i < numPDEs_; i++ ) {
        sID.rowInd = i;
        for ( UInt j = 0; j < numPDEs_; j++ ) {
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
		  BaseMatrix::EntryType entryType ) {


    // STEP 1: Generate empty SBM_Matrix
    SBM_Matrix *retMat = NULL;
    retMat = new SBM_Matrix( numPDEs_, numPDEs_, sbmSymm_ );
    if ( retMat == NULL ) {
      EXCEPTION( "SBM_System::GenerateSBM_Matrix: "
               << "This is the end my friend!\n"
               << "Generation of empty SBM_Matrix failed!" );
    }

    // STEP 2: Populate with sub-matrices
    std::set<SubMatrixID,SortSubMatrixID>::iterator sIt;
    BaseGraph *graph;
    for ( sIt = feSubMatrices_[matType].begin();
          sIt != feSubMatrices_[matType].end(); sIt++ ) {

      // Determine associated PDE identifiers
      FeFctIdType pde1 = (*sIt).rowInd;
      FeFctIdType pde2 = (*sIt).colInd;

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
          retMat->SetSubMatrix ( pde1, pde2, entryType, BaseMatrix::SPARSE_SYM,
                                 nrows, ncols, graph->GetNNE() );
        }
        else {
          retMat->SetSubMatrix ( pde1, pde2, entryType, BaseMatrix::SPARSE_NONSYM,
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

  void SBM_System::SetElementMatrix( FEMatrixType matrix_id, 
                                          const Matrix<Double>& elemMat,
                                          FeFctIdType idPDE1,
                                          const StdVector<Integer>& eqnNrs1,
                                          FeFctIdType idPDE2,
                                          const StdVector<Integer>& eqnNrs2,
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
        EXCEPTION("SBM_System::SetElementMatrix:\n"
                  << " idPDE1 ........... " << idPDE1 << '\n'
                  << " idPDE2 ........... " << idPDE2 << '\n'
                  << " setCounterPart ... " << std::boolalpha << setCounterPart
                  << " sbmSymm .......... true\n"
                  << " leads to loss of element matrix!" );
      }

      // Matrix in lower part and counter part should be set,
      // so we change to the counterPart
      else {
        setCounterPart = false;
        reverse = true;
      }
    }

    // Extract associated sub-matrix
    StdMatrix *subMat1 = &((*sysMat_[matrix_id])( idPDE1, idPDE2 ));
    StdMatrix *subMat2 = &((*sysMat_[matrix_id])( idPDE2, idPDE1 ));

    // Delegate remaining work to assemble
    if ( reverse == false ) {
      assemble_->SetElementMatrix( matrix_id, idPDE1, idPDE2,
                                   subMat1, subMat2,
                                   idbcHandler_, elemMat,
                                   eqnNrs1,
                                   eqnNrs2,  
                                   numLastFreeDof_[idPDE1],
                                   numLastFreeDof_[idPDE2],
                                   setCounterPart);
    }
    else {
      assemble_->SetCounterPartOnly( matrix_id, idPDE1, idPDE2,
                                     subMat1, subMat2,
                                     idbcHandler_, elemMat,
                                     eqnNrs1,
                                     eqnNrs2,
                                     numLastFreeDof_[idPDE1],
                                     numLastFreeDof_[idPDE2] );
    }
  }
  void SBM_System::SetElementMatrix(FEMatrixType matrix_id, 
                                    const Matrix<Complex>& elemmat,
                                    FeFctIdType identifierPDE1,
                                    const StdVector<Integer>& eqnNrs1,
                                    FeFctIdType identifierPDE2,
                                    const StdVector<Integer>& eqnNrs2,
                                    bool setCounterPar ) {
    WARN( "Adapt me SetElementMatrix(Complex)");
  }

  // **************
  //   InitMatrix
  // **************
  void SBM_System::InitMatrix( FEMatrixType matrixID,
                               const FeFctIdType pdeID ) {


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

    PtrParamNode setupNode;
    std::string idbcHandlingStr = "penalty";
    bool usePenaltyFormulation = true;
    
    setupNode = xml_->Get("setup", ParamNode::INSERT);
    setupNode->GetValue("idbcHandling", idbcHandlingStr, ParamNode::INSERT);
    if(idbcHandlingStr == "penalty") {
      usePenaltyFormulation = true;
    } else {
      usePenaltyFormulation = false;
    }

    // Do we need to re-insert penalty values into matrix?
    if ( ( matrixID == NOTYPE || matrixID == SYSTEM ) && usePenaltyFormulation ) {
      assembleDirichletToSysMat_ = true;
    }
  }


  // ***********************
  //   InitRHS (Version 1)
  // ***********************
  void SBM_System::InitRHS( FeFctIdType pdeID ) {


    if ( pdeID == NO_FCT_ID ) {
      rhs_->Init();
    }
    else {
      rhs_->GetPointer( pdeID )->Init();
    }
  }


  // ***********
  //   InitSol
  // ***********
  void SBM_System::InitSol( FeFctIdType pdeID ) {


    if ( pdeID == NO_FCT_ID ) {
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
        WARN("SBM_System::ConstructEffectiveMatrix: "
             << "Map with factors is empty, but there are "
             << matrixTypes_.size() << " FE matrices in the game!");
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
  void SBM_System::SetElementRHS( const Vector<Double>& elemRHS, 
                                  const FeFctIdType idPDE,
                                  StdVector<Integer>& eqnNrs ) {
    WARN( "adapt me" );
//
//    // Delegate work to EntryManipulator
//    assemble_->SetElementRHS( rhs_->GetPointer(idPDE), elemRHS, connect,
//                              length, numLastFreeDof_[idPDE] );
  }
  void SBM_System::SetElementRHS( const Vector<Complex>& elemRHS, 
                                  const FeFctIdType idPDE,
                                  StdVector<Integer>& eqnNrs ) {
    WARN( "adapt me");
    //
    //    // Delegate work to EntryManipulator
    //    assemble_->SetElementRHS( rhs_->GetPointer(idPDE), elemRHS, connect,
    //                              length, numLastFreeDof_[idPDE] );
  }


  // *************
  //   UpdateRHS
  // *************
  void SBM_System::UpdateRHS( FEMatrixType matrix_id, const BaseVector& fup ) {
    EXCEPTION( "Method not yet implemented! -> UpdateRHS" );
  }


  // *************************
  //   SetDirichlet (Double)
  // *************************
  void SBM_System::SetDirichlet(const FeFctIdType pdeID,
                                Integer eqnNum, const Double &val ) {


    // Perform some consistency checks
#ifdef DEBUG_SBMSYSTEM

    // Some situation dependend modifiers
    UInt maxVal = usingPenalty_ == true ? sizePerPDE_[pdeID] :
      numLastFreeDof_[pdeID];
    UInt minVal = usingPenalty_ == true ? 1 : numLastFreeDof_[pdeID] + 1;

    std::string sysName;
    xml_->Get("name", sysName, false);

    if ( eqnNum > maxVal || eqnNum < minVal  ) {
      EXCEPTION( "SBMSystem::SetDirichlet: Inconsistency detected:"
                 << "\n pdeID           = " << pdeID
                 << "\n eqnNum          = " << eqnNum
                 << "\n val             = " << val
                 << "\n numBC           = " << numDirichletValues_
                 << "\n size            = " << sizePerPDE_[pdeID]
                 << "\n numLastFreeDof_ = " << numLastFreeDof_[pdeID]
                 << "\n minVal          = " << minVal
                 << "\n maxVal          = " << maxVal
                 << "\n SystemName is '"
                 << sysName << "'";
    }
#endif

    // Delegate work to IDBC handler
    idbcHandler_->SetIDBC( pdeID, eqnNum,  val );
  }


  // **************************
  //   SetDirichlet (Complex)
  // **************************
  void SBM_System::SetDirichlet( const FeFctIdType pdeID,
                                 Integer eqnNum, const Complex &val ) {


    // Perform some consistency checks
#ifdef DEBUG_SBMSYSTEM

    // Some situation dependend modifiers
    UInt maxVal = sizePerPDE_[pdeID];
    UInt minVal = usingPenalty_ == true ? 1 : numLastFreeDof_[pdeID] + 1;

    std::string sysName;
    xml_->Get("name", sysName, false);

    if ( eqnNum > maxVal || eqnNum < minVal ) {
      EXCEPTION( "SBMSystem::SetDirichlet: Inconsistency detected:"
                 << "\n pdeID           = " << pdeID
                 << "\n eqnNum          = " << eqnNum
                 << "\n val             = " << val
                 << "\n numBC           = " << numDirichletValues_
                 << "\n sizePerPDE      = " << sizePerPDE_[pdeID]
                 << "\n numLastFreeDof_ = " << numLastFreeDof_[pdeID]
                 << "\n minVal          = " << minVal
                 << "\n maxVal          = " << maxVal
                 << "\n SystemName is '"
                 << sysName << "'" )
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
  void SBM_System::Solve(PtrParamNode analysis_id) {


    // If the penalty formulation is used and we have inhomogeneous
    // Dirichlet boundary conditions, then the righ-hand side is
    // "contaminated" with penalty terms
    if ( numDirichletValues_ > 0 && usingPenalty_ ) {
      solver_->SetUsingPenalty( false );
    }


    // Iterative solvers require an initial guess and in the penalty case
    // we should insert the Dirichlet values into it
    if ( dynamic_cast<BaseIterativeSolver*>(solver_) != NULL &&
         usingPenalty_ ) {
    idbcHandler_->SetDofsToIDBC( sol_ );
    (*cla) << " Inserted Dirichlet values into initial guess"
           << std::endl;
    }


    // Assume that everything will go well
    PtrParamNode out = olasInfo_->Get(ParamNode::PROCESS)->Get("solver");
    out->Get("solutionIsOkay")->SetValue(true);

    // Now modifiy the right-hand side vector
    idbcHandler_->AddIDBCToRHS( rhs_ );

    // Remove the export linear system stuff, it has changed in standardsys.cc an as below
    // the solve part is commentet out, I see no reason to export linsys also here, it would
    // require a generalization anyway. Fabian 16.11.07
    // check if we do export stuff
     PtrParamNode els = xml_->Get("exportLinSys", ParamNode::PASS );
     std::string file;
     std::string base;

     // TODO: This is most ugly copy & paste from standardsys.cc -> Generalize common parts!!
     // need it common even when exclusive solution
     if(els) {
       std::ostringstream os;
       std::string name = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
       os << name;
       std::string id = analysis_id->Get("analysis_id")->As<std::string>();
       boost::replace_all(id, ":", "_");
       os << "_" << id;
       base = os.str();
     }

     // check if we do not only want the solution
     if(els && els->Get("solution")->As<std::string>() != "exclusive")
     {
       // two formats. The harwell-boing format includes the rhs!
       if(els->Get("format")->As<std::string>()  == "harwell-boeing")
       {
         EXCEPTION( "Harwell-Boeing Format not implemented for SBM-case" );
       }
       else // classical (default) matrix-market
       {
         sysMat_[SYSTEM]->Export((base+".mtx").c_str() );

         if(els->HasByVal("damping", true) && sysMat_[DAMPING] != NULL)
           sysMat_[DAMPING]->Export((base+"_damping.mtx").c_str() );

         if(els->HasByVal("auxiliary", true) && sysMat_[AUXILIARY] != NULL)
           sysMat_[AUXILIARY]->Export((base+"_aux.mtx").c_str() );

         // rhs is only in harwell-boing included
         rhs_->Export((base+".vec").c_str() );
       }
     }
     if(els && els->HasByVal("initialGuess", true))
           sol_->Export((base+"_intial_guess.vec").c_str());

    // Trigger solution
    solver_->Solve( *sysMat_[SYSTEM], *precond_, *rhs_, *sol_ );

    // Now de-modifiy the right-hand side vector
    idbcHandler_->RemoveIDBCFromRHS( rhs_ );

    // Check that solution went fine, if not issue a warning
    if ( out->Get("solutionIsOkay")->As<bool>() == false ) {
      WARN("Solver reports a problem! Consult .las file for "
           << "further diagnostics!");
    }
  }


  // ************************
  //   AddToDiagMatrixEntry
  // ************************
  void SBM_System::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                         const FeFctIdType pdeID,
                                         Integer eqnNum,
                                         Double val  ) {
    WARN( "Adapt me");
//    // Determine sub-matrix
//    StdMatrix *stdMat = sysMat_[matrixID]->GetPointer( pdeID, pdeID );
//
//    // Delegate work to implementation in assemble class
//    assemble_->AddToDiagMatrixEntry( stdMat, eqnNum,  val );
  }
  
  void SBM_System::AddToDiagMatrixEntry( FEMatrixType matrixID,
                                          const FeFctIdType pdeID,
                                          Integer eqnNum,
                                          Complex val  ) {
     WARN( "Adapt me");
 //    // Determine sub-matrix
 //    StdMatrix *stdMat = sysMat_[matrixID]->GetPointer( pdeID, pdeID );
 //
 //    // Delegate work to implementation in assemble class
 //    assemble_->AddToDiagMatrixEntry( stdMat, eqnNum,  val );
   }
  

 // ******************
  //   GetMatrixEntry
  // ******************
  void SBM_System::GetMatrixEntry( FEMatrixType matrixID,
                                       const FeFctIdType rowPdeID,
                                       Integer rowEqnNum,
                                       const FeFctIdType colPdeID,
                                       Integer colEqnNum,
                                       Double & val ) {


    EXCEPTION( "Not implemented" );

  }

  void SBM_System::GetMatrixEntry( FEMatrixType matrixID,
                                       const FeFctIdType rowPdeID,
                                       Integer rowEqnNum,
                                       const FeFctIdType colPdeID,
                                       Integer colEqnNum2,
                                       Complex & val ) {

    EXCEPTION( "Not implemented" );
  }

  // ******************
  //   GetMatrixEntry
  // ******************

  void SBM_System::SetMatrixEntry( FEMatrixType matrixID,
                                   const FeFctIdType rowPdeID,
                                   Integer rowEqnNum,
                                   const FeFctIdType colPdeID,
                                   Integer colEqnNum,
                                   Double val, bool setCounterPart ) {
    EXCEPTION( "Not implemented" );
  }

  void SBM_System::SetMatrixEntry( FEMatrixType matrixID,
                                   const FeFctIdType rowPdeID,
                                   Integer rowEqnNum,
                                   const FeFctIdType colPdeID,
                                   Integer colEqnNum,
                                   Complex val, bool setCounterPart ) {
    EXCEPTION( "Not implemented" );
  }

  // ********************
  //   SetMatrixRowVals
  // ********************
  void SBM_System::SetMatrixRowVals( FEMatrixType matrixID,
                                         const FeFctIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Double val ) {
    EXCEPTION( "Not implemented" );
  }
  void SBM_System::SetMatrixRowVals( FEMatrixType matrixID,
                                         const FeFctIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Complex val ) {
    EXCEPTION( "Not implemented" );
  }

  // ********************
  //   SetMatrixColVals
  // ********************
  void SBM_System::SetMatrixColVals( FEMatrixType matrixID,
                                         const FeFctIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Double val ) {
    EXCEPTION( "Not implemented" );
  }

  void SBM_System::SetMatrixColVals( FEMatrixType matrixID,
                                         const FeFctIdType pdeID,
                                         Integer eqnNum, UInt dof,
                                         Complex val ) {
    EXCEPTION( "Not implemented" );
  }


  // **********
  //   Export
  // **********
  void SBM_System::Export( FEMatrixType matrixID, char *filename,
                           char *comment ) const {
    sysMat_[matrixID]->Export( filename, comment );
  }


  // ****************
  //   CreateSolver
  // ****************
  void SBM_System::CreateSolver(){

    // HARD CODED: Create conjugate gradient solver
    WARN( "At the moment we use a hard-coded CG-solver" );
    solver_ = GenerateSolverObject( *(sysMat_[SYSTEM]), xml_, olasInfo_);
  }


  void SBM_System::CreatePrecond(){
    precond_ = GenerateSBMPrecondObject( *(sysMat_[SYSTEM]), xml_, olasInfo_ );
  }

  // ***************
  //   SetupSolver
  // ***************
  void SBM_System::SetupSolver(PtrParamNode analysis_id) {
    solver_->Setup( *sysMat_[SYSTEM]);
  }


  void SBM_System::CreateEigenSolver() {
    WARN("SBM_System::CreateEigenSolver not yet implemented!");
  }

  void SBM_System::SetupEigenSolver( UInt numFreq, Double shift, bool quadratic ) {
    WARN("SBM_System::SetupEigenSolver not yet implemented!");
  }
  
  void SBM_System::CalcEigenFrequencies( Vector<Complex>& frequencies,
                                         Vector<Double>& err ) {
    WARN("SBM_System::CalcEigenFrequencies not yet implemented!");
  }

  void SBM_System::CalcEigenFrequencies( Vector<Double>& frequencies,
                                         Vector<Double>& err ) {
    WARN("SBM_System::CalcEigenFrequencies not yet implemented!");
  }

  void SBM_System::CalcEigenMode( UInt numMode ) {
    WARN("SBM_System::CalcEigenMode not yet implemented!");
  }

}
