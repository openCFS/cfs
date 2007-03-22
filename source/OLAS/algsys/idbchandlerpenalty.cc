// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "algsys/idbchandlerpenalty.hh"
#include "algsys/standardsys.hh"
#include "algsys/baseentrymanipulator.hh"
#include "algsys/generateentrymanipulator.hh"
#include "matvec/matvec.hh"


namespace OLAS {


  // ***********************************
  //   Constructor (StdMatrix systems)
  // ***********************************
  template <typename T>
  IDBC_HandlerPenalty<T>::IDBC_HandlerPenalty( UInt numIDBC, UInt blockSize ){

    ENTER_FCN( "IDBC_HandlerPenalty::IDBC_HandlerPenalty" );

    // Log some information
    (*cla) << "\n IDBC_HandlerPenalty: Administrating "
           << numIDBC << " inhom. Dirichlet BCs.\n\n";

    // This is the StdMatrix case
    sbmCase_ = false;
    offset_  = NULL;

    // Initialise penalty term
    penaltyTerm_ = 0.0;

    // Allocate memory for internal arrays
    numIDBC_ = numIDBC;
    NewArray( dirichletEQN_      , UInt, numIDBC_ );
    NewArray( dirichletComponent_, UInt, numIDBC_ );

#ifdef DEBUG_IDBCHANDLERPENALTY
    // Initialise arrays
    for ( UInt i = 1; i <= numIDBC_; i++ ) {
      dirichletEQN_[i] = 0;
      dirichletComponent_[i] = 0;
    }
#endif

    // Determine entry type from template parameter
    eType_ = NOENTRYTYPE;
    if ( assocType<T>::tagM == assocType<Double>::tagM ) {
      eType_ = DOUBLE;
    }
    else if ( assocType<T>::tagM == assocType<Complex>::tagM ) {
      eType_ = COMPLEX;
    }
    else {
      (*error) << "Internal template error! No swearing please!";
      Error( __FILE__, __LINE__ );
    }

    // Generate vector for storing Dirchlet values
    // NOTE: We use SPARSE_NONSYM as matrix storage type in order to
    //       obtain a Vector<T> object
    dirichletValue_ = GenerateStdVectorObject( SPARSE_NONSYM, eType_,
                                               blockSize, numIDBC_ );

    // Generate our personal assemble object
    assembler_ = GenerateEntryManipulatorObject( eType_, blockSize );
  }


  // *****************************
  //   Constructor (SBM systems)
  // *****************************
  template <typename T>
  IDBC_HandlerPenalty<T>::IDBC_HandlerPenalty( UInt numIDBC, UInt numPDEs,
                                               UInt *bcOffsets ) {

    ENTER_FCN( "IDBC_HandlerPenalty::IDBC_HandlerPenalty" );

    // Log some information
    (*cla) << "\n IDBC_HandlerPenalty: Administrating "
           << numIDBC << " inhom. Dirichlet BCs.\n\n";

    // This is the SBM_Matrix case
    sbmCase_ = true;

    // Initialise penalty term
    penaltyTerm_ = 0.0;

    // Allocate memory for internal arrays
    numIDBC_ = numIDBC;
    NewArray( dirichletEQN_      , UInt, numIDBC_    );
    NewArray( dirichletComponent_, UInt, numIDBC_    );
    NewArray( offset_            , UInt, numPDEs + 1 );

    // Generate array containing index pointers
    //
    // Note: In the algebraic system these are real offsets, while here the
    //       values actually are index pointers, so we must add one to them.
    for ( UInt i = 1; i <= numPDEs + 1; i++ ) {
      offset_[i] = bcOffsets[i];
    }

    // Determine entry type from template parameter
    eType_ = NOENTRYTYPE;
    if ( assocType<T>::tagM == assocType<Double>::tagM ) {
      eType_ = DOUBLE;
    }
    else if ( assocType<T>::tagM == assocType<Complex>::tagM ) {
      eType_ = COMPLEX;
    }
    else {
      (*error) << "Internal template error! No swearing please!";
      Error( __FILE__, __LINE__ );
    }

    // Generate vector for storing Dirchlet values
    StdVector *stdVec = NULL;
    BaseVector *bVec = NULL;
    SBM_Vector *sbmVec = New SBM_Vector( numPDEs );
    UInt aux = 0;
    for ( UInt i = 1; i <= numPDEs; i++ ) {

      // Determine number of IDBCs for i-th PDE
      aux = offset_[i+1] - offset_[i];

      // If there are IDBCs generate and insert sub-vector
      //
      // NOTE: We use SPARSE_NONSYM as matrix storage type in order to
      //       obtain a Vector<T> object
      if ( aux > 0 ) {
        bVec = GenerateStdVectorObject( SPARSE_NONSYM, eType_, 1, aux );
        stdVec = dynamic_cast<StdVector*>(bVec);
        sbmVec->SetSubVector( stdVec, i );
      }
    }
    dirichletValue_ = sbmVec;

#ifdef DEBUG_IDBCHANDLERPENALTY
    // Initialise arrays
    for ( UInt i = 1; i <= numIDBC_; i++ ) {
      dirichletEQN_[i] = 0;
      dirichletComponent_[i] = 0;
    }
    dirichletValue_->Init();
#endif

    // Generate our personal assemble object
    assembler_ = GenerateEntryManipulatorObject( eType_, 1 );
  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  IDBC_HandlerPenalty<T>::~IDBC_HandlerPenalty() {

    ENTER_FCN( "IDBC_HandlerPenalty::~IDBC_HandlerPenalty" );

    // Free memory for internal arrays
    DeleteArray( dirichletEQN_       );
    DeleteArray( dirichletComponent_ );
    DeleteArray( offset_             );

    // Delete vector of Dirichlet values
    delete dirichletValue_;

    // Delete personal assemble object
    delete assembler_;
  }


  // ***********************
  //   InitDirichletValues
  // ***********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::InitDirichletValues() {
    ENTER_FCN( "IDBC_HandlerPenalty::InitDirichletValues" );
    dirichletValue_->Init();
  }


  // *********************
  //   AdaptSystemMatrix
  // *********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::AdaptSystemMatrix( BaseMatrix &sysMat ) {

    ENTER_FCN( "IDBC_HandlerPenalty::AdaptSystemMatrix" );

    // Compute penalty term
    penaltyTerm_ = sysMat.GetMaxDiag() * PENALTY;

    // Check that each IDBC was initialised
#ifdef DEBUG_IDBCHANDLERPENALTY
    for ( UInt i = 1; i <= numIDBC_; i++ ) {
      if ( dirichletEQN_[i] == 0 ) {
        PrintInternalState( std::cerr );
        (*error) << "Found unset IDBC (eqnNo = 0): # = " << i;
        Error( __FILE__, __LINE__ );
      }
      if ( dirichletComponent_[i] == 0 ) {
        PrintInternalState( std::cerr );
        (*error) << "Found unset IDBC (comp = 0): # = " << i;
        Error( __FILE__, __LINE__ );
      }
    }
    PrintInternalState( *debug );
#endif

    // StdMatrix case
    if ( sbmCase_ == false ) {

      TRY_CAST {
        RefCast( sysMat, StdMatrix, stdMat );
        assembler_->AdaptSystemMatrix( stdMat, dirichletEQN_,
                                       dirichletComponent_, numIDBC_,
                                       penaltyTerm_ );
      } CATCH_CAST;
    }

    // SBM_Matrix case
    else {

      TRY_CAST {
        RefCast( sysMat, SBM_Matrix, sbmMat );
        StdMatrix *stdMat = NULL;
        UInt aux = 0;
        for ( UInt i = 1; i <= (UInt)sbmMat.GetNrows(); i++ ) {
          stdMat = sbmMat.GetPointer(i,i);
          aux = offset_[i+1] - offset_[i];
          if ( aux > 0 && stdMat != NULL ) {
            assembler_->AdaptSystemMatrix( *stdMat,
                                           dirichletEQN_ + offset_[i],
                                           dirichletComponent_ + offset_[i],
                                           aux,
                                           penaltyTerm_ );
          }
        }
      } CATCH_CAST;
    }

    // Be verbose
    (*cla) << "\n IDBC_HandlerPenalty:\n"
           << " -> Adapted system matrix\n"
           << " -> Penalty term = "
           << std::scientific << penaltyTerm_ << '\n'
           << std::endl;
    cla->unsetf( std::ios::scientific );
  }


  // ****************
  //   AddIDBCToRHS
  // ****************
  template <typename T>
  void IDBC_HandlerPenalty<T>::AddIDBCToRHS( BaseVector *rhs ) {

    ENTER_FCN( "IDBC_HandlerPenalty::AddIDBCToRHS" );

    // --------------
    // StdMatrix case
    // --------------
    if ( sbmCase_ == false ) {

      TRY_CAST {
        PtrCast( rhs, StdVector, stdVec );
        PtrCast( dirichletValue_, StdVector, stdVal );
        assembler_->AdaptRHSForIDBC( *stdVec,
                                     dirichletEQN_,
                                     dirichletComponent_,
                                     *stdVal,
                                     penaltyTerm_,
                                     numIDBC_ );
      } CATCH_CAST;
    }

    // ---------------
    // SBM_Matrix case
    // ---------------
    else {

      TRY_CAST {

        // First down-cast vectors from Base to SBM
        PtrCast( rhs, SBM_Vector, sbmVec );
        PtrCast( dirichletValue_, SBM_Vector, sbmVal );

        // We need pointers to sub-vectors
        StdVector *stdVec = NULL;
        StdVector *stdVal = NULL;

        // Now loop over all PDEs / sub-vectors and
        // set the Dirichlet values
        UInt aux = 0;
        for ( UInt i = 1; i <= sbmVec->GetSize() ; i++ ) {

          // obtain sub-vectors
          stdVec = sbmVec->GetPointer( i );
          stdVal = sbmVal->GetPointer( i );

          // get number of IDBCs for i-th PDE
          aux = offset_[i+1] - offset_[i];

          // Now adapt the sub-vector
          if ( aux > 0 && stdVec != NULL && stdVal != NULL ) {
            assembler_->AdaptRHSForIDBC( *stdVec,
                                         dirichletEQN_ + offset_[i],
                                         dirichletComponent_ + offset_[i],
                                         *stdVal,
                                         penaltyTerm_,
                                         aux );
          }
        }
      } CATCH_CAST;
    }
  }


  // ***********
  //   SetIDBC
  // ***********
  template <typename T>
  void IDBC_HandlerPenalty<T>::SetIDBC( PdeIdType pdeID, UInt eqnNo,
                                        UInt comp, UInt numBC,
                                        const T &val ) {

    ENTER_FCN( "IDBC_HandlerPenalty::SetIDBC" );

    // StdMatrix case
    if ( sbmCase_ == false ) {
      dirichletValue_->SetVectorEntry( numBC, val );
      dirichletComponent_[numBC] = comp;
      dirichletEQN_      [numBC] = eqnNo;
    }

    // SBM_Matrix case
    else {
      SBM_Vector *sbmVec = dynamic_cast<SBM_Vector*>( dirichletValue_ );
      StdVector *stdVec = sbmVec->GetPointer( pdeID );
      stdVec->SetVectorEntry( numBC, val );
      dirichletComponent_[ numBC + offset_[pdeID] ] = comp;
      dirichletEQN_      [ numBC + offset_[pdeID] ] = eqnNo;
    }
  }


  // *****************
  //   SetDofsToIDBC
  // *****************
  template <typename T>
  void IDBC_HandlerPenalty<T>::SetDofsToIDBC( BaseVector *vec ) {

    ENTER_FCN( "IDBC_HandlerPenalty::SetDofsToIDBC" );

    // -----------------
    // CASE 1: StdMatrix
    // -----------------
    if ( sbmCase_ == false ) {

      // Down-cast vector
      StdVector *stdVec = dynamic_cast<StdVector*>( vec );
      if ( stdVec == NULL ) {
        Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
      }

      // Insert values
      T aux;
      for ( UInt i = 1; i <= numIDBC_; i++ ) {
        dirichletValue_->GetVectorEntry( i, aux );
        assembler_->SetVectorEntry( stdVec, dirichletEQN_[i],
                                    dirichletComponent_[i], aux );
      }
    }

    // ------------------
    // CASE 2: SBM_Matrix
    // ------------------
    else {

      // Down-cast vector
      SBM_Vector *sbmVec = dynamic_cast<SBM_Vector*>( vec );
      SBM_Vector *sbmVal = dynamic_cast<SBM_Vector*>( dirichletValue_ );
      if ( sbmVec == NULL || sbmVal == NULL ) {
        Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
      }

      // Loop over all sub-vectors
      StdVector *stdVec = NULL;
      StdVector *stdVal = NULL;
      T aux;
      UInt idx;

      for ( UInt i = 1; i <= sbmVec->GetSize(); i++ ) {

#ifdef DEBUG_IDBCHANDLERPENALTY
        (*debug) << " IDBC_HandlerPenalty::SetDofsToIDBC:"
                 << " Setting " << offset_[i+1] - offset_[i]
                 << " IDBCs for PDE #" << i << std::endl;
#endif

        stdVec = sbmVec->GetPointer( i );
        stdVal = sbmVal->GetPointer( i );

        // Insert values
        for ( UInt j = 1; j <= offset_[i+1] - offset_[i]; j++ ) {
          stdVal->GetVectorEntry( j, aux );
          idx = offset_[i] + j;
          assembler_->SetVectorEntry( stdVec, dirichletEQN_[idx],
                                      dirichletComponent_[idx], aux );
        }
      }
    }
  }


  // ****************************
  //   InstantiatePublicMethods
  // ****************************
  template <typename T>
  void IDBC_HandlerPenalty<T>::InstantiatePublicMethods() {

    ENTER_FCN( "IDBC_HandlerPenalty::InstantiatePublicMethods" );

    Error( "This function should never be called", __FILE__, __LINE__ );

    // Some dummy variables we need
    BaseVector *dummyVec = NULL;
    std::map<FEMatrixType,Double> dummyMap;
    BaseMatrix *dummyMat = NULL;
    T aux(0.0);

    // Public methods in alphabetical order
    AdaptSystemMatrix( *dummyMat );
    AddIDBCToRHS( dummyVec );
    BuiltSystemMatrix( dummyMap );
    InitDirichletValues();
    InitMatrix( SYSTEM );
    RemoveIDBCFromRHS( dummyVec );
    SetDofsToIDBC( dummyVec );
    SetIDBC( NO_PDE_ID, 1, 0, 0, aux );
  }


  // **********************
  //   PrintInternalState
  // **********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::PrintInternalState( std::ostream &os ) {

    ENTER_FCN( "IDBC_HandlerPenalty::PrintInternalState" );

    // ----------------
    // CASE1: StdSystem
    // ----------------
    if ( sbmCase_ == false ) {

      os << " IDBC_HandlerPenalty --- Internal State:\n"
         << " --------------------------------------\n"
         << " numIDBC_    = " << numIDBC_ << std::endl
         << " # | eqnNo | comp | value\n";

      T val;
      PtrCast( dirichletValue_, StdVector, stdVal );

      for ( UInt i = 1; i <= numIDBC_; i++ ) {
        stdVal->GetVectorEntry( i, val );
        os << i << " | "
           << dirichletEQN_[i] << " | "
           << dirichletComponent_[i] << " | "
           << val << std::endl;
      }
    }


    // -----------------
    // CASE2: SBM_System
    // -----------------
    else {

      os << " IDBC_HandlerPenalty --- Internal State:\n"
         << " --------------------------------------\n"
         << " numIDBC_ = " << numIDBC_ << '\n' << std::endl;

      StdVector *stdVal = NULL;
      PtrCast( dirichletValue_, SBM_Vector, sbmVec );
      UInt aux = 0;
      UInt idx = 0;
      T val;

      // Compute column widths for pretty-printing
      aux = numIDBC_ + 1;
      UInt cw1 = aux > 0 ? (UInt)std::log10( (float)aux ) + 1 : 1;
      UInt cw2 = 5;
      UInt cw3 = 4;

      for ( UInt i = 1; i <= dirichletValue_->GetSize(); i++ ) {

        // Compute number of IDBCs
        aux = offset_[i+1] - offset_[i];

        // Get hold of sub-vector
        stdVal = sbmVec->GetPointer(i);

        // Print IDBCs for this PDE
        os << " PDE # " << i << " / numIDBCs = " << aux << std::endl;
        if ( aux > 0 ) {
          os << " # | eqnNo | comp | value" << std::endl;
          for ( UInt j = 1; j <= aux; j++ ) {

            stdVal->GetVectorEntry( j, val );
            idx = offset_[i] + j;

            os << std::setfill( ' ' )
               << std::setw(cw1) << idx << " | "
               << std::setw(cw2) << dirichletEQN_[idx] << " | "
               << std::setw(cw3) << dirichletComponent_[idx] << " | "
               << val << std::endl;
          }
        }
        os << std::endl;
      }
    }

    os << " --------------------------------------\n" << std::endl;

  }

}
