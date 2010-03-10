// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/sbmmatrix.hh"
#include "MatVec/sbmvector.hh"

#include "OLAS/algsys/standardsys.hh"
#include "OLAS/algsys/baseentrymanipulator.hh"

#include "idbchandlerpenalty.hh"

namespace CoupledField {


  // ***********************************
  //   Constructor (StdMatrix systems)
  // ***********************************
  template <typename T>
  IDBC_HandlerPenalty<T>::IDBC_HandlerPenalty( UInt numIDBC, UInt blockSize ){


    // Log some information
    (*cla) << "\n IDBC_HandlerPenalty: Administrating "
           << numIDBC << " inhom. Dirichlet BCs.\n\n";

    // This is the StdMatrix case
    sbmCase_ = false;
    offset_  = NULL;

    // Initialise next index for boundary condition array
    nextIndex_ = 0;

    // Initialise penalty term
    penaltyTerm_ = 0.0;

    // Allocate memory for internal arrays
    numIDBC_ = numIDBC;
    NEWARRAY( dirichletEQN_, UInt, numIDBC_ );

#ifdef DEBUG_IDBCHANDLERPENALTY
    // Initialise arrays
    for ( UInt i = 0; i < numIDBC_; i++ ) {
      dirichletEQN_[i] = 0;
    }
#endif

    // Determine entry type from template parameter
    eType_ = BaseMatrix::NOENTRYTYPE;
    if ( AssocType<T>::tagM == AssocType<Double>::tagM ) {
      eType_ = BaseMatrix::DOUBLE;
    }
    else if ( AssocType<T>::tagM == AssocType<Complex>::tagM ) {
      eType_ = BaseMatrix::COMPLEX;
    }
    else {
      EXCEPTION( "Internal template error! No swearing please!");
    }

    // Generate vector for storing Dirchlet values
    // NOTE: We use SPARSE_NONSYM as matrix storage type in order to
    //       obtain a Vector<T> object
    dirichletValue_ = GenerateSingleVectorObject( BaseMatrix::SPARSE_NONSYM, eType_,
                                                  numIDBC_ );

    // Generate our personal assemble object
    assembler_ = new BaseEntryManipulator();
  }


  // *****************************
  //   Constructor (SBM systems)
  // *****************************
  template <typename T>
  IDBC_HandlerPenalty<T>::IDBC_HandlerPenalty( UInt numIDBC, UInt numPDEs,
                                               UInt *bcOffsets ) {


    // Log some information
    (*cla) << "\n IDBC_HandlerPenalty: Administrating "
           << numIDBC << " inhom. Dirichlet BCs.\n\n";

    // This is the SBM_Matrix case
    sbmCase_ = true;

    // Initialise penalty term
    penaltyTerm_ = 0.0;

    // Allocate memory for internal arrays
    numIDBC_ = numIDBC;
    NEWARRAY( dirichletEQN_      , UInt, numIDBC_    );
    NEWARRAY( offset_            , UInt, numPDEs + 1 );

    // Generate array containing index pointers
    //
    // Note: In the algebraic system these are real offsets, while here the
    //       values actually are index pointers, so we must add one to them.
    for ( UInt i = 0; i < numPDEs + 1; i++ ) {
      offset_[i] = bcOffsets[i];
    }

    // Determine entry type from template parameter
    eType_ = BaseMatrix::NOENTRYTYPE;
    if ( AssocType<T>::tagM == AssocType<Double>::tagM ) {
      eType_ = BaseMatrix::DOUBLE;
    }
    else if ( AssocType<T>::tagM == AssocType<Complex>::tagM ) {
      eType_ = BaseMatrix::COMPLEX;
    }
    else {
      EXCEPTION( "Internal template error! No swearing please!" );
    }

    // Generate vector for storing Dirchlet values
    SingleVector *stdVec = NULL;
    BaseVector *bVec = NULL;
    SBM_Vector *sbmVec = new SBM_Vector( numPDEs );
    UInt aux = 0;
    for ( UInt i = 0; i < numPDEs; i++ ) {

      // Determine number of IDBCs for i-th PDE
      aux = offset_[i+1] - offset_[i];

      // If there are IDBCs generate and insert sub-vector
      //
      // NOTE: We use SPARSE_NONSYM as matrix storage type in order to
      //       obtain a Vector<T> object
      if ( aux > 0 ) {
        bVec = GenerateSingleVectorObject( BaseMatrix::SPARSE_NONSYM, eType_, aux );
        stdVec = dynamic_cast<SingleVector*>(bVec);
        sbmVec->SetSubVector( stdVec, i );
      }
    }
    dirichletValue_ = sbmVec;

#ifdef DEBUG_IDBCHANDLERPENALTY
    // Initialise arrays
    for ( UInt i = 0; i < numIDBC_; i++ ) {
      dirichletEQN_[i] = 0;
    }
    dirichletValue_->Init();
#endif

    // Generate our personal assemble object
    assembler_ = new BaseEntryManipulator();
  }


  // **************
  //   Destructor
  // **************
  template <typename T>
  IDBC_HandlerPenalty<T>::~IDBC_HandlerPenalty() {


    // Free memory for internal arrays
    delete [] ( dirichletEQN_ );
    delete [] ( offset_ );

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
    dirichletValue_->Init();
  }


  // *********************
  //   AdaptSystemMatrix
  // *********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::AdaptSystemMatrix( BaseMatrix &sysMat ) {


    // Compute penalty term
    penaltyTerm_ = sysMat.GetMaxDiag() * 1e12;

    // Check that each IDBC was initialised
#ifdef DEBUG_IDBCHANDLERPENALTY
    for ( UInt i = 0; i < numIDBC_; i++ ) {
      if ( dirichletEQN_[i] == 0 ) {
        PrintInternalState( std::cerr );
        EXCEPTION("Found unset IDBC (eqnNo = 0): # = " << i);
      }
    }
//    PrintInternalState( *debug );
#endif

    // StdMatrix case
    if ( sbmCase_ == false ) {
      StdMatrix& stdMat = dynamic_cast<StdMatrix&>(sysMat);
      assembler_->AdaptSystemMatrix( stdMat, dirichletEQN_,
          numIDBC_,
          penaltyTerm_ );
    }

    // SBM_Matrix case
    else {
      SBM_Matrix& sbmMat = dynamic_cast<SBM_Matrix&>(sysMat);
        StdMatrix *stdMat = NULL;
        UInt aux = 0;
        const unsigned int srows(sbmMat.GetNumRows());
        for ( UInt i = 0; i < srows; i++ ) {
          stdMat = sbmMat.GetPointer(i,i);
          aux = offset_[i+1] - offset_[i];
          if ( aux > 0 && stdMat != NULL ) {
            assembler_->AdaptSystemMatrix( *stdMat,
                                           dirichletEQN_ + offset_[i],
                                           aux,
                                           penaltyTerm_ );
          }
        }
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


    // --------------
    // StdMatrix case
    // --------------
    if ( sbmCase_ == false ) {
        SingleVector* stdVec = dynamic_cast<SingleVector*>(rhs);
        SingleVector* stdVal = dynamic_cast<SingleVector*>(dirichletValue_);
        if(stdVec == NULL || stdVal == NULL) EXCEPTION("Invalid cast attempt!");
        assembler_->AdaptRHSForIDBC( *stdVec,
            dirichletEQN_,
            *stdVal,
            penaltyTerm_,
            numIDBC_ );
    }

    // ---------------
    // SBM_Matrix case
    // ---------------
    else {
        // First down-cast vectors from Base to SBM
      SBM_Vector* sbmVec = dynamic_cast<SBM_Vector*>(rhs);
      SBM_Vector* sbmVal = dynamic_cast<SBM_Vector*>(dirichletValue_);
      if(sbmVec == NULL || sbmVal == NULL) EXCEPTION("Invalid cast attempt!");

        // We need pointers to sub-vectors
        SingleVector *stdVec = NULL;
        SingleVector *stdVal = NULL;

        // Now loop over all PDEs / sub-vectors and
        // set the Dirichlet values
        UInt aux = 0;
        for ( UInt i = 0; i < sbmVec->GetSize() ; i++ ) {

          // obtain sub-vectors
          stdVec = sbmVec->GetPointer( i );
          stdVal = sbmVal->GetPointer( i );

          // get number of IDBCs for i-th PDE
          aux = offset_[i+1] - offset_[i];

          // Now adapt the sub-vector
          if ( aux > 0 && stdVec != NULL && stdVal != NULL ) {
            assembler_->AdaptRHSForIDBC( *stdVec,
                                         dirichletEQN_ + offset_[i],
                                         *stdVal,
                                         penaltyTerm_,
                                         aux );
          }
        }
    }// else
  }


  // ***********
  //   SetIDBC
  // ***********
  template <typename T>
  void IDBC_HandlerPenalty<T>::SetIDBC( PdeIdType pdeID, UInt eqnNo,
                                        const T &val ) {


    // StdMatrix case
    UInt index = 0;
    if ( sbmCase_ == false ) {
      if( bcIndices_[pdeID].find( eqnNo ) ==
          bcIndices_[pdeID].end() ) {
        index = nextIndex_;
        bcIndices_[pdeID][eqnNo] = nextIndex_++;
      } else {
        index = bcIndices_[pdeID][eqnNo];
      }

      dirichletValue_->SetEntry( index, val );
      dirichletEQN_      [index] = eqnNo;
    }

    // SBM_Matrix case
    else {
      if( bcIndices_[pdeID].find( eqnNo ) ==
               bcIndices_[pdeID].end() ) {
        index = bcIndices_[pdeID].size()+1;
        bcIndices_[pdeID][eqnNo] = index;
      } else  {
        index = bcIndices_[pdeID][eqnNo];
      }

      SBM_Vector *sbmVec = dynamic_cast<SBM_Vector*>( dirichletValue_ );
      SingleVector *stdVec = sbmVec->GetPointer( pdeID );
      stdVec->SetEntry( index, val );
      dirichletEQN_      [ index+offset_[pdeID] ] = eqnNo;
    }
  }


  // *****************
  //   SetDofsToIDBC
  // *****************
  template <typename T>
  void IDBC_HandlerPenalty<T>::SetDofsToIDBC( BaseVector *vec ) {


    // -----------------
    // CASE 1: StdMatrix
    // -----------------
    if ( sbmCase_ == false ) {

      // Down-cast vector
      SingleVector *stdVec = dynamic_cast<SingleVector*>( vec );
      if ( stdVec == NULL ) {
        EXCEPTION( WRONG_CAST_MSG );
      }

      // Insert values
      T aux;
      for ( UInt i = 0; i < numIDBC_; i++ ) {
        dirichletValue_->GetEntry( i, aux );
        assembler_->SetVectorEntry( stdVec, dirichletEQN_[i],
                                    aux );
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
        EXCEPTION( WRONG_CAST_MSG );
      }

      // Loop over all sub-vectors
      SingleVector *stdVec = NULL;
      SingleVector *stdVal = NULL;
      T aux;
      UInt idx;

      for ( UInt i = 0; i < sbmVec->GetSize(); i++ ) {
/*
        (*debug) << " IDBC_HandlerPenalty::SetDofsToIDBC:"
                 << " Setting " << offset_[i+1] - offset_[i]
                 << " IDBCs for PDE #" << i << std::endl;
*/
        stdVec = sbmVec->GetPointer( i );
        stdVal = sbmVal->GetPointer( i );

        // Insert values
        for ( UInt j = 0; j < offset_[i+1] - offset_[i]; j++ ) {
          stdVal->GetEntry( j, aux );
          idx = offset_[i] + j;
          assembler_->SetVectorEntry( stdVec, dirichletEQN_[idx],
                                       aux );
        }
      }
    }
  }

  // **********************
  //   PrintInternalState
  // **********************
  template <typename T>
  void IDBC_HandlerPenalty<T>::PrintInternalState( std::ostream &os ) {


    // ----------------
    // CASE1: StdSystem
    // ----------------
    if ( sbmCase_ == false ) {

      os << " IDBC_HandlerPenalty --- Internal State:\n"
         << " --------------------------------------\n"
         << " numIDBC_    = " << numIDBC_ << std::endl
         << " # | eqnNo | comp | value\n";

      T val;
      SingleVector* stdVal = dynamic_cast<SingleVector*>(dirichletValue_);
      if(stdVal == NULL) EXCEPTION("Invalid cast attempt!");

      for ( UInt i = 0; i < numIDBC_; i++ ) {
        stdVal->GetEntry( i, val );
        os << i << " | "
           << dirichletEQN_[i] << " | "
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

      SingleVector *stdVal = NULL;
      SBM_Vector* sbmVec = dynamic_cast<SBM_Vector*>(dirichletValue_);
      if(sbmVec == NULL) EXCEPTION("Invalid cast attempt!");
      UInt aux = 0;
      UInt idx = 0;
      T val;

      // Compute column widths for pretty-printing
      aux = numIDBC_ + 1;
      UInt cw1 = aux > 0 ? (UInt)std::log10( (float)aux ) + 1 : 1;
      UInt cw2 = 5;
      // UInt cw3 = 4; // TODO: Check if this is still needed

      for ( UInt i = 0; i < dirichletValue_->GetSize(); i++ ) {

        // Compute number of IDBCs
        aux = offset_[i+1] - offset_[i];

        // Get hold of sub-vector
        stdVal = sbmVec->GetPointer(i);

        // Print IDBCs for this PDE
        os << " PDE # " << i << " / numIDBCs = " << aux << std::endl;
        if ( aux > 0 ) {
          os << " # | eqnNo | comp | value" << std::endl;
          for ( UInt j = 0; j < aux; j++ ) {

            stdVal->GetEntry( j, val );
            idx = offset_[i] + j;

            os << std::setfill( ' ' )
               << std::setw(cw1) << idx << " | "
               << std::setw(cw2) << dirichletEQN_[idx] << " | "
               << val << std::endl;
          }
        }
        os << std::endl;
      }
    }

    os << " --------------------------------------\n" << std::endl;

  }

  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class IDBC_HandlerPenalty<Double>;
    template class IDBC_HandlerPenalty<Complex>;
  #endif  
}
