// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "algsys/idbchandler.hh"
#include "graph/graph.hh"
#include "matvec/matvec.hh"


namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  template <typename T>
  IDBC_Handler<T>::IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                                 BaseGraphManager *graphManager,
                                 UInt numPDEs, bool sbmCase ) {

    ENTER_FCN( "IDBC_Handler::IDBC_Handler" );

    // -----------------------
    // Obtain basic dimensions
    // -----------------------

    // Determine number of dofs fixed by inhomogeneous Dirichlet BCs
    // and number of free dofs
    if ( sbmCase == false ) {
      NewArray( numFixedDofs_, UInt, 1 );
      NewArray( numFreeDofs_ , UInt, 1 );
      numFreeDofs_[1]  = graphManager->GetIDBCGraph()->GetSize();
      numFixedDofs_[1] = graphManager->GetIDBCGraph()->GetNumColsMat();
    }
    else {

      // Down-cast graph manager
      GraphManagerSBMMat *gmSBM = NULL;
      gmSBM = dynamic_cast<GraphManagerSBMMat *>(graphManager);
      if ( gmSBM == NULL ) {
        Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
      }

      NewArray( numFixedDofs_, UInt, numPDEs );
      NewArray( numFreeDofs_ , UInt, numPDEs );

      // Determine number of free dofs per PDE
      // and initialise number of fixed dofs to zero
      for ( UInt i = 1; i <= numPDEs; i++ ) {
        numFreeDofs_[i]  = gmSBM->GetGraph(i)->GetSize();
        numFixedDofs_[i] = 0;
      }

      // Set number of fixed dofs per PDE
      for ( UInt i = 1; i <= numPDEs; i++ ) {
        for ( UInt j = 1; j <= numPDEs; j++ ) {
          if ( gmSBM->IDBCGraphExists(i,j) == true ) {
            numFixedDofs_[i] = gmSBM->GetIDBCGraph(i,j)->GetNumColsMat();
          }
        }
      }
    }

    // -----------------------------------
    // Process set of required FE matrices
    // -----------------------------------

    // Check size of set
    if ( usedFEMatrices.empty() == true ) {
      (*error) << "IDBC_Handler::IDBC_Handler: Call to constructor with "
               << "empty usedFEMatrices set!";
      Error( __FILE__, __LINE__ );
    }
    else if ( usedFEMatrices.size() > MAX_NUM_FE_MATRICES ) {
      (*error) << "IDBC_Handler::IDBC_Handler: Detected serious "
               << "inconsistency! Size of usedFEMatrices set exceeds "
               << "value of MAX_NUM_FE_MATRICES macro";
      Error( __FILE__, __LINE__ );
    }

    // Since the set will only have a few entries we can afford to
    // generate a copy
    std::set<FEMatrixType>::const_iterator it;
    for ( it = usedFEMatrices.begin(); it != usedFEMatrices.end(); it++ ) {
      myMatrices_.insert( *it );
    }

    // ----------------------------------------------------
    // Generate empty matrices for storing IDBC information
    // ----------------------------------------------------

    // Determine entry type from template parameter
    MatrixEntryType eType = NOENTRYTYPE;
    if ( assocType<T>::tagM == assocType<Double>::tagM ) {
      eType = DOUBLE;
    }
    else if ( assocType<T>::tagM == assocType<Complex>::tagM ) {
      eType = COMPLEX;
    }
    else {
      (*error) << "Internal template error! No swearing please!";
      Error( __FILE__, __LINE__ );
    }

    // Initialise basic pointer structure
    NewArray( auxMat_, BaseMatrix*, MAX_NUM_FE_MATRICES );
    for ( UInt i = 1; i <= MAX_NUM_FE_MATRICES; i++ ) {
      auxMat_[i] = NULL;
    }

    // STDMATRIX case:
    if ( sbmCase == false ) {

      // Note: We always use CRS_Matrices currently. The entry type depends
      //       on the class' template parameter and the blocksize is
      //       restricted to the scalar setting.

      // Get hold of matrix graph
      BaseGraph *graph = graphManager->GetIDBCGraph();

      // Generation of required StdMatrices
      StdMatrix *stdMat = NULL;
      for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {

        // First generate empty object
        stdMat = GenerateStdMatrixObject( eType,
                                          SPARSE_NONSYM,
                                          1,
                                          numFreeDofs_[1],
                                          numFixedDofs_[1],
                                          graph->GetNNE() );

        // Now set sparsity pattern
        stdMat->SetSparsityPattern( *graph );

        // Insert pointer into BaseMatrix array
        auxMat_[*it] = stdMat;
      }
    }

    // SBMMATRIX case:
    else {

      SBM_Matrix *sbmMat = NULL;
      BaseGraph *graph = NULL;

      // Down-cast graph manager
      GraphManagerSBMMat *gmSBM = NULL;
      gmSBM = dynamic_cast<GraphManagerSBMMat *>(graphManager);
      if ( gmSBM == NULL ) {
        Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
      }

      // Generation of required SBM_Matrices
      for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {

        // Generate empty SBM_Matrix
        sbmMat = New SBM_Matrix( numPDEs, numPDEs, false );
        if ( sbmMat == NULL ) {
          (*error) << "IDBC_Handler::IDBC_Handler: "
                   << "Failed to generate SBM companion matrix for "
                   << "FE-Matrix '" << Enum2String( *it ) << "'";
          Error( __FILE__, __LINE__ );
        }

        // Populate matrix with sub-matrices for storing free-fixed-dof
        // connectivity
        for ( UInt i = 1; i <= numPDEs; i++ ) {
          for ( UInt k = 1; k <= numPDEs; k++ ) {

            if ( gmSBM->IDBCGraphExists(i,k) == true ) {

              graph = gmSBM->GetIDBCGraph(i,k);

              sbmMat->SetSubMatrix( i, k, eType, SPARSE_NONSYM, 1,
                                    graph->GetSize(),
                                    graph->GetNumColsMat(),
                                    graph->GetNNE() );

              sbmMat->GetPointer(i,k)->SetSparsityPattern( *graph );
            }
          }
        }

        // Insert pointer into BaseMatrix array
        auxMat_[*it] = sbmMat;
      }
    }


    // ------------------------------------------------
    // Generate vector for storing the Dirichlet values
    // ------------------------------------------------
    vecIDBC_ = GenerateVectorObject( *auxMat_[*(myMatrices_.begin())] );
    vecIDBC_->Init();


    // -------------------------
    // Set internal status flags
    // -------------------------
    addIDBCPossible_ = false;
    remIDBCPossible_ = false;
    sbmCase_         = sbmCase;


    // --------------------
    // Log some information
    // --------------------
    if ( sbmCase == false ) {
      (*cla) << "\n IDBC_Handler: Administrating "
             << numFixedDofs_[1] << " inhom. Dirichlet BCs.\n\n";
    }
    else {
    }

  }


  // **************
  //   Destructor
  // **************
  template <typename T> IDBC_Handler<T>::~IDBC_Handler() {

    ENTER_FCN( "IDBC_Handler::~IDBC_Handler" );

    // Delete vector of Dirichlet values
    delete vecIDBC_;

    // Delete equation number splitting arrays
    DeleteArray( numFreeDofs_  );
    DeleteArray( numFixedDofs_ );

    // Delete internal FE matrices
    mySetIterator it;
    for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {
      delete auxMat_[*it];
    }
    DeleteArray( auxMat_ );
  }


  // *********************
  //   BuiltSystemMatrix
  // *********************
  template <typename T> void IDBC_Handler<T>::
  BuiltSystemMatrix( const std::map<FEMatrixType, Double> &factors ) {

    ENTER_FCN( "IDBC_Handler::BuiltSystemMatrix" );

    BaseMatrix *sys = auxMat_[SYSTEM];

    std::map<FEMatrixType,Double>::const_iterator it;

    for ( it = factors.begin(); it != factors.end(); it++ ) {

      // Check that we have a factor and a FE matrix
      if ( auxMat_[(*it).first] != NULL  && (*it).second != 0.0 ) {
        sys->Add( (*it).second, *auxMat_[(*it).first] );
      }
    }

    // Adapt internal status flags
    addIDBCPossible_ = true;
  }


  // ***********
  //   SetIDBC
  // ***********
  template <typename T>
  void IDBC_Handler<T>::SetIDBC( PdeIdType pdeID, UInt eqnNo, UInt comp,
                                 UInt bcNum, const T &val ) {

    ENTER_FCN( "IDBC_Handler::SetIDBC" );

    // CASE 1: StdVector
    if ( sbmCase_ == false ) {
      vecIDBC_->SetVectorEntry( eqnNo - numFreeDofs_[1], val );
    }

    // CASE 2: SBM_Vector
    else {
      PtrCast( vecIDBC_, SBM_Vector, sbmVec );
      sbmVec->GetPointer(pdeID)->SetVectorEntry( eqnNo - numFreeDofs_[pdeID],
                                                 val );
    }

    // Adapt status flags
    remIDBCPossible_ = false;
  }

 // ************************
  //   AddWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::AddWeightFixedToFree( FEMatrixType matID,
                                              PdeIdType pdeID1,
                                              PdeIdType pdeID2,
                                              UInt rowInd,
                                              UInt colInd,
                                              Double realPart,
                                              Double imagPart ) {

    ENTER_FCN( "IDBC_Handler::AddWeightFixedToFree" );

    // ---------------------------
    // Obtain pointer to StdMatrix
    // ---------------------------
    StdMatrix *stdMat = NULL;

    // CASE 1: StdMatrix
    if ( sbmCase_ == false ) {
      stdMat = dynamic_cast<StdMatrix*>( auxMat_[matID] );
    }

    // CASE 2: SBM_Matrix
    else {

      SBM_Matrix *sbmMat = dynamic_cast<SBM_Matrix*>( auxMat_[matID] );
      stdMat = sbmMat->GetPointer( pdeID1, pdeID2 );

#ifdef DEBUG_IDBCHANDLER
      if ( stdMat == NULL ) {
        (*error) << "IDBC_Handler::SetWeightFixedToFree: Invalid access "
                 << "to non-existant sub-matrix (" << pdeID1
                 << ", " << pdeID2;
        Error( __FILE__, __LINE__ );
      }
#endif
    }

    // -----------------
    // Now add the value
    // -----------------

    // In the real-valued case find correct sub-matrix and insert value
    if ( assocType<T>::tagM == assocType<Double>::tagM ) {
      stdMat->AddToMatrixEntry( rowInd, colInd, realPart );
    }

    // We can savely assume that we are in the complex case now, since
    // otherwise an error would already have occured in the constructor.
    // Here we do the same as above, but must assemble the complex value
    // first
    else {
      Complex aux(realPart, imagPart);
      stdMat->AddToMatrixEntry( rowInd, colInd, aux );
    }
  }


  // ************************
  //   SetWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::SetWeightFixedToFree( FEMatrixType matID,
                                              PdeIdType pdeID1,
                                              PdeIdType pdeID2,
                                              UInt rowInd,
                                              UInt colInd,
                                              Double realPart,
                                              Double imagPart ) {

    ENTER_FCN( "IDBC_Handler::SetWeightFixedToFree" );

    // ---------------------------
    // Obtain pointer to StdMatrix
    // ---------------------------
    StdMatrix *stdMat = NULL;

    // CASE 1: StdMatrix
    if ( sbmCase_ == false ) {
      stdMat = dynamic_cast<StdMatrix*>( auxMat_[matID] );
    }

    // CASE 2: SBM_Matrix
    else {

      SBM_Matrix *sbmMat = dynamic_cast<SBM_Matrix*>( auxMat_[matID] );
      stdMat = sbmMat->GetPointer( pdeID1, pdeID2 );

#ifdef DEBUG_IDBCHANDLER
      if ( stdMat == NULL ) {
        (*error) << "IDBC_Handler::SetWeightFixedToFree: Invalid access "
                 << "to non-existant sub-matrix (" << pdeID1
                 << ", " << pdeID2;
        Error( __FILE__, __LINE__ );
      }
#endif
    }

    // -----------------
    // Now set the value
    // -----------------

    // In the real-valued case find correct sub-matrix and insert value
    if ( assocType<T>::tagM == assocType<Double>::tagM ) {
      stdMat->SetMatrixEntry( rowInd, colInd, realPart, false );
    }

    // We can savely assume that we are in the complex case now, since
    // otherwise an error would already have occured in the constructor.
    // Here we do the same as above, but must assemble the complex value
    // first
    else {
      Complex aux(realPart, imagPart);
      stdMat->SetMatrixEntry( rowInd, colInd, aux, false );
    }
  }

  // ************************
  //   GetWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::GetWeightFixedToFree( FEMatrixType matID, 
                                              PdeIdType pdeID1,
                                              PdeIdType pdeID2, 
                                              UInt rowInd, UInt colInd,
                                              Double & realPart,
                                              Double & imagPart ) const {

    ENTER_FCN( "IDBC_Handler::GetWeightFixedToFree" );

      // ---------------------------
    // Obtain pointer to StdMatrix
    // ---------------------------
    StdMatrix *stdMat = NULL;

    // CASE 1: StdMatrix
    if ( sbmCase_ == false ) {
      stdMat = dynamic_cast<StdMatrix*>( auxMat_[matID] );
    }

    // CASE 2: SBM_Matrix
    else {

      SBM_Matrix *sbmMat = dynamic_cast<SBM_Matrix*>( auxMat_[matID] );
      stdMat = sbmMat->GetPointer( pdeID1, pdeID2 );

#ifdef DEBUG_IDBCHANDLER
      if ( stdMat == NULL ) {
        (*error) << "IDBC_Handler::GetWeightFixedToFree: Invalid access "
                 << "to non-existant sub-matrix (" << pdeID1
                 << ", " << pdeID2;
        Error( __FILE__, __LINE__ );
      }
#endif
    }

    // -----------------
    // Now get the value
    // -----------------

    // In the real-valued case find correct sub-matrix and get value
    if ( assocType<T>::tagM == assocType<Double>::tagM ) {
      stdMat->GetMatrixEntry( rowInd, colInd, realPart );
    }

    // We can savely assume that we are in the complex case now, since
    // otherwise an error would already have occured in the constructor.
    // Here we do the same as above, but must disassemble the complex value
    // afterwards
    else {
      Complex aux;
      stdMat->GetMatrixEntry( rowInd, colInd, aux );
      realPart = aux.real();
      imagPart = aux.imag();
    }
  }
    

  // *****************
  //   SetRowWeights
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetRowWeights( FEMatrixType matID, 
                                       PdeIdType pdeID, UInt rowInd,
                                       Double realPart, 
                                       Double imagPart ) {

    ENTER_FCN( "IDBC_Handler::SetRowWeights" );

  }


  // *****************
  //   SetColWeights
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetColWeights( FEMatrixType matID, 
                                       PdeIdType pdeID,UInt colInd,
                                       Double realPart, 
                                       Double imagPart ) {

    ENTER_FCN( "IDBC_Handler::SetColWeights" );
  }
  
  // *****************
  //   SetDofsToIDBC
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetDofsToIDBC( BaseVector *vec ) {

    ENTER_FCN( "IDBC_Handler::SetDofsToIDBC" );

    // ------------------
    // CASE 1: SBM_Vector
    // ------------------
    if ( sbmCase_ == true ) {

      // Down-cast vector
      SBM_Vector *sbmVec = dynamic_cast<SBM_Vector*>( vec );
      SBM_Vector *sbmVal = dynamic_cast<SBM_Vector*>( vecIDBC_ );
      if ( sbmVec == NULL || sbmVal == NULL ) {
        Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
      }

      // Loop over all sub-vectors
      StdVector *stdVec = NULL;
      StdVector *stdVal = NULL;
      for ( UInt i = 1; i <= sbmVec->GetSize(); i++ ) {

        stdVec = sbmVec->GetPointer( i );
        stdVal = sbmVal->GetPointer( i );

        // Insert values
        for ( UInt j = 1; j <= numFixedDofs_[i]; j++ ) {
          T aux;
          stdVal->GetVectorEntry( j, aux );
          stdVec->SetVectorEntry( j + numFreeDofs_[i], aux );
        }
      }
    }

    // -----------------
    // CASE 2: StdVector
    // -----------------
    else {

      // Down-cast vector
      StdVector *stdVec = dynamic_cast<StdVector*>( vec );
      if ( stdVec == NULL ) {
        Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
      }

      // Insert values
      for ( UInt i = 1; i <= numFixedDofs_[1]; i++ ) {
        T aux;
        vecIDBC_->GetVectorEntry( i, aux );
        vec->SetVectorEntry( i + numFreeDofs_[1], aux );
      }
    }
  }


  // ****************************
  //   InstantiatePublicMethods
  // ****************************
  template <typename T>
  void IDBC_Handler<T>::InstantiatePublicMethods() {

    ENTER_FCN( "IDBC_Handler::InstantiatePublicMethods" );

    Error( "This function should never be called", __FILE__, __LINE__ );

    // Some dummy variables we need
    BaseVector *dummyVec = NULL;
    std::map<FEMatrixType,Double> dummyMap;
    BaseMatrix *dummyMat = NULL;

    // Public methods in alphabetical order
    AdaptSystemMatrix( *dummyMat );
    AddIDBCToRHS( dummyVec );
    BuiltSystemMatrix( dummyMap );
    InitDirichletValues();
    InitMatrix( SYSTEM );
    RemoveIDBCFromRHS( dummyVec );
    SetDofsToIDBC( dummyVec );
    SetIDBC( NO_PDE_ID, 1, 0, 0, 0.0 );
    SetWeightFixedToFree( SYSTEM, 1, 2, 3, 4, 1.0, 2.0 );
  }

}
