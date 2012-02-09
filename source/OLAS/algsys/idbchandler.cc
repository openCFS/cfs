// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <stddef.h>
#include <string>
#include <utility>

#include "MatVec/SingleVector.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/generatematvec.hh"
#include "MatVec/sbmmatrix.hh"
#include "MatVec/sbmvector.hh"
#include "MatVec/stdmatrix.hh"
#include "MatVec/typedefs.hh"
#include "OLAS/algsys/idbchandler.hh"
#include "OLAS/graph/basegraph.hh"
#include "OLAS/graph/basegraphmanager.hh"
#include "OLAS/graph/graphmanagersbmmat.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template <typename T>
  IDBC_Handler<T>::IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                                 BaseGraphManager *graphManager,
                                 UInt numPDEs, bool sbmCase ) {


    // -----------------------
    // Obtain basic dimensions
    // -----------------------

    // Determine number of dofs fixed by inhomogeneous Dirichlet BCs
    // and number of free dofs
    if ( sbmCase == false ) {
      NEWARRAY( numFixedDofs_, UInt, 1 );
      NEWARRAY( numFreeDofs_ , UInt, 1 );
      numFreeDofs_[0]  = graphManager->GetIDBCGraph()->GetSize();
      numFixedDofs_[0] = graphManager->GetIDBCGraph()->GetNumColsMat();
    }
    else {

      // Down-cast graph manager
      GraphManagerSBMMat *gmSBM = NULL;
      gmSBM = dynamic_cast<GraphManagerSBMMat *>(graphManager);
      if ( gmSBM == NULL ) {
        EXCEPTION( WRONG_CAST_MSG );
      }

      NEWARRAY( numFixedDofs_, UInt, numPDEs );
      NEWARRAY( numFreeDofs_ , UInt, numPDEs );

      // Determine number of free dofs per PDE
      // and initialise number of fixed dofs to zero
      for ( UInt i = 0; i < numPDEs; i++ ) {
        numFreeDofs_[i]  = gmSBM->GetGraph(i)->GetSize();
        numFixedDofs_[i] = 0;
      }

      // Set number of fixed dofs per PDE
      for ( UInt i = 0; i < numPDEs; i++ ) {
        for ( UInt j = 0; j < numPDEs; j++ ) {
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
      EXCEPTION( "IDBC_Handler::IDBC_Handler: Call to constructor with "
               << "empty usedFEMatrices set!");
    }
    else if ( usedFEMatrices.size() > MAX_NUM_FE_MATRICES ) {
      EXCEPTION( "IDBC_Handler::IDBC_Handler: Detected serious "
               << "inconsistency! Size of usedFEMatrices set exceeds "
               << "value of MAX_NUM_FE_MATRICES macro");
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
    BaseMatrix::EntryType eType = BaseMatrix::NOENTRYTYPE;
    if ( AssocType<T>::tagM == AssocType<Double>::tagM ) {
      eType = BaseMatrix::DOUBLE;
    }
    else if ( AssocType<T>::tagM == AssocType<Complex>::tagM ) {
      eType = BaseMatrix::COMPLEX;
    }
    else {
      EXCEPTION( "Internal template error! No swearing please!" );
    }

    // Initialise basic pointer structure
    NEWARRAY( auxMat_, BaseMatrix*, MAX_NUM_FE_MATRICES );
    for ( UInt i = 0; i < MAX_NUM_FE_MATRICES; i++ ) {
      auxMat_[i] = NULL;
    }

    // SPARSE_MATRIX case:
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
                                          BaseMatrix::SPARSE_NONSYM,
                                          numFreeDofs_[0],
                                          numFixedDofs_[0],
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
        EXCEPTION( WRONG_CAST_MSG );
      }

      // Generation of required SBM_Matrices
      for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {

        // Generate empty SBM_Matrix
        sbmMat = new SBM_Matrix( numPDEs, numPDEs, false );
        if ( sbmMat == NULL ) {
          std::string tmp;
          Enum2String( *it, tmp );
          EXCEPTION( "IDBC_Handler::IDBC_Handler: "
                   << "Failed to generate SBM companion matrix for "
                   << "FE-Matrix '" << tmp << "'");
        }

        // Populate matrix with sub-matrices for storing free-fixed-dof
        // connectivity
        for ( UInt i = 0; i < numPDEs; i++ ) {
          for ( UInt k = 0; k < numPDEs; k++ ) {

            if ( gmSBM->IDBCGraphExists(i,k) == true ) {

              graph = gmSBM->GetIDBCGraph(i,k);

              sbmMat->SetSubMatrix( i, k, eType, BaseMatrix::SPARSE_NONSYM,
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
             << numFixedDofs_[0] << " inhom. Dirichlet BCs.\n\n";
    }
    else {
    }

  }


  // **************
  //   Destructor
  // **************
  template <typename T> IDBC_Handler<T>::~IDBC_Handler() {


    // Delete vector of Dirichlet values
    delete vecIDBC_;

    // Delete equation number splitting arrays
    delete [] ( numFreeDofs_  );
    delete [] ( numFixedDofs_ );

    // Delete internal FE matrices
    mySetIterator it;
    for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {
      delete auxMat_[*it];
    }
    delete [] ( auxMat_ );
  }


  // *********************
  //   BuiltSystemMatrix
  // *********************
  template <typename T> void IDBC_Handler<T>::
  BuiltSystemMatrix( const std::map<FEMatrixType, Double> &factors ) {


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
  void IDBC_Handler<T>::SetIDBC( PdeIdType pdeID, UInt eqnNo, const T &val ) {


    // CASE 1: SingleVector
    if ( sbmCase_ == false ) {
      vecIDBC_->SetEntry( eqnNo - numFreeDofs_[0] - 1, val );
    }

    // CASE 2: SBM_Vector
    else {
      SBM_Vector* sbmVec = dynamic_cast<SBM_Vector*>(vecIDBC_);
      if(sbmVec == NULL) EXCEPTION("Invalid cast attempt!");
      sbmVec->GetPointer(pdeID)->SetEntry( eqnNo - numFreeDofs_[pdeID] - 1,
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
                                              const T& val ) {


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
        EXCEPTION("IDBC_Handler::SetWeightFixedToFree: Invalid access "
                 << "to non-existant sub-matrix (" << pdeID1
                 << ", " << pdeID2);
      }
#endif
    }

    // -----------------
    // Now add the value
    // -----------------
    stdMat->AddToMatrixEntry( rowInd, colInd, val );
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
                                              const T& val ) {


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
        EXCEPTION("IDBC_Handler::SetWeightFixedToFree: Invalid access "
                 << "to non-existant sub-matrix (" << pdeID1
                 << ", " << pdeID2);
      }
#endif
    }

    // -----------------
    // Now set the value
    // -----------------
    stdMat->SetMatrixEntry( rowInd, colInd, val, false );
  }

  // ************************
  //   GetWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::GetWeightFixedToFree( FEMatrixType matID, 
                                              PdeIdType pdeID1,
                                              PdeIdType pdeID2, 
                                              UInt rowInd, UInt colInd,
                                              T& val ) const {


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
        EXCEPTION("IDBC_Handler::GetWeightFixedToFree: Invalid access "
                 << "to non-existant sub-matrix (" << pdeID1
                 << ", " << pdeID2);
      }
#endif
    }

    // -----------------
    // Now get the value
    // -----------------
    stdMat->GetMatrixEntry( rowInd, colInd, val );
    
  }
    

  // *****************
  //   SetRowWeights
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetRowWeights( FEMatrixType matID, 
                                       PdeIdType pdeID, UInt rowInd,
                                       const T& val ) {


  }


  // *****************
  //   SetColWeights
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetColWeights( FEMatrixType matID, 
                                       PdeIdType pdeID,UInt colInd,
                                       const T& val ) {

  }
  
  // *****************
  //   SetDofsToIDBC
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetDofsToIDBC( BaseVector *vec ) {


    // ------------------
    // CASE 1: SBM_Vector
    // ------------------
    if ( sbmCase_ == true ) {

      // Down-cast vector
      SBM_Vector *sbmVec = dynamic_cast<SBM_Vector*>( vec );
      SBM_Vector *sbmVal = dynamic_cast<SBM_Vector*>( vecIDBC_ );
      if ( sbmVec == NULL || sbmVal == NULL ) {
        EXCEPTION( WRONG_CAST_MSG );
      }

      // Loop over all sub-vectors
      SingleVector *stdVec = NULL;
      SingleVector *stdVal = NULL;
      for ( UInt i = 0; i < sbmVec->GetSize(); i++ ) {

        stdVec = sbmVec->GetPointer( i );
        stdVal = sbmVal->GetPointer( i );

        // Insert values
        for ( UInt j = 0; j < numFixedDofs_[i]; j++ ) {
          T aux;
          stdVal->GetEntry( j, aux );
          stdVec->SetEntry( j + numFreeDofs_[i], aux );
        }
      }
    }

    // -----------------
    // CASE 2: SingleVector
    // -----------------
    else {

      // Down-cast vector
      SingleVector *stdVec = dynamic_cast<SingleVector*>( vec );
      if ( stdVec == NULL ) {
        EXCEPTION( WRONG_CAST_MSG );
      }

      // Insert values
      for ( UInt i = 0; i < numFixedDofs_[0]; i++ ) {
        T aux;
        vecIDBC_->GetEntry( i, aux );
        vec->SetEntry( i + numFreeDofs_[0], aux );
      }
    }
  }

  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class IDBC_Handler<Double>;
    template class IDBC_Handler<Complex>;
  #endif
  
}
