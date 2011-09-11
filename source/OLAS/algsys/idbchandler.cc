// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "MatVec/stdmatrix.hh"
#include "OLAS/graph/graphmanager.hh"
#include "OLAS/algsys/idbchandler.hh"

#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(idbcElim)
DEFINE_LOG(idbcElim, "idbcElim")

namespace CoupledField {
  


  // ***************
  //   Constructor
  // ***************
  template <typename T>
  IDBC_Handler<T>::IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                                 GraphManager *gM,
                                 UInt numBlocks ) {

    LOG_TRACE(idbcElim) << "Generating IDBC_Handler for " << numBlocks << "blocks";
    
    // -----------------------
    // Obtain basic dimensions
    // -----------------------

    // Determine number of dofs fixed by inhomogeneous Dirichlet BCs
    // and number of free dofs
    numFixedDofs_.Resize( numBlocks );
    numFreeDofs_.Resize( numBlocks );

    // Determine number of free dofs per block
    // and initialise number of fixed dofs to zero
    for ( UInt i = 0; i < numBlocks; i++ ) {
      numFreeDofs_[i]  = gM->GetGraph(i)->GetSize();
      numFixedDofs_[i] = 0;
    }

    // Set number of fixed dofs per block
    for ( UInt i = 0; i < numBlocks; i++ ) {
      for ( UInt j = 0; j < numBlocks; j++ ) {
        LOG_DBG(idbcElim) << "Block (" << i  << ", " << j << "):";
        LOG_DBG(idbcElim) << "\t#free dofs: " << numFreeDofs_[i];

        if ( gM->IDBCGraphExists(i,j) == true ) {
          numFixedDofs_[i] = gM->GetIDBCGraph(i,j)->GetNumColsMat();
          
          LOG_DBG(idbcElim) << "\t#fixed dofs: " << numFixedDofs_[i];
        } else {
          LOG_DBG(idbcElim) << "\t#fixed dofs: -";
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

    // SBMMATRIX case:
    SBM_Matrix *sbmMat = NULL;
    BaseGraph *graph = NULL;
    
    // Generation of required SBM_Matrices
    LOG_DBG(idbcElim) << "Generating Auxiliary matrices for:";
    for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {

      LOG_DBG(idbcElim) << "\t-" << feMatrixType.ToString(*it);
      
      // Generate empty SBM_Matrix
      sbmMat = new SBM_Matrix( numBlocks, numBlocks, false );
      if ( sbmMat == NULL ) {
        std::string tmp;
        Enum2String( *it, tmp );
        EXCEPTION( "IDBC_Handler::IDBC_Handler: "
            << "Failed to generate SBM companion matrix for "
            << "FE-Matrix '" << tmp << "'");
      }

      // Populate matrix with sub-matrices for storing free-fixed-dof
      // connectivity
      for ( UInt i = 0; i < numBlocks; i++ ) {
        for ( UInt k = 0; k < numBlocks; k++ ) {
          
          if ( gM->IDBCGraphExists(i,k) == true ) {
            graph = gM->GetIDBCGraph(i,k);

            sbmMat->SetSubMatrix( i, k, eType, BaseMatrix::SPARSE_NONSYM,
                                  graph->GetSize(),
                                  graph->GetNumColsMat(),
                                  graph->GetNNE() );
            LOG_DBG3(idbcElim) << "\t\tblock(" << i << ", " << k << "), size " 
                                << graph->GetSize() << " x " << graph->GetNumColsMat()
                                << " with " <<  graph->GetNNE()  << " entries";
            sbmMat->GetPointer(i,k)->SetSparsityPattern( *graph );
          }
        }
      }

      // Insert pointer into BaseMatrix array
      auxMat_[*it] = sbmMat;
    }


    // ------------------------------------------------
    // Generate vector for storing the Dirichlet values
    // ------------------------------------------------
    vecIDBC_ = dynamic_cast<SBM_Vector*>(
        GenerateVectorObject( *auxMat_[*(myMatrices_.begin())] ));
    vecIDBC_->Init();

    // -------------------------
    // Set internal status flags
    // -------------------------
    addIDBCPossible_ = false;
    remIDBCPossible_ = false;
  }


  // **************
  //   Destructor
  // **************
  template <typename T> IDBC_Handler<T>::~IDBC_Handler() {

    // Delete vector of Dirichlet values
    delete vecIDBC_;

    // Delete internal FE matrices
    mySetIterator it;
    for ( it = myMatrices_.begin(); it != myMatrices_.end(); it++ ) {
      delete auxMat_[*it];
    }
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
  // ****************
  //   AddIDBCToRHS
  // ****************
  //! @copydoc BaseIDBC_Handler::AddIDBCToRHS()
  template <typename T> void IDBC_Handler<T>::
  AddIDBCToRHS( SBM_Vector *rhs ) {
    
    if ( addIDBCPossible_ == false ) {
      EXCEPTION( "IDBCHandler::AddIDBCToRHS: Internal error! Refusing to "
          << "add Dirichlet BCs, since addIDBCPossible_ = false! "
          << "Did you call BuiltSystemMatrix()?");
    }
    
    auxMat_[SYSTEM]->MultSub( *vecIDBC_, *rhs );
    remIDBCPossible_ = true;
  }

  // **********************
  //   RemoveIDBCFromRHS
  // **********************
  template <typename T> void IDBC_Handler<T>::
  RemoveIDBCFromRHS( SBM_Vector *rhs ) {
    
    if ( remIDBCPossible_ == false ) {
      EXCEPTION( "IDBCHandler::RemoveIDBCFromRHS: Internal error! "
          << "Refusing to remove Dirichlet BCs, since "
          << "remIDBCPossible_ = false! "
          << "Either FE matrices or Dirichlet values have changed "
          << "since last call to AddIDBCToRHS()!");
    }
    auxMat_[SYSTEM]->MultAdd( *vecIDBC_, *rhs );
  }


  // ***********
  //   SetIDBC
  // ***********
  template <typename T>
  void IDBC_Handler<T>::SetIDBC( UInt blockNum, UInt index, const T &val ) {

    vecIDBC_->GetPointer(blockNum)->
        SetEntry( index - numFreeDofs_[blockNum] - 1, val);

    // Adapt status flags
    remIDBCPossible_ = false;
  }
  
  // ***********
  //   GetIDBC
  // ***********
  template <typename T>
  void IDBC_Handler<T>::GetIDBC( UInt blockNum, UInt index, T &val ) {

    vecIDBC_->GetPointer(blockNum)->
        GetEntry( index - numFreeDofs_[blockNum] - 1, val);
  }

 // ************************
  //   AddWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::AddWeightFixedToFree( FEMatrixType matID,
                                              UInt rowBlock,
                                              UInt colBlock,
                                              UInt rowInd,
                                              UInt colInd,
                                              const T& val ) {

    StdMatrix *stdMat = auxMat_[matID]->GetPointer(rowBlock, colBlock);
    stdMat->AddToMatrixEntry( rowInd, colInd, val );
  }


  // ************************
  //   SetWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::SetWeightFixedToFree( FEMatrixType matID,
                                              UInt rowBlock,
                                              UInt colBlock,
                                              UInt rowInd,
                                              UInt colInd,
                                              const T& val ) {

    StdMatrix *stdMat = auxMat_[matID]->GetPointer( rowBlock, colBlock );
    stdMat->SetMatrixEntry( rowInd, colInd, val, false );
  }

  // ************************
  //   GetWeightFixedToFree
  // ************************
  template <typename T>
  void IDBC_Handler<T>::GetWeightFixedToFree( FEMatrixType matID,
                                              UInt rowBlock,
                                              UInt colBlock,
                                              UInt rowInd,
                                              UInt colInd,
                                              T& val )  {
    
    StdMatrix *stdMat = auxMat_[matID]->GetPointer( rowBlock, colBlock );
    stdMat->GetMatrixEntry( rowInd, colInd, val );
  }
    
  // *****************
  //   SetDofsToIDBC
  // *****************
  template <typename T>
  void IDBC_Handler<T>::SetDofsToIDBC( SBM_Vector *vec ) {

    // Loop over all sub-vectors
    SingleVector *stdVec = NULL;
    SingleVector *stdVal = NULL;
    for ( UInt i = 0; i < vec->GetSize(); i++ ) {

      stdVec = vec->GetPointer( i );
      stdVal = vecIDBC_->GetPointer( i );

      // Insert values
      for ( UInt j = 0; j < numFixedDofs_[i]; j++ ) {
        T aux;
        stdVal->GetEntry( j, aux );
        stdVec->SetEntry( j + numFreeDofs_[i], aux );
      }
    }
   
  }
  
  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class IDBC_Handler<Double>;
    template class IDBC_Handler<Complex>;
  #endif
  
}
