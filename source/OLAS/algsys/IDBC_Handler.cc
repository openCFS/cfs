#include "MatVec/StdMatrix.hh"
#include "OLAS/graph/GraphManager.hh"
#include "OLAS/algsys/IDBC_Handler.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include <boost/type_traits/is_complex.hpp>

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
    BaseMatrix::EntryType eType = boost::is_complex<T>() ? BaseMatrix::COMPLEX : BaseMatrix::DOUBLE;

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
      
    } // loop over matrix types


    // ------------------------------------------------
    // Generate vector for storing the Dirichlet values
    // ------------------------------------------------
    vecIDBC_ = dynamic_cast<SBM_Vector*>(
        GenerateVectorObject( *auxMat_[*(myMatrices_.begin())] ));
    vecIDBC_->Init();

    vecOldIDBC_ = dynamic_cast<SBM_Vector*>(
        GenerateVectorObject( *auxMat_[*(myMatrices_.begin())] ));
    vecOldIDBC_->Init();

    // ------------------------------------------------------------------
    // Generate vector for storing temporary values for AddFixedToFreeRHS
    // ------------------------------------------------------------------
    auxVec_ = dynamic_cast<SBM_Vector*>(
        GenerateVectorObject( *auxMat_[*(myMatrices_.begin())] ));
    auxVec_->Init();
    
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
    vecIDBC_ = NULL;
    
    delete vecOldIDBC_;
    vecOldIDBC_ = NULL;

    delete auxVec_;
    auxVec_ = NULL;

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
  BuildSystemMatrix( const std::map<FEMatrixType, Double> &factors,
                     std::map<UInt, std::set<UInt> >& colInd ) {

    SBM_Matrix *sys = auxMat_[SYSTEM];

    std::map<FEMatrixType,Double>::const_iterator it;

    for ( it = factors.begin(); it != factors.end(); it++ ) {

      // Check that we have a factor and a FE matrix
      if ( auxMat_[(*it).first] != NULL  && (*it).second != 0.0 ) {

        // generate empty set of indices (= all rows)
        std::map<UInt, std::set<UInt> > rowInd;

        sys->Add( (*it).second , *auxMat_[(*it).first], 
                  rowInd, colInd );
      }
    }

    // Adapt internal status flags
    addIDBCPossible_ = true;
  }
  // ****************
  //   AddIDBCToRHS
  // ****************
  template <typename T> void IDBC_Handler<T>::
  AddIDBCToRHS( SBM_Vector *rhs, bool deltaIDBC ) {
    
    if ( addIDBCPossible_ == false ) {
      EXCEPTION( "IDBCHandler::AddIDBCToRHS: Internal error! Refusing to "
          << "add Dirichlet BCs, since addIDBCPossible_ = false! "
          << "Did you call BuildSystemMatrix()?");
    }
    
    auxMat_[SYSTEM]->MultSub( *vecIDBC_, *rhs );
    if(deltaIDBC){
      auxMat_[SYSTEM]->MultAdd( *vecOldIDBC_, *rhs );
    }
    remIDBCPossible_ = true;
  }

  template <typename T> void IDBC_Handler<T>::
  SetOldDirichletValues() {
    /*
     * Copy entries of vecIDBC to vecOldIDBC
     */
    vecOldIDBC_->Init();
    vecOldIDBC_->Add(dynamic_cast<const SBM_Vector&>(*vecIDBC_));
  }

  template <typename T> void IDBC_Handler<T>::
  ToString() {
    std::cout <<"Old IDBC Values: " << vecOldIDBC_->ToString() << std::endl;
    std::cout <<"New IDBC Values: " << vecIDBC_->ToString() << std::endl;
  }

  // **********************
  //   RemoveIDBCFromRHS
  // **********************
  template <typename T> void IDBC_Handler<T>::
  RemoveIDBCFromRHS( SBM_Vector *rhs, bool deltaIDBC ) {
    
    if ( remIDBCPossible_ == false ) {
      EXCEPTION( "IDBCHandler::RemoveIDBCFromRHS: Internal error! "
          << "Refusing to remove Dirichlet BCs, since "
          << "remIDBCPossible_ = false! "
          << "Either FE matrices or Dirichlet values have changed "
          << "since last call to AddIDBCToRHS()!");
    }
    auxMat_[SYSTEM]->MultAdd( *vecIDBC_, *rhs );
    if(deltaIDBC){
      auxMat_[SYSTEM]->MultSub( *vecOldIDBC_, *rhs );
    }
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
  void IDBC_Handler<T>::GetIDBC( UInt blockNum, UInt index, T &val, bool deltaIDBC ) {

    vecIDBC_->GetPointer(blockNum)->
        GetEntry( index - numFreeDofs_[blockNum] - 1, val);

    if(deltaIDBC){
      T aux;
      vecOldIDBC_->GetPointer(blockNum)->
          GetEntry( index - numFreeDofs_[blockNum] - 1, aux);
      val -= aux;
    }
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

    LOG_DBG2(idbcElim) << "AddWeightFixedToFree:";
    LOG_DBG2(idbcElim) << "\tmatID:    " << matID;
    LOG_DBG2(idbcElim) << "\trowBlock: " << rowBlock;
    LOG_DBG2(idbcElim) << "\tcolBlock: " << colBlock;
    LOG_DBG2(idbcElim) << "\trowInd:   " << rowInd;
    LOG_DBG2(idbcElim) << "\tcolInd:   " << colInd;

    StdMatrix *stdMat = auxMat_[matID]->GetPointer(rowBlock, colBlock);
    stdMat->AddToMatrixEntry( rowInd, colInd, val );
  }


    // ************************
    //   AddFixedToFreeRHS
    // ************************
    template <typename T>
    void IDBC_Handler<T>::AddFixedToFreeRHS( FEMatrixType matID,
                                                  UInt colBlock,
                                                  UInt colInd,
                                                  SBM_Vector *rhs,
                                                  const T& val ) {

      LOG_DBG2(idbcElim) << "AddFixedToFreeRHS:";
      LOG_DBG2(idbcElim) << "\tmatID:    " << matID;
      LOG_DBG2(idbcElim) << "\tcolBlock: " << colBlock;
      LOG_DBG2(idbcElim) << "\tcolInd:   " << colInd;
      LOG_DBG2(idbcElim) << "\value:     " << val;

      auxVec_->Init();
      SingleVector * curVec = auxVec_->GetPointer(colBlock);
      curVec->SetEntry(colInd - numFreeDofs_[colBlock] - 1,val);

      auxMat_[matID]->MultAdd(*auxVec_,*rhs);
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
  void IDBC_Handler<T>::SetDofsToIDBC( SBM_Vector *vec, bool deltaIDBC ) {

    // Loop over all sub-vectors
    SingleVector *stdVec = NULL;
    SingleVector *stdVal = NULL;

    if(!deltaIDBC){
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
    } else {
      SingleVector *stdVal2 = NULL;

      for ( UInt i = 0; i < vec->GetSize(); i++ ) {

        stdVec = vec->GetPointer( i );
        stdVal = vecIDBC_->GetPointer( i );
        stdVal2 = vecOldIDBC_->GetPointer( i );

        // Insert values
        for ( UInt j = 0; j < numFixedDofs_[i]; j++ ) {
          T aux, aux2;
          stdVal->GetEntry( j, aux );
          stdVal2->GetEntry( j, aux2 );
          stdVec->SetEntry( j + numFreeDofs_[i], aux-aux2 );
        }
      }
    }
   
  }
  
  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class IDBC_Handler<Double>;
    template class IDBC_Handler<Complex>;
  #endif
  
}
