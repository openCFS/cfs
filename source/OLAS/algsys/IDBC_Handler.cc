#include "MatVec/StdMatrix.hh"
#include "OLAS/graph/GraphManager.hh"
#include "OLAS/algsys/IDBC_Handler.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include <boost/type_traits/is_complex.hpp>

DEFINE_LOG(idbcElim, "idbcElim")

namespace CoupledField {
  


  // ***************
  //   Constructor
  // ***************
  template <typename T>
  IDBC_Handler<T>::IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                                 GraphManager *gM,
                                 UInt numBlocks ) {

    LOG_TRACE(idbcElim) << "Generating IDBC_Handler for " << numBlocks << " blocks";

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
          } // k loop
        } // i loop

      // Insert pointer into BaseMatrix array
      auxMat_[*it] = sbmMat;

      LOG_DBG3(idbcElim) << "auxMat_[*it]->ToString()"<< std::endl;
      LOG_DBG3(idbcElim) << auxMat_[*it]->ToString();
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


  // ***************
  //   Constructor for MultiHarmonic
  // ***************
  template <typename T>
  IDBC_Handler<T>::IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                                 GraphManager *gM,
                                 UInt numBlocks, UInt M, UInt a ) {

    LOG_TRACE(idbcElim) << "Generating IDBC_Handler";

    // -----------------------
    // Obtain basic dimensions
    // -----------------------

    //Numbers of Frequency used
    UInt N = 2*M+1 ;
    // N*N - (M*M +M) = Numbers of occupied entries in multiharmonic system
    UInt octupiedEntries = N*N - (M*M+M);

    // Determine number of dofs fixed by inhomogeneous Dirichlet BCs
    // and number of free dofs (for each harmonic)
    numFixedDofs_.Resize( octupiedEntries );
    numFreeDofs_.Resize( octupiedEntries, 0 );

    // Set number of fixed dofs per block
    //===================================
    // counter to add these in the list numFreeDofs/numFixedDofs for each occupied entry
    // Not shure if numFreeDofs_ is needed lateron.
    // Added transposed ones the same order as later on. Hopefully than the order will fit...

    UInt counter = 0;
    for ( UInt sbmRow = 0; sbmRow < N; ++sbmRow ) {
      for ( UInt sbmCol = sbmRow; sbmCol < sbmRow +a; ++sbmCol ) {
        if( sbmCol < N){
          LOG_DBG(idbcElim) << "(sbmRow,sbmCol)=(" << sbmRow  << ", " << sbmCol << "):";
          if ( gM->IDBCGraphExists(sbmRow,sbmCol) == true ) {

            numFreeDofs_[counter]  = gM->GetGraph(sbmRow,sbmCol)->GetSize();
            numFixedDofs_[counter] = gM->GetIDBCGraph(sbmRow,sbmCol)->GetNumColsMat();
            LOG_DBG(idbcElim) << "\t#fixed dofs: " << numFixedDofs_[counter];
            if(sbmRow !=sbmCol){
              //Transposed ones
              counter++;
              numFreeDofs_[counter]  = gM->GetGraph(sbmCol,sbmRow)->GetSize();
              numFixedDofs_[counter] = gM->GetIDBCGraph(sbmCol,sbmRow)->GetNumColsMat();
              LOG_DBG(idbcElim) << "\t#fixed dofs: " << numFixedDofs_[counter];
            }
          } else {
            LOG_DBG(idbcElim) << "\t#fixed dofs: -";
          }

          counter++;
//          std::cout << counter << std::endl;

        } //if
      } //sbmCol
    } //sbmRow

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

      // Generate empty SBM_Matrix for Multiharmonic
        sbmMat = new SBM_Matrix( N, N, false );
        if ( sbmMat == NULL ) {
          std::string tmp;
          Enum2String( *it, tmp );
          EXCEPTION( "IDBC_Handler::IDBC_Handler: "
              << "Failed to generate SBM companion matrix for "
              << "FE-Matrix '" << tmp << "'");
        }
        // Populate matrix with sub-matrices for storing free-fixed-dof
        // connectivity

        for ( UInt sbmRow = 0; sbmRow < N; ++sbmRow ) {
          for ( UInt sbmCol = sbmRow; sbmCol < sbmRow +a; ++sbmCol ) {
            if( sbmCol < N){
              LOG_DBG(idbcElim) << "(sbmRow,sbmCol)=(" << sbmRow  << ", " << sbmCol << "):";
              if ( gM->IDBCGraphExists(sbmRow,sbmCol) == true ) {
                graph = gM->GetIDBCGraph(sbmRow,sbmCol);

                sbmMat->SetSubMatrix( sbmRow, sbmCol, eType, BaseMatrix::SPARSE_NONSYM,
                                      graph->GetSize(),
                                      graph->GetNumColsMat(),
                                      graph->GetNNE() );
                LOG_DBG3(idbcElim) << "\t\t (sbmRow,sbmCol)=(" << sbmRow << ", " << sbmCol << "), size "
                                    << graph->GetSize() << " x " << graph->GetNumColsMat()
                                    << " with " <<  graph->GetNNE()  << " entries";
                sbmMat->GetPointer(sbmRow,sbmCol)->SetSparsityPattern( *graph );
                // Do not forget on the tranposed entries
                if (sbmRow !=sbmCol){
                  if ( gM->IDBCGraphExists(sbmCol,sbmRow) == true ) {
                    graph = gM->GetIDBCGraph(sbmCol,sbmRow);

                    sbmMat->SetSubMatrix( sbmCol, sbmRow, eType, BaseMatrix::SPARSE_NONSYM,
                                          graph->GetSize(),
                                          graph->GetNumColsMat(),
                                          graph->GetNNE() );
                    LOG_DBG3(idbcElim) << "Transposed Entry";
                    LOG_DBG3(idbcElim) << "\t\t (sbmCol,sbmRow)=(" << sbmCol << ", " << sbmRow << "), size "
                                        << graph->GetSize() << " x " << graph->GetNumColsMat()
                                        << " with " <<  graph->GetNNE()  << " entries";
                    sbmMat->GetPointer(sbmCol,sbmRow)->SetSparsityPattern( *graph );
                  }
                } // if transpose
              } // if exists
            }//if
          } // sbmCol loop
        } // sbmRow loop

      // Insert pointer into BaseMatrix array
      auxMat_[*it] = sbmMat;
      LOG_DBG3(idbcElim) << "auxMat_[*it]->ToString()"<< std::endl;
      LOG_DBG3(idbcElim) << auxMat_[*it]->ToString();

    } // loop over matrix types

    //Those here are getting set to 0. It must be elsewhere
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

  template <typename T>
  void IDBC_Handler<T>::BuildSystemMatrix( const std::map<FEMatrixType, Complex> &factors, std::map<UInt, std::set<UInt> >& colInd ) {
    EXCEPTION("BuildSystemMatrix not implemented for complex factors");
  }

  // ****************
  //   AddIDBCToRHS
  // ****************
  template <typename T> void IDBC_Handler<T>::
  AddIDBCToRHS( SBM_Vector *rhs, bool deltaIDBC ) {
    
    LOG_DBG3(idbcElim) << "vecIDBC_"<< std::endl;
    LOG_DBG3(idbcElim)<< vecIDBC_->ToString();
    if ( addIDBCPossible_ == false ) {
      EXCEPTION( "IDBCHandler::AddIDBCToRHS: Internal error! Refusing to "
          << "add Dirichlet BCs, since addIDBCPossible_ = false! "
          << "Did you call BuildSystemMatrix()?");
    }
    LOG_DBG3(idbcElim) << "AddIDBCToRHS";
    LOG_DBG3(idbcElim) << "Old IDBC Values: " << vecOldIDBC_->ToString() ;
    LOG_DBG3(idbcElim) << "New IDBC Values: " << vecIDBC_->ToString();
    LOG_DBG3(idbcElim) << "old RHS: " << rhs->ToString();

    auxMat_[SYSTEM]->MultSub( *vecIDBC_, *rhs );
    if(deltaIDBC){
      auxMat_[SYSTEM]->MultAdd( *vecOldIDBC_, *rhs );
    }

    LOG_DBG3(idbcElim) << "RHS with IDBC: " << rhs->ToString();

    remIDBCPossible_ = true;


// This is just for exporting the auxMat_
//    UInt size = 7;
//    // same as ComputeIndex method in GraphManager, here with a lambda function
//    auto ComputeIndex = [size](UInt a, UInt b ) { return (size) * a + b;};
//
//    // store the sbm-indices of the blocks, which have to be assembled
//    // NOTE: If we are using the optimized version, we only consider odd harmonics,
//    // therefore the column isn't anymore iRow+harmonic !!!
//    // =========== This is difficult to describe here, see notes (K.Roppert 2018) =====
//    StdVector<UInt> sbmInd;
//    bool isFullSys = true;
//    for (UInt harm = 0; harm < 7; ++harm){
//      for( UInt iRow = 0; iRow < size; ++iRow ) {
//        if(harm == 0){
//          // sbm-blocks on the main diagonal
//          sbmInd.Push_back( ComputeIndex(iRow, iRow) );
//        }else{
//          Integer col = (isFullSys)? (iRow - harm) : (iRow - harm/2);
//          if( col < (Integer)size && col >= 0){
//            // change row/columns
//            sbmInd.Push_back( ComputeIndex(col, iRow) );
//          }
//        }
//      }
//    }
//    std::string path = "/home/alex/Documents/Doktorarbeit/PiezoHysteresis/Testcases/DirecletBC/SystemMatrix/auxMat";
//    auxMat_[SYSTEM]->Export_MultHarm(path, BaseMatrix::MATRIX_MARKET, 3, sbmInd);
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

  template <typename T>
  void IDBC_Handler<T>::PrintIDBCvec(){
    std::cout << "IDBC-Vector" << std::endl;
    std::cout << vecIDBC_->ToString() << std::endl;
  }

  
  // Explicit template instantiation
  template class IDBC_Handler<Double>;
  template class IDBC_Handler<Complex>;
  
}
