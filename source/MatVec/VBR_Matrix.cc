#include "VBR_Matrix.hh"

#include "opdefs.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include "MatVec/BLASLAPACKInterface.hh"
#include "MatVec/Matrix.hh"

DEFINE_LOG(vbrMat, "vbrMatrix")

namespace CoupledField {


  void TestVBRMatrix();
  
  template<typename T>
  VBR_Matrix<T>::VBR_Matrix() {
    
    // initialize derived variables
    this->nnz_   = 0;
    this->ncols_ = 0;
    this->nrows_ = 0;
    
    nNzEff_  = 0;
    nbRows_  = 0;
    nbCols_  = 0;
    nBlocks_ = 0;
    oneOverBlockSize_ = 0.0;
    numNonDiagBlockEntires_ = 0;
    
    bRow_         = NULL;
    bCol_         = NULL;
    valPtr_       = NULL;
    colInd_       = NULL;
    rowPtr_       = NULL;
    diagPtr_      = NULL;
    diagBlockPtr_ = NULL;
    data_         = NULL;
  
  }

  template<typename T>
  VBR_Matrix<T>::VBR_Matrix( const VBR_Matrix<T> &origMat ) {
    
    this->nnz_   = origMat.nnz_;
    this->ncols_ = origMat.ncols_;
    this->nrows_ = origMat.nrows_;

    nNzEff_  = origMat.nNzEff_;
    nbRows_  = origMat.nbRows_;
    nbCols_  = origMat.nbCols_;
    nBlocks_ = origMat.nBlocks_;
    oneOverBlockSize_ = origMat.oneOverBlockSize_;
    numNonDiagBlockEntires_ = origMat.numNonDiagBlockEntires_;

    // Create empty data structurs
    NEWARRAY( bRow_,   UInt, nbRows_ + 1 );
    NEWARRAY( rowPtr_, UInt, nbRows_ + 1 );
    
    NEWARRAY( diagBlockPtr_, UInt, nbRows_ + 1 );
    for(UInt i=0; i < nbRows_+1; ++i ) {
      bRow_[i] = origMat.bRow_[i];
      rowPtr_[i] = origMat.rowPtr_[i];
      diagBlockPtr_[i] = origMat.diagBlockPtr_[i];
    }
    
    NEWARRAY( bCol_,   UInt, nbCols_ + 1 );
    for(UInt i=0; i < nbCols_+1; ++i ) {
      bCol_[i] = origMat.bCol_[i];
    }

    NEWARRAY( valPtr_, UInt, nBlocks_ + 1);
    for(UInt i=0; i < nBlocks_+1; ++i ) {
      valPtr_[i] = origMat.valPtr_[i];
    }

    NEWARRAY( colInd_, UInt, nBlocks_);
    for(UInt i=0; i < nBlocks_; ++i ) {
      colInd_[i] = origMat.colInd_[i];
    }

    NEWARRAY( data_, T, this->nNzEff_ );
    for(UInt i=0; i < this->nNzEff_ ; ++i ) {
      data_[i] = origMat.data_[i];
    }
    
    // Only copy diagonal pointer, if matrix is quadratic
    if( this->ncols_ == this->nrows_ ) {
      NEWARRAY( diagPtr_, UInt, this->nrows_ );
      for(UInt i=0; i < this->nrows_ ; ++i ) {
        diagPtr_[i] = origMat.diagPtr_[i];
      }
    }
  }

  template<typename T>
    VBR_Matrix<T>::VBR_Matrix(UInt nrows, UInt ncols, UInt nnz ) {
    
    this->nnz_   = nnz;
    this->ncols_ = ncols;
    this->nrows_ = nrows;

    nNzEff_  = 0;
    nbRows_  = 0;
    nbCols_  = 0;
    nBlocks_ = 0;
    oneOverBlockSize_ = 0.0;
    numNonDiagBlockEntires_ = 0;
    
    bRow_    = NULL;
    bCol_    = NULL;
    valPtr_  = NULL;
    colInd_  = NULL;
    rowPtr_  = NULL;
    diagPtr_ = NULL;
  }
  
  template<typename T>
  VBR_Matrix<T>::~VBR_Matrix( ) {
    delete[] bRow_;
    bRow_    = NULL;

    delete[] bCol_;
    bCol_    = NULL;

    delete[] valPtr_;
    valPtr_  = NULL;
    
    delete[] colInd_;
    colInd_  = NULL;
    
    delete[] rowPtr_;
    rowPtr_  = NULL;
    
    delete[] diagPtr_;
    diagPtr_ = NULL;
    
    delete[] diagBlockPtr_;
    diagBlockPtr_ = NULL;
    
    delete[] data_;
    data_    = NULL;
  }


  template<typename T>
  void VBR_Matrix<T>::SetSparsityPattern( BaseGraph &graph ) {

    // Obtain block definition from graph object
    StdVector<std::pair<UInt,UInt> >& rBlockInfo 
    = graph.GetRowBlockDefinition();

    StdVector<std::pair<UInt,UInt> >& cBlockInfo 
    = graph.GetColBlockDefinition();
    
    LOG_DBG(vbrMat) << "================================\n";
    LOG_DBG(vbrMat) << "VBR - Block Definition\n";
    LOG_DBG(vbrMat) << "================================\n";
    
    // Check row blocks
    // If no definition is available, we assume block size 1 for row blocks
    if( rBlockInfo.GetSize() > 0 ) {
      nbRows_ = rBlockInfo.GetSize();
      NEWARRAY( bRow_,   UInt, nbRows_ + 1 );
      NEWARRAY( rowPtr_, UInt, nbRows_ + 1 );
      NEWARRAY( diagBlockPtr_, UInt, nbRows_ + 1 );
      LOG_DBG(vbrMat) << "1) ROW-BLOCKS:";
      for( UInt i = 0; i < nbRows_; ++i ) {
        bRow_[i] = rBlockInfo[i].first;
        LOG_DBG(vbrMat) << "\t#" << i << ": " << rBlockInfo[i].first << ", " 
            << rBlockInfo[i].second << std::endl;
      }
    }else {
      // generate block-size 1 system
      nbRows_ = this->nrows_;
      NEWARRAY( bRow_,   UInt, nbRows_ + 1 );
      NEWARRAY( rowPtr_, UInt, nbRows_ + 1 );
      NEWARRAY( diagBlockPtr_, UInt, nbRows_ + 1 );
      for( UInt i = 0; i < nbRows_; ++i ) {
        bRow_[i] = i;
      }
    }

    // Check column blocks:
    // If no definition is available, we assume block size 1 for col blocks
    if( cBlockInfo.GetSize() > 0 ) {
      nbCols_ = cBlockInfo.GetSize();  
      NEWARRAY( bCol_,   UInt, nbCols_ + 1 );
      LOG_DBG(vbrMat) << "2) COL-BLOCKS:";
      for( UInt i = 0; i < nbCols_; ++i ) {
        bCol_[i] = cBlockInfo[i].first;
        LOG_DBG(vbrMat) << "\t#" << i << ": " << cBlockInfo[i].first << ", " 
            << cBlockInfo[i].second << std::endl;
      }
    } else {
      // generate block-size 1 system
      nbCols_ = this->ncols_;  
      NEWARRAY( bCol_,   UInt, nbCols_ + 1 );
      for( UInt i = 0; i < nbCols_; ++i ) {
        bCol_[i] = i;
      }
    }
    
    bRow_[nbRows_] = this->nrows_;
    bCol_[nbCols_] = this->ncols_;

    // As we do not know the absolute number of blocks yet, we
    // have to create temporary arrays for valPtr array.
    StdVector<UInt> valPtr, colInd;
    valPtr.Reserve(this->nrows_ * 10 );
    colInd.Reserve(this->nrows_ * 10 );

    // create counters
    UInt bRowStart = 0, bRowEnd = 0, bRowSize = 0;
    UInt bColEnd = 0;
    nBlocks_ = 0;
    UInt actInd = 0; // pointer for current entry in valPtr array
    UInt accBlockSize = 0;

    // do a loop over all blockRows
    for( UInt ibRow = 0; ibRow < nbRows_; ++ibRow ) {

      // set pointer to start of this block
      rowPtr_[ibRow] = colInd.GetSize();

      // determine current row block size
      bRowStart = bRow_[ibRow];
      bRowEnd = bRow_[ibRow+1]-1;
      bRowSize =  bRowEnd - bRowStart+1;
      accBlockSize += bRowSize;

      // create temporary set for used column blocks in this row
      std::set<UInt> usedColBlocks;

      // loop over all rows within this block
      for( UInt actRow = bRowStart; actRow <= bRowEnd; ++actRow ) {

        // counter for current block column
        UInt ibCol = 0;

        // determine current row block size
        bColEnd  = bCol_[ibCol+1]-1;

        // get edges  for current row from graph
        UInt * cols= graph.GetGraphRow(actRow);
        UInt numColEntries = graph.GetRowSize(actRow);

        // loop over all column entries of this row
        for( UInt iCol = 0; iCol < numColEntries; ++iCol ) {
          const UInt & actCol = cols[iCol];
          while( actCol > bColEnd ) {
            bColEnd   = bCol_[++ibCol+1]-1;
          }
          usedColBlocks.insert(ibCol);
          
          // count number of off-block-diagonal entries
          if( ibCol != ibRow )
            numNonDiagBlockEntires_++;
            
        } // loop over column entries
      } // loop over rows of block


      // now we know all referenced column blocks for this block-row
      std::set<UInt>::iterator it = usedColBlocks.begin();
      for( ; it != usedColBlocks.end(); ++it ) {
        colInd.Push_back(*it);
        valPtr.Push_back(actInd);

        // check, if we have a diagonal block and remeber it for later use
        if( *it == ibRow) {
          diagBlockPtr_[ibRow] = actInd;
        }
        actInd += bRowSize * (bCol_[(*it)+1]-bCol_[*it]);
      }
    }

    // Now we know the total number of blocks...
    nBlocks_ = colInd.GetSize();
    rowPtr_[nbRows_] = nBlocks_;
    
    // .. as well as the "effective" number of nonzeros 
    // (including padding zeros)
    this->nNzEff_ =  actInd;
    valPtr.Push_back(this->nnz_);
    
    // move temporary arrays to final ones
    NEWARRAY( valPtr_, UInt, nBlocks_ + 1);
    NEWARRAY( colInd_, UInt, nBlocks_);
    for( UInt i = 0; i < nBlocks_; ++i ) {
      valPtr_[i] = valPtr[i];
      colInd_[i] = colInd[i];
    }
    valPtr_[nBlocks_] = this->nNzEff_;
    
    // Empty temporary vector (before creating memory for the
    // data array)
    valPtr.Clear();
    colInd.Clear();

    // Reserve memory for nonzeros and initialize it
    NEWARRAY( data_, T, this->nNzEff_ );
    Init();

    // Store average block size.
    // This value is used to accelerate searching for rows in the
    // FindBlock()-method
    oneOverBlockSize_ = float(1.0) / (float(accBlockSize) / nbRows_);

    // find diagonal entries, if this is a square matrix
    if( this->nrows_ == this->ncols_ ) {
      NEWARRAY( diagPtr_, UInt, this->nrows_ );
      FindDiagonalEntries();
    }
  }

  template<typename T>
  void VBR_Matrix<T>::SetSparsityPattern( const StdMatrix &srcMat ) {

    EXCEPTION( "VBR_Matrix::SetSparsityPattern: Not implemented" );
  }
  
  template<typename T>
  void VBR_Matrix<T>::SetSize(UInt nrows, UInt ncols, UInt nnz) {
    EXCEPTION("Implement me");
  }
  
  template<typename T>
  void VBR_Matrix<T>::Init() {
    for ( UInt i = 0; i < this->nNzEff_; i++ ) {
      data_[i] = 0.0;
    }
  }

  template<typename T>
  inline void VBR_Matrix<T>::Mult( const Vector<T> &mvec,
                                   Vector<T> &rvec ) const {

    UInt rStart; // start index of row
    UInt cStart; // start index of col
    UInt rbs;    // row block size
    UInt cbs;    // col block size
    UInt colNum; // column number
    UInt ind;    //index to data array


    rvec.Init();
    //     std::cerr << "MAT-VEC-Mulitplication\n";
    //     std::cerr << "=======================\n";
    // loop over rowblock
    for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
      rStart = bRow_[ibr];
      rbs = bRow_[ibr+1] -  bRow_[ibr];
      //      std::cerr << "* rowBlock #" << ibr << std::endl;
      //      std::cerr << "\trStart: " << rStart << std::endl;
      //      std::cerr << "\trow bs: " << rbs << std::endl;


      // loop over col blocks
      for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {

        colNum = colInd_[ibc];
        cStart = bCol_[colNum];
        cbs = bCol_[colNum+1] - bCol_[colNum];
        //        std::cerr << "\t colBlock #" << ibc << std::endl;
        //        std::cerr << "\t\tccolNum: " << colNum << std::endl;
        //        std::cerr << "\t\tcStart: " << cStart << std::endl;
        //        std::cerr << "\t\tcol bs: " << cbs << std::endl;

        // perform mat-vec multiplication on sub-block
        ind = valPtr_[ibc];

        for( UInt i = rStart; i < rStart+rbs; ++i ) {
          for( UInt j = cStart; j <cStart+cbs; ++j ) {
            //            std::cerr << "\t\t\t revc[" << i << "] += data_["
            //                      << ind << "] * mvec[" << j << "]\n";
            //            std::cerr << "\t\t\t revc[" << i << "] += " << data_[ind]
            //                      << "* " << mvec[j] << "\n\n";
              rvec[i] += data_[ind++] * mvec[j];
              //              std::cerr << "\t\t\t revc[" << i << "] += " << data_[ind]
              //                        << "* " << mvec[j] << "\n\n";
          }  // loop over cols within block
        } // loop over rows within block
      } //block cols
    } // block rows
  
  }
  
  template<typename T>
  inline void VBR_Matrix<T>::MultAdd( const Vector<T> & mvec,
                                      Vector<T> & rvec ) const {
    EXCEPTION("General case not implemented");
    
  }

  template<>
  inline void VBR_Matrix<Double>::MultAdd( const Vector<Double> & mvec,
                                      Vector<Double> & rvec ) const {
      int rStart; // start index of row
      int cStart; // start index of col
      int rbs;    // row block size
      int cbs;    // col block size
      UInt colNum; // column number
      UInt ind;    //index to data array
      
#ifdef NDEBUG
      char trans = 'T';
      Double alpha = 1.0;
      Double beta = 1.0;
      Integer inc = 1;
#endif
      
      // loop over row blocks
      for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
        rStart = bRow_[ibr];
        rbs = bRow_[ibr+1] -  bRow_[ibr];

        // loop over col blocks
        for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {
          colNum = colInd_[ibc];
          cStart = bCol_[colNum];
          cbs = bCol_[colNum+1] - bCol_[colNum];
          ind = valPtr_[ibc];

          // perform mat-vec multiplication on dense sub-block
#ifdef NDEBUG
          dgemv( &trans, &cbs, &rbs, &alpha, &(data_[ind]), &cbs, 
                 &mvec[cStart], &inc, &beta, &rvec[rStart], &inc);
          ind+=cbs*rbs;
#else
          for( int i = rStart; i < rStart+rbs; ++i ) {
            for( int j = cStart; j <cStart+cbs; ++j ) {
              rvec[i] += data_[ind++] * mvec[j];
            }  // loop over cols within block
          } // loop over rows within block

#endif
        } // block cols
      } // block rows
    }

  template<typename T>
    inline void VBR_Matrix<T>::MultSub( const Vector<T> &mvec,
                                        Vector<T> &rvec ) const {
    EXCEPTION("General case not implemented");
  }

  template<>
  inline void VBR_Matrix<Double>::MultSub( const Vector<Double> &mvec,
                                           Vector<Double> &rvec ) const {

    int rStart; // start index of row
    int cStart; // start index of col
    int rbs;    // row block size
    int cbs;    // col block size
    UInt colNum; // column number
    UInt ind;    //index to data array

#ifdef NDEBUG
    char trans = 'T';
    Double alpha = -1.0;
    Double beta = 1.0;
    Integer inc = 1;
#endif

    // loop over row blocks
    for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
      rStart = bRow_[ibr];
      rbs = bRow_[ibr+1] -  bRow_[ibr];

      // loop over col blocks
      for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {
        colNum = colInd_[ibc];
        cStart = bCol_[colNum];
        cbs = bCol_[colNum+1] - bCol_[colNum];
        ind = valPtr_[ibc];
        
        // perform mat-vec multiplication on dense sub-block
#ifdef NDEBUG
        dgemv( &trans, &cbs, &rbs, &alpha, &(data_[ind]), &cbs, 
               &mvec[cStart], &inc, &beta, &rvec[rStart], &inc);
        ind+=cbs*rbs;
#else
        for( int i = rStart; i < rStart+rbs; ++i ) {
          for( int j = cStart; j <cStart+cbs; ++j ) {
            rvec[i] -= data_[ind++] * mvec[j];
          }  // loop over cols within block
        } // loop over rows within block

#endif
      } // block cols
    } // block rows
  }


  template<typename T>
  inline void VBR_Matrix<T>::CompRes( Vector<T> &r, const Vector<T> &x,
                                      const Vector<T> &b ) const {
    
    EXCEPTION("General case not implemented");
  }
  
  template<>
  inline void VBR_Matrix<Double>::CompRes( Vector<Double> &r, 
                                           const Vector<Double> &x,
                                           const Vector<Double> &b ) const {
    int rStart; // start index of row
    int cStart; // start index of col
    int rbs;    // row block size
    int cbs;    // col block size
    UInt colNum; // column number
    UInt ind;    //index to data array
    
#ifdef NDEBUG
    char trans = 'T';
    Double alpha = -1.0;
    Double beta = 1.0;
    Integer inc = 1;
#endif

    // initialize return vector
    r = b;
    
    // loop over row blocks
    for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
      rStart = bRow_[ibr];
      rbs = bRow_[ibr+1] -  bRow_[ibr];
      
      // loop over col blocks
      for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {
        colNum = colInd_[ibc];
        cStart = bCol_[colNum];
        cbs = bCol_[colNum+1] - bCol_[colNum];
        ind = valPtr_[ibc];
        
        // perform mat-vec multiplication on dense sub-block
#ifdef NDEBUG
        dgemv( &trans, &cbs, &rbs, &alpha, &(data_[ind]), &cbs, 
               &x[cStart], &inc, &beta, &r[rStart], &inc);
        ind+=cbs*rbs;
#else
        for( int i = rStart; i < rStart+rbs; ++i ) {
          for( int j = cStart; j <cStart+cbs; ++j ) {
            r[i] -= data_[ind++] * x[j];
          }  // loop over cols within block
        } // loop over rows within block
#endif
      } //block cols
    } // block rows
  }


  template<typename T>
  inline void VBR_Matrix<T>::MultT( const Vector<T> & mvec,
                                    Vector<T> & rvec ) const {
    EXCEPTION("Implement me");
  }

  template<typename T>
  inline void VBR_Matrix<T>::MultTAdd( const Vector<T> &mvec,
                                       Vector<T> &rvec ) const {
    EXCEPTION("Not implemented");
  }
    

  template<>
  inline void VBR_Matrix<Double>::MultTAdd( const Vector<Double> &mvec,
                                            Vector<Double> &rvec ) const {
    int rStart; // start index of row
    int cStart; // start index of col
    int rbs;    // row block size
    int cbs;    // col block size
    UInt colNum; // column number
    UInt ind;    //index to data array

#ifdef NDEBUG
    char trans = 'N';
    Double alpha = 1.0;
    Double beta = 1.0;
    Integer inc = 1;
#endif

    // loop over row blocks
    for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
      rStart = bRow_[ibr];
      rbs = bRow_[ibr+1] -  bRow_[ibr];

      // loop over col blocks
      for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {
        colNum = colInd_[ibc];
        cStart = bCol_[colNum];
        cbs = bCol_[colNum+1] - bCol_[colNum];
        ind = valPtr_[ibc];

        //perform mat-vec multiplication on dense sub-block
#ifdef NDEBUG
        dgemv( &trans, &cbs, &rbs, &alpha, &(data_[ind]), &cbs, 
               &mvec[rStart], &inc, &beta, &rvec[cStart], &inc);
        ind+=cbs*rbs;
#else
        for( int i = rStart; i < rStart+rbs; ++i ) {
          for( int j = cStart; j <cStart+cbs; ++j ) {
            rvec[j] += data_[ind++] * mvec[i];
          }  // loop over cols within block
        } // loop over rows within block
#endif

      } // block cols
    } // block rows
  }


  template<typename T>
  inline void VBR_Matrix<T>::MultTSub( const Vector<T> &mvec,
                                       Vector<T> &rvec ) const {
    EXCEPTION("Not implemented")

  }
    
  template<>
  inline void VBR_Matrix<Double>::MultTSub( const Vector<Double> &mvec,
                                            Vector<Double> &rvec ) const {
    int rStart; // start index of row
   int cStart; // start index of col
   int rbs;    // row block size
   int cbs;    // col block size
   UInt colNum; // column number
#ifdef NDEBUG
   UInt ind;    //index to data array

   char trans = 'N';
   Double alpha = -1.0;
   Double beta = 1.0;
   Integer inc = 1;
#endif

   // loop over row blocks
   for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
     rStart = bRow_[ibr];
     rbs = bRow_[ibr+1] -  bRow_[ibr];

     // loop over col blocks
     for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {
       colNum = colInd_[ibc];
       cStart = bCol_[colNum];
       cbs = bCol_[colNum+1] - bCol_[colNum];
#ifdef NDEBUG
       ind = valPtr_[ibc];
       
       //perform mat-vec multiplication on dense sub-block
       dgemv( &trans, &cbs, &rbs, &alpha, &(data_[ind]), &cbs, 
              &mvec[rStart], &inc, &beta, &rvec[cStart], &inc);
       ind+=cbs*rbs;
#else
        for( int i = rStart; i < rStart+rbs; ++i ) {
          for( int j = cStart; j <cStart+cbs; ++j ) {
          }  // loop over cols within block
        } // loop over rows within block
#endif
       // alternative: BLAS version

     } // block cols
   } // block rows
  }

  template<typename T>
  void VBR_Matrix<T>::Transpose(UInt* col_ptr, UInt* row_ptr, T* data_ptr) const
  {
    EXCEPTION( "Not implemented" );
  }


  template<typename T>
  std::string VBR_Matrix<T>::ToString( char colSeparator,
                                       char rowSeparator) const {
    EXCEPTION("Implement me");
  }

  template<typename T>
  void VBR_Matrix<T>::ExportMatrixMarket(const std::string& fname, const std::string& comment ) const{
    // Open output file and check for errors
    FILE *fp = fopen( fname.c_str(), "w" );
    if ( fp == NULL ) {
      EXCEPTION( "Cannot open file " << fname << " for writing!" );
    }

    // ---------------------
    //   Write file header
    // ---------------------

    // Matrix Market Format Specification
    if ( this->GetEntryType() == BaseMatrix::DOUBLE ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
    }
    else if ( this->GetEntryType() == BaseMatrix::COMPLEX ) {
      fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
    }
    else {
      EXCEPTION( "Expected either DOUBLE or COMPLEX as matrix entry type" );
    }

    // User-supplied private comment
    if ( comment != "" ) {
      fprintf( fp, "%%\n%% %s\n%%\n", comment.c_str() );
    }
    else {
      fprintf( fp, "%%\n%% Matrix exported by openCFS\n%%\n" );
    }

    // Information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", this->nrows_ , this->ncols_ ,
             this->nNzEff_ );


    // ------------------------
    //   Write matrix entries
    // ------------------------
    
    UInt rStart, cStart, rbs, cbs, colNum, ind;
    // loop over rowblock
    for( UInt ibr = 0; ibr < nbRows_; ++ibr ) {      
      rStart = bRow_[ibr];
      rbs = bRow_[ibr+1] -  bRow_[ibr];

      // loop over col blocks
      for( UInt ibc = rowPtr_[ibr]; ibc < rowPtr_[ibr+1]; ++ibc) {
        colNum = colInd_[ibc];
        cStart = bCol_[colNum];
        cbs = bCol_[colNum+1] - bCol_[colNum];
        ind = valPtr_[ibc];

        // write out dense sub-block
        for( UInt i = rStart; i < rStart+rbs; ++i ) {
          for( UInt j = cStart; j <cStart+cbs; ++j ) {

            // store row and column index
            fprintf( fp, "%6d\t%6d\t", i  + 1, j + 1);
            
            // store non-zero entry
            OpType<T>::ExportEntry( this->data_[ind++],
                                          0, 0, fp );
            fprintf( fp, "\n" );
        }  // loop over cols within block
      } // loop over rows within block
    } //block cols
  } // block rows
    
       // ------------------
    //   Finish writing 
    // ------------------

    // close output file
    if ( fclose( fp ) == EOF ) {
      WARN("Could not close file " << fname << " after writing!");
    }
  }


  template<typename T> 
  void VBR_Matrix<T>::FindBlock( UInt row, UInt col, UInt& rowB, 
                                UInt& colB, UInt& offset ) const {
    UInt l = 0;
    UInt u = nbRows_-1; 
    UInt k;
    bool found = false;
    offset = 0;
    UInt colSize;
    
    // use initial guess for k: multiply row by avg. blockwidth to get good starting value
    k = UInt(oneOverBlockSize_ * row );
    
    // look for row
    while( l <= u ) {
      if( bRow_[k] <=  row && bRow_[k+1] > row) {
        rowB = k;
        offset = (row-bRow_[k]);
        break;
      } else if (bRow_[k] > row ){
        u = k - 1;
      } else {
        l = k + 1;
      }
      k = (l+u) >> 1; // half the interval
    }
    // look for col
    l = rowPtr_[rowB];
    u = rowPtr_[rowB+1]-1;
    while( l <= u ) {
      k = (l+u) >> 1; // half the interval
      if( bCol_[colInd_[k]] <=  col && bCol_[colInd_[k]+1] > col) {
        found = true;
        colB = colInd_[k];
        colSize = bCol_[colInd_[k]+1] - bCol_[colInd_[k]];

        // determine offset within current block 
        offset = (col - bCol_[colInd_[k]])+offset*colSize; 
        
        // add offset of block start
        offset += valPtr_[k];
        break;
      } else if (bCol_[colInd_[k]] > col ){
        u = k - 1;
      } else {
        l = k + 1;
      }
    }
    if( !found ) {
      EXCEPTION("Index pair (" << row << ", " << col << ") not found");
    }
    
  }

  
  template<typename T>
  void VBR_Matrix<T>::AddToMatrixEntry( UInt i, UInt j, const T& v ){

    
    UInt rowB, colB, offset;
    // get correct block (double search)
    FindBlock(i,j,rowB,colB, offset);
    data_[offset] += v;
  }


  template<typename T>
  void VBR_Matrix<T>::SetMatrixEntry( UInt i, UInt j, const T &v,
                                      bool setCounterPart ) {
    
    UInt rowB, colB, offset;
    // get correct block (double search)
    FindBlock(i,j,rowB,colB, offset);
    data_[offset] = v;
  }
  
  template<typename T>
  void VBR_Matrix<T>::GetMatrixEntry( UInt i, UInt j,
                                      T &v ) const {
    UInt rowB, colB, offset;
    // get correct block (double search)
    FindBlock(i,j,rowB,colB, offset);
    v =  data_[offset];
  }


  template<typename T>
  void VBR_Matrix<T>::Insert(UInt row, UInt col,  UInt pos){
    
    
    EXCEPTION("Implement me");
  }
  
  template<typename T>
  void VBR_Matrix<T>::GetNumBlocks(UInt& nRowBlocks, 
                                   UInt& nColBlocks, 
                                   UInt& numBlocks ) const {
    nRowBlocks = nbRows_;
    nColBlocks = nbCols_;
    numBlocks = nBlocks_;
  }
  
  
  template <typename T>
  Double VBR_Matrix<T>::GetMemoryUsage() const {
    
    // sum up contributions
    Double sum = 0;
    
    // indices
    sum += ( (nbRows_ + 1) * 3 // bRow_, rowPtr_, diagBlockPtr
           + (nbCols_ + 1)     // bCol_
           + (nBlocks_ + 1)    // valPtr
           + (nBlocks_) )       // colInd_
           * sizeof(UInt);
    if( this->ncols_ == this->nrows_ ) { 
      sum += this->nrows_ * sizeof(UInt);
    }
    
    // data array
    sum += this->nNzEff_ * sizeof(T);
    
    return sum;
  }
  
  template<typename T>
  void VBR_Matrix<T>::GetDiagBlock( UInt blockRow, 
                                    DenseMatrix& diagBlock ) const {

    Matrix<T> & ret = static_cast<Matrix<T>& >(diagBlock);
    const T* intPt = &data_[diagBlockPtr_[blockRow]];
    UInt size = bRow_[blockRow+1] - bRow_[blockRow]; 
    ret.Resize( size );
    for( UInt i = 0; i < size; ++ i ) {
      for( UInt j = 0; j < size; ++ j ) {
        ret[i][j] = intPt[size*i+j];
      }
    }
  }
  
  template<typename T>
  void VBR_Matrix<T>::SetDiagBlock( UInt blockRow, 
                                    const DenseMatrix& diagBlock ) {

    const Matrix<T> & ret = static_cast<const Matrix<T>& >(diagBlock);
    T* intPt = &data_[diagBlockPtr_[blockRow]];
    // check for same size
    
    UInt size = bRow_[blockRow+1] - bRow_[blockRow]; 
    for( UInt i = 0; i < size; ++ i ) {
      for( UInt j = 0; j < size; ++ j ) {
         intPt[size*i+j] = ret[i][j];
      }
    }
  }
  
  template<typename T>
  void VBR_Matrix<T>::Scale( Double factor ) {
    EXCEPTION("Not implemented");
  }
  
  template<typename T>
  void VBR_Matrix<T>::Scale( Double factor,
                             const std::set<UInt>& rowIndices,
                             const std::set<UInt>& colIndices ) {
    EXCEPTION("Implement me");
  }


  template<typename T>
  void VBR_Matrix<T>::Add( const Double alpha, const StdMatrix& mat ) {
    EXCEPTION("Implement me");
  }

  template<typename T>
  void VBR_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat ) {
    EXCEPTION("Complex Add-subset is not implemented for VBR_Matrix.");
  }

  template<>
  void VBR_Matrix<Complex>::Add( const Double alpha, const StdMatrix& mat ) {

    EXCEPTION("Implement me");
  }

  template<>
  void VBR_Matrix<Complex>::Add( const Complex alpha, const StdMatrix& mat ) {
    EXCEPTION("Complex Add-subset is not implemented for VBR_Matrix.");
  }

  template<typename T>
  void VBR_Matrix<T>::Add( const Double alpha, const StdMatrix& mat,
                           const std::set<UInt>& rowIndices,
                           const std::set<UInt>& colIndices ) {
    EXCEPTION("Not implemented");
  }

  template<typename T>
  void VBR_Matrix<T>::Add( const Complex alpha, const StdMatrix& mat, const std::set<UInt>& rowIndices, const std::set<UInt>& colIndices ) {
    EXCEPTION("Complex Add-subset is not implemented for VBR_Matrix.");
  }

  template <typename T>
  void VBR_Matrix<T>::FindDiagonalEntries() {

    // loop over all rows
    UInt bRow = 0;
    UInt bCol = 0;
    UInt offset = 0;
    for( UInt i = 0; i < this->nrows_; ++i ) {
      FindBlock( i, i, bRow, bCol, offset );
      diagPtr_[i] = offset;
    }
  }

// Explicit template instantiation
  template class VBR_Matrix<Double>;
  template class VBR_Matrix<Complex>;

} // end of namespace

  // ==========================================================================
  //   TEST SECTION: VERIFY MATRIX IMPLEMENTATION
  // ==========================================================================
  // This test takes the reference matrix for the VBR-format from the guide
  //    "A reference model implementation of the Sparse BLAS in Fortran 95"
  //    http://dl.acm.org/citation.cfm?doid=567806.567811
  
#ifdef VBR_TEST_SECTION
#include "MatVec/Matrix.hh"
#include "MatVec/crs_Matrix.hh"
#include "MatVec/BLASLAPACKInterface.hh"

namespace CoupledField {
template<typename T>
 void FillMat( T& mat ) {
   // set entries
     mat.Init();
     
     // (0,0)
     mat.AddToMatrixEntry(0,0,4);
     mat.AddToMatrixEntry(0,1,2);
     mat.AddToMatrixEntry(1,0,1);
     mat.AddToMatrixEntry(1,1,5);
     
     // (0,2)
     mat.AddToMatrixEntry(0,5,1);
     mat.AddToMatrixEntry(1,5,2);
     
     // (0,4)
     mat.AddToMatrixEntry(0,9,-1);
     mat.AddToMatrixEntry(0,10,1);
     mat.AddToMatrixEntry(1,10,-1);
     
     // (1,1)
     mat.AddToMatrixEntry(2,2,6);
     mat.AddToMatrixEntry(2,3,1);
     mat.AddToMatrixEntry(2,4,2);
     mat.AddToMatrixEntry(3,2,2);
     mat.AddToMatrixEntry(3,3,7);
     mat.AddToMatrixEntry(3,4,1);
     mat.AddToMatrixEntry(4,2,-1);
     mat.AddToMatrixEntry(4,3,2);
     mat.AddToMatrixEntry(4,4,9);
     
     // (1,2)
     mat.AddToMatrixEntry(2,5,2);
     mat.AddToMatrixEntry(4,5,3);
     
     // (2,0)
     mat.AddToMatrixEntry(5,0,2);
     mat.AddToMatrixEntry(5,1,1);
     
     // (2,1)
     mat.AddToMatrixEntry(5,2,3);
     mat.AddToMatrixEntry(5,3,4);
     mat.AddToMatrixEntry(5,4,5);
     
     // (2,2)
     mat.AddToMatrixEntry(5,5,10);

     // (2,3)
     mat.AddToMatrixEntry(5,6,4);
     mat.AddToMatrixEntry(5,7,3);
     mat.AddToMatrixEntry(5,8,2);

     // (3,2)
     mat.AddToMatrixEntry(6,5,4);
     mat.AddToMatrixEntry(7,5,3);
     
     // (3,3)
     mat.AddToMatrixEntry(6,6,13);
     mat.AddToMatrixEntry(6,7,4);
     mat.AddToMatrixEntry(6,8,2);
     mat.AddToMatrixEntry(7,6,3);
     mat.AddToMatrixEntry(7,7,11);
     mat.AddToMatrixEntry(7,8,3);
     mat.AddToMatrixEntry(8,6,2);
     mat.AddToMatrixEntry(8,8,7);
     
     // (4,0)
     mat.AddToMatrixEntry(9,0,8);
     mat.AddToMatrixEntry(9,1,4);
     mat.AddToMatrixEntry(10,0,-2);
     mat.AddToMatrixEntry(10,1,3);
     
     // (4,4)
     mat.AddToMatrixEntry(9,9,25);
     mat.AddToMatrixEntry(9,10,3);
     mat.AddToMatrixEntry(10,9,8);
     mat.AddToMatrixEntry(10,10,12);
 }
 
 
 void TestVBRMatrix() {
   
   std::cerr << "==============================\n";
   std::cerr << "    Perform VBR matrix  test      \n";
   std::cerr << "==============================\n";
   
   // Generate new basegraph#
   
   std::cerr << "\t filling graph...\n";
   BaseGraph graph (11, 11, BaseOrdering::NOREORDERING);
   
   std::vector<UInt> row(1), edges;

   // === row 0 ===
   edges.clear();
   edges.push_back(0);
   edges.push_back(1);
   edges.push_back(5);
   edges.push_back(9);
   edges.push_back(10);
   row[0] = 0;
   graph.AddVertexNeighbours(row, edges);

   // === row 1 ===
   edges.clear();
   edges.push_back(0);
   edges.push_back(1);
   edges.push_back(5);
   edges.push_back(10);
   row[0] = 1;
   graph.AddVertexNeighbours(row, edges);


   // === row 2===
   edges.clear();
   edges.push_back(2);
   edges.push_back(3);
   edges.push_back(4);
   edges.push_back(5);
   row[0] = 2;
   graph.AddVertexNeighbours(row, edges);


   // === row 3 ===
   edges.clear();
   edges.push_back(2);
   edges.push_back(3);
   edges.push_back(4);
   row[0] = 3;
   graph.AddVertexNeighbours(row, edges);

   // === row 4 ===
   edges.clear();
   edges.push_back(2);
   edges.push_back(3);
   edges.push_back(4);
   edges.push_back(5);
   row[0] = 4;
   graph.AddVertexNeighbours(row, edges);

   // === row 5 ===
   edges.clear();
   edges.push_back(0);
   edges.push_back(1);
   edges.push_back(2);
   edges.push_back(3);
   edges.push_back(4);
   edges.push_back(5);
   edges.push_back(6);
   edges.push_back(7);
   edges.push_back(8);
   row[0] = 5;
   graph.AddVertexNeighbours(row, edges);

   // === row 6 ===
   edges.clear();
   edges.push_back(5);
   edges.push_back(6);
   edges.push_back(7);
   edges.push_back(8);
   row[0] = 6;
   graph.AddVertexNeighbours(row, edges);

   // === row 7 ===
   edges.clear();
   edges.push_back(5);
   edges.push_back(6);
   edges.push_back(7);
   edges.push_back(8);
   row[0] = 7;
   graph.AddVertexNeighbours(row, edges);

   // === row 8 ===
   edges.clear();
   edges.push_back(6);
   edges.push_back(8);
   row[0] = 8;
   graph.AddVertexNeighbours(row, edges);

   // === row 9 ===
   edges.clear();
   edges.push_back(0);
   edges.push_back(1);
   edges.push_back(9);
   edges.push_back(10);
   row[0] = 9;
   graph.AddVertexNeighbours(row, edges);
   
   // === row 10 ===
   edges.clear();
   edges.push_back(0);
   edges.push_back(1);
   edges.push_back(9);
   edges.push_back(10);
   row[0] = 10;
   graph.AddVertexNeighbours(row, edges);
   
   // Define blocks
   std::cerr << "\t setting blockinfo ....\n";
   StdVector<StdVector<UInt> > blocks;
   blocks.Resize(5);
   blocks[0] = 0,1;
   blocks[1] = 2,3,4;
   blocks[2] = 5;
   blocks[3] = 6,7,8;
   blocks[4] = 9,10;
   graph.SetBlockInfo(&blocks);
   
   // perform reordering
   StdVector<UInt> order;
   std::cerr << "\t finalising assembly ....\n";
   graph.FinaliseAssembly( false, &order, &order);
   std::cerr << "\t order: " << order.ToString() << std::endl;
   
   // print graph
   std::cerr << "\t Graph output:\n";
   graph.Print(std::cerr);
   
   // generate VBR-matrix and fill entries
   std::cerr << "\t generate VBR matrix & fill pattern ...\n";
   VBR_Matrix<Double> vbrMat(11,11,52);
   vbrMat.SetSparsityPattern(graph);
   FillMat( vbrMat);
   
   // generate CRS matrix as well for comparison
   CRS_Matrix<Double> crsMat(11,11,52);
   crsMat.SetSparsityPattern(graph);
   FillMat( crsMat);
   
   // Get maximum diagonal
   std::cerr << "VBR maxDiag = " << vbrMat.GetMaxDiag() << std::endl;
   std::cerr << "CRS maxDiag = " << crsMat.GetMaxDiag() << std::endl;
   
   // create rhs row vector and compare
   Vector<Double> solVec(11), rhsVec(11), retVec(11);
   for(UInt i = 0; i < 11; ++i ) {
     solVec[i] = i+1;
     rhsVec[i] = (i+1)*2;
   }
   
   // multiply mat * solVec
   std::cerr << "\nResult of mat * solVec:\n";
   vbrMat.Mult(solVec, retVec);
   std::cerr << "VBR: " << retVec.ToString() << std::endl;
   crsMat.Mult(solVec, retVec);
   std::cerr << "CRS: " << retVec.ToString() << std::endl;
   
   
   
   // compute residual:
   std::cerr << "\nRedidual r = rhsVec - mat*solVec:\n";
   vbrMat.CompRes(retVec, solVec, rhsVec);
   std::cerr << "VBR: " << retVec.ToString() << std::endl;
   crsMat.CompRes(retVec, solVec, rhsVec);
   std::cerr << "CRS: " << retVec.ToString() << std::endl;
   
   // use MultAdd
   Vector<Double> temp;
   temp = retVec;
   std::cerr << "\n\nMultAdd r += mat * solVec\n";
   vbrMat.MultAdd(solVec, retVec);
   std::cerr << "VBR: " << retVec.ToString() << std::endl;
   retVec = temp;
   crsMat.MultAdd(solVec, retVec);
   std::cerr << "CRS: " << retVec.ToString() << std::endl;
   
   // use MultAdd-Transposed
   temp = retVec;
   std::cerr << "\n\nMultAdd r += mat^T * solVec\n";
   vbrMat.MultTAdd(solVec, retVec);
   std::cerr << "VBR: " << retVec.ToString() << std::endl;
   retVec = temp;
   crsMat.MultTAdd(solVec, retVec);
   std::cerr << "CRS: " << retVec.ToString() << std::endl;
   
   // use MultSub-Transposed
   temp = retVec;
   std::cerr << "\n\nMultSub  r -= mat^T * solVec\n";
   vbrMat.MultTSub(solVec, retVec);
   std::cerr << "VBR: " << retVec.ToString() << std::endl;
   retVec = temp;
   crsMat.MultTSub(solVec, retVec);
   std::cerr << "CRS: " << retVec.ToString() << std::endl;
   
   
   
   // ========================
   //  INVERT GENERAL MATRIX
   // ========================
   std::cerr << "===================================\n";
   std::cerr << " TEST INVERSION OF DIAGONAL BLOCKS \n";
   std::cerr << "===================================\n";

   UInt size, rowStart;
   
   Double * block = NULL;

   for( UInt diagRow = 0; diagRow < 5; ++diagRow  ) {
     block = vbrMat.GetDiagBlock( diagRow, size, rowStart);
     std::cerr << "\ndiagBlock [" << diagRow << "] is\n"; 
     for( UInt i = 0; i < size; ++i ) {
       for(UInt j = 0; j < size; ++j ) {
         std::cerr << block[i*size+j] << ", ";
       }
       std::cerr << std::endl;
     }
     
     double * fac = new double[size*size];
     std::copy(block,block+size*size, fac);
     int N = size;
     int *IPIV = new int[N+1];
     int LWORK = N*N;
     double *WORK = new double[LWORK];
     int INFO;

     dgetrf(&N,&N,fac,&N,IPIV,&INFO);
     dgetri(&N,fac,&N,IPIV,WORK,&LWORK,&INFO);

     delete IPIV;
     delete WORK;

     
     
     // print inverse
     std::cerr << "inverse is :\n";
     for( UInt i = 0; i < size; ++i ) {
       for(UInt j = 0; j < size; ++j ) {
         std::cerr << fac[i*size+j] << ", ";
       }
       std::cerr << std::endl;
     }    
     
   }
   
   
//   // ========================
//   //  TEST BLAS DGEMMV
//   // ========================
//   
//   // test 
//   {
//     Matrix<Double> a(2,2);
//     a[0][0] = 1.0;
//     a[0][1] = 2.0;
//     a[1][0] = 3.0;
//     a[1][1] = 4.0;
//     Vector<Double> r(4), s(4);
//     s.Init(1.0);
//     r.Init(0.0);
//     
//     char trans = 'T';
//     int size = 2;
//     double alpha = 1.0;
//     double beta = 0.0;
//     int inc = 1;
//     
//     DGEMV(&trans, 
//           &size,
//           &size,
//           &alpha,
//           a[0],
//           &size,
//           &s[1], 
//           &inc,
//           &beta,  &r[1], &inc);
//     std::cerr << "r = " << r.ToString() << std::endl;
//
//   }
  
 }
}

#endif
