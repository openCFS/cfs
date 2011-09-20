// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// the following headers are required for Export()
#include <cstdio>
#include <sstream>


#include "MatVec/vbr_matrix.hh"
#include "opdefs.hh"
#include "DataInOut/Logging/cfslog.hh"

DECLARE_LOG(vbrMat)
DEFINE_LOG(vbrMat, "vbrMatrix")

namespace CoupledField {


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
    
    EXCEPTION("Implement me");
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
    StdVector<std::pair<UInt,UInt> >& blockInfo 
    = graph.GetBlockDefinition();
    
    

    // we know already numbers of block rows / cols
    nbRows_ = blockInfo.GetSize();
    nbCols_ = blockInfo.GetSize();

    // we know already the sizes of the blocks. Assuming, that
    // every diagonal block is present, we can initialize
    // bRow_, bCols_ 
    NEWARRAY( bRow_,   UInt, nbRows_ + 1 );
    NEWARRAY( bCol_,   UInt, nbCols_ + 1 );
    NEWARRAY( rowPtr_, UInt, nbRows_ + 1 );
    NEWARRAY( diagBlockPtr_, UInt, nbRows_ + 1 );

    
    LOG_DBG(vbrMat) << "================================\n";
    LOG_DBG(vbrMat) << "VBR - Block Definition\n";
    LOG_DBG(vbrMat) << "================================\n";

    for( UInt i = 0; i < nbRows_; ++i ) {
      bRow_[i] = blockInfo[i].first;
      bCol_[i] = blockInfo[i].first;
      LOG_DBG(vbrMat) << "\t#" << i << ": " << blockInfo[i].first << ", " 
          << blockInfo[i].second << std::endl;
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
    oneOverBlockSize_ = 1.0 / (accBlockSize / nbRows_);

    // find diagonal entries
    NEWARRAY( diagPtr_, UInt, this->nrows_ );
    FindDiagonalEntries();
  }

  template<typename T>
  void VBR_Matrix<T>::CopySparsityPattern( StdMatrix &mat ) const {

    EXCEPTION( "VBR_Matrix::CopySparsityPattern: Not implemented" );

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
      UInt rStart; // start index of row
      UInt cStart; // start index of col
      UInt rbs;    // row block size
      UInt cbs;    // col block size
      UInt colNum; // column number
      UInt ind;    //index to data array

      // initialize return vector
      //rvec = mvec;

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
          for( UInt i = rStart; i < rStart+rbs; ++i ) {
            for( UInt j = cStart; j <cStart+cbs; ++j ) {
              rvec[i] += data_[ind++] * mvec[j];
            }  // loop over cols within block
          } // loop over rows within block
        } //block cols
      } // block rows
    }


  template<typename T>
  inline void VBR_Matrix<T>::MultSub( const Vector<T> &mvec,
                                      Vector<T> &rvec ) const {

    EXCEPTION("Implement me");
  }

  template<typename T>
  inline void VBR_Matrix<T>::CompRes( Vector<T> &r, const Vector<T> &x,
                                      const Vector<T> &b ) const {
    UInt rStart; // start index of row
    UInt cStart; // start index of col
    UInt rbs;    // row block size
    UInt cbs;    // col block size
    UInt colNum; // column number
    UInt ind;    //index to data array

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
        for( UInt i = rStart; i < rStart+rbs; ++i ) {
          for( UInt j = cStart; j <cStart+cbs; ++j ) {
            r[i] -= data_[ind++] * x[j];
          }  // loop over cols within block
        } // loop over rows within block
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
    EXCEPTION("Implement me");
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
  void VBR_Matrix<T>::Export( const char *fname,
                              const char *comment ) const{
    // Open output file and check for errors
    FILE *fp = fopen( fname, "w" );
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
    if ( comment != NULL ) {
      fprintf( fp, "%%\n%% %s\n%%\n", comment );
    }
    else {
      fprintf( fp, "%%\n%% Matrix exported by OLAS\n%%\n" );
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
  void VBR_Matrix<T>::FindBlock(UInt row, UInt col, UInt& rowB, UInt& colB, UInt& offset) {
    UInt l = 0;
    UInt u = nbRows_-1; 
    UInt k;
    bool found = false;
    offset = 0;
    UInt colSize;
    
    UInt numIter = 0;
    
    // use initial guess for k: multiply row by avg. blockwidth to get good starting value
    k = UInt(oneOverBlockSize_ * row );
    
    // look for row
    while( l <= u ) {
      numIter++;
      
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
//    std::cerr << "sarch-iter I: " << numIter << std::endl;
    // look for col
    l = rowPtr_[rowB];
    u = rowPtr_[rowB+1]-1;
    while( l <= u ) {
      numIter++;
      k = (l+u) >> 1; // half the interval
      if( bCol_[colInd_[k]] <=  col && bCol_[colInd_[k]+1] > col) {
        found = true;
        colB = colInd_[k];
        colSize = bCol_[colInd_[k]+1] - bCol_[colInd_[k]];
        offset = (col - bCol_[colInd_[k]])+offset*colSize; 
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
    
//   std::cerr << "sarch-iter: " << numIter << std::endl;
  }

  
  template<typename T>
  void VBR_Matrix<T>::AddToMatrixEntry( UInt i, UInt j, const T& v ){

    
    // loop over bRow_ until block row is found (log search)
    UInt rowB, colB, offset;
//    std::cerr << "\n\n===============================\n";
    FindBlock(i,j,rowB,colB, offset);
//    std::cerr << "rowB of row " << i << " is " << rowB << std::endl;
//    std::cerr << "colB of col " << j << " is " << colB << std::endl;
//    std::cerr << "offset = " << offset << std::endl;
//    
    data_[offset] += v;
    // loop over bCol until block col is found (log search)
   
    // calculate diff blockRow_[r]-i and blockCol_[c]-j to find entry
    
    // set entry in data_ array
    
    
  }


  template<typename T>
  void VBR_Matrix<T>::SetMatrixEntry( UInt i, UInt j, const T &v,
                                      bool setCounterPart ) {
    
    EXCEPTION("Implement me");
  }
  
  template<typename T>
  void VBR_Matrix<T>::GetMatrixEntry( UInt i, UInt j,
                                      T &v ) const {

    EXCEPTION("Implement me");
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
  
  template<typename T>
   T* VBR_Matrix<T>::GetDiagBlock( UInt blockRow, UInt& size, UInt& rowStart ) {
   
    T* ret = &data_[diagBlockPtr_[blockRow]];
    size = bRow_[blockRow+1] - bRow_[blockRow]; 
    rowStart = bRow_[blockRow];
    return ret;
  }
  
  

  template<typename T>
  void VBR_Matrix<T>::Scale( Double factor ) {
    
EXCEPTION("Not implemented");
  }


  template<typename T>
  Double VBR_Matrix<T>::GetMaxDiag() const {
    
    double maxDiag = 0;
    double current = 0;
    UInt i = 0;
    for( i = 0; i < this->nrows_; ++i ) {
      current = OpType<T>::MaxDiag( data_[ diagPtr_[i] ] );
      maxDiag = maxDiag > current ? maxDiag : current; 
    }
    
    return maxDiag;
  }


  template<typename T>
  void VBR_Matrix<T>::Add( const Double alpha, const StdMatrix& mat ) {
    EXCEPTION("Implement me");
  }

  template<>
  void VBR_Matrix<Complex>::Add( const Double alpha, const StdMatrix& mat ) {

    EXCEPTION("Implement me");
  }

  template <typename T>
  void VBR_Matrix<T>::FindDiagonalEntries() {

    // loop over all rows
    UInt bRow, bCol, offset;
    for( UInt i = 0; i < this->nrows_; ++i ) {
      FindBlock( i, i, bRow, bCol, offset );
      diagPtr_[i] = offset;
    }
  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class VBR_Matrix<Double>;
  template class VBR_Matrix<Complex>;
#endif

} // end of namespace

  // ==========================================================================
  //   TEST SECTION: VERIFY MATRIX IMPLEMENTATION
  // ==========================================================================
  // This test takes the reference matrix for the VBR-format from the guide
  //    "A reference model implementation of the Sparse BLAS in Fortran 95"
  //    http://dl.acm.org/citation.cfm?doid=567806.567811
  
#ifdef VBR_TEST_SECTION
#include "MatVec/crs_matrix.hh"
#include "MatVec/matrixLapackSupport.hh"
#include "MatVec/matrixBLASSupport.hh"

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
   
   WARN("Remove temporary inlcludes for matrix");
   
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
   
   // create rhs row vecrtor and compare
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
   
   
   // Get diagonal blocks and invert them
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
     
     // now try to invert block
//      // ------------------------
//      //defines if the upper (U) or lower (L) triangle of the solution is stored
//      char uplo;
//      int n, nrhs, lda, ldb, info;
//      //A is symmetric, so it makes no difference is the data are accessed row or columnwise
//      double * A = block;
//      //B is not necessary symmetric. so you have to give the function the transposed matrix
//      double * B = inv;
//
//      //set up parameters according to the netlib reference
//      //we want to store the lower triangle
//      uplo = 'U';
//      //n is the dimension of A
//      n = size;;
//      //nrhs is th enumber of colums of B and therefore the number of rows of B_trans
//      nrhs = size;;
//      //lda is the leading dimension of A -> max(1,n);
//      lda = n;
//      //ldb is the leading dimension of B -> max(1,n);
//      ldb = n;
//      StdVector<Integer> piv(size);
//      //solve system with dposv_
//      //dposv_(&uplo, &n , &nrhs, A,&lda, B, &ldb, &info);
//      dgesv_(&n , &nrhs, A,&lda, &piv[0], B, &ldb, &info);
//      std::cerr << "info: " << info << std::endl;
//      
     
     // ========================
     //  INVERT GENERAL MATRIX
     // ========================
     double * fac = new double[size*size];
     std::copy(block,block+size*size, fac);
     int N = size;
     int *IPIV = new int[N+1];
     int LWORK = N*N;
     double *WORK = new double[LWORK];
     int INFO;

     LP_DGETRF(&N,&N,fac,&N,IPIV,&INFO);
     LP_DGETRI(&N,fac,&N,IPIV,WORK,&LWORK,&INFO);

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
   
   
   // ========================
   //  TEST BLAS DGEMMV
   // ========================
   
   // test 
   {
     Matrix<Double> a(2,2);
     a[0][0] = 1.0;
     a[0][1] = 2.0;
     a[1][0] = 3.0;
     a[1][1] = 4.0;
     Vector<Double> r(4), s(4);
     s.Init(1.0);
     r.Init(0.0);
     
     char trans = 'T';
     int size = 2;
     double alpha = 1.0;
     double beta = 0.0;
     int inc = 1;
     
     DGEMV(&trans, 
           &size,
           &size,
           &alpha,
           a[0],
           &size,
           &s[1], 
           &inc,
           &beta,  &r[1], &inc);
     std::cerr << "r = " << r.ToString() << std::endl;

   }
  
 }
}

#endif
