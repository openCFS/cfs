// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "General/environment.hh"
#include "General/exception.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/opdefs.hh"
#include "MatVec/crs_matrix.hh" // IWYU pragma: keep
#include "MatVec/stdmatrix.hh"
#include "OLAS/precond/ilu0precond.hh"

namespace CoupledField {


// =====================================================
//   Constructor (for use in GenerateStdPrecondObject)
// =====================================================
template <typename T> class Vector;

template <typename T>
ILU0Precond<T>::ILU0Precond( const StdMatrix& mat, PtrParamNode solverNode,
                             PtrParamNode olasInfo ) {


  // Set pointers to communication objects
  this->xml_ = solverNode;
  this->olasInfo_ = olasInfo;

  // Set size information
  size_ = mat.GetNumRows();

  // Initialise internal data arrays
  NEWARRAY( diagPos_, UInt, size_ );
  NEWARRAY( ilu_rptr_, UInt, size_ + 1 );
  ilu_cidx_ = NULL;
  ilu_data_ = NULL;

}


// ==============
//   Destructor
// ==============
template <typename T>
ILU0Precond<T>::~ILU0Precond() {
  // De-allocate internal data arrays
  delete [] ( diagPos_  );
  delete [] ( ilu_rptr_ );
  delete [] ( ilu_cidx_ );
  delete [] ( ilu_data_ );
}


// =================================
//   Application of Preconditioner
// =================================
template <typename T>
void ILU0Precond<T>::Apply( const CRS_Matrix<T> &mat,
                            const Vector<T> & f, Vector<T> & u ) const {



  UInt i, j;
  T_Mtype sum;

  // Step 1: Solve L * u~ = f by forward substitution
  for ( i = 0; i < size_; i++ ) {
    sum = 0.0;
    for ( j = ilu_rptr_[i]; j < diagPos_[i]; j++ ) {
      sum += ilu_data_[j] * u[ilu_cidx_[j]];
    }
    u[i] = f[i] - sum;
  }

  // Step 2: Solve U * u = u~ by backward substitution
  i = size_;
  while ( i > 0 ) {
    i--;
    sum = 0.0;
    for ( j = diagPos_[i] + 1; j < ilu_rptr_[i+1]; j++ ) {
      sum += ilu_data_[j] * u[ilu_cidx_[j]];
    }
    u[i] = (u[i] - sum) / ilu_data_[diagPos_[i]];
  }

}

// ===========================
//   Setup of Preconditioner
// ===========================
template <typename T>
void ILU0Precond<T>::Setup( CRS_Matrix<T> &mat ) {

  // get correct ParamNode
  PtrParamNode pNode = this->xml_->Get("ILU0", ParamNode::INSERT );
  pNode->GetValue("logging", logging_, ParamNode::INSERT ) ;

  // Assure that matrix is in proper format for current implementation
  // NOTE: By default while using CFS++ the matrix will be in LEX format,
  // so it would be more efficient to re-write the code to take advantage
  // of that format, or simply to use a diagPtr_ for acessing diagonal
  // entries!
  //    mat.ChangeLayout( CRS_Matrix<T>::LEX_DIAG_FIRST );

  // Get the data of the matrix. Note that we use only CRS here
  const UInt *a_rptr = mat.GetRowPointer();
  const UInt *a_cidx = mat.GetColPointer();
  const T_Mtype *a_val  = mat.GetDataPointer();

  // How large is the matrix to be factored?
  UInt a_nnz  = mat.GetNnz();
  dim_ = mat.GetNumRows();
  UInt a_size = mat.GetNumRows();

  // Allocate memory for the ILU entries
  NEWARRAY( ilu_data_, T_Mtype, a_nnz );
  NEWARRAY( ilu_cidx_, UInt, a_nnz );

  //set ilu-matrix structure according to system matrx

  UInt row;
  // first we copy the row pointers since they remain the same
  for ( row = 0; row < a_size + 1; row++ ) {
    ilu_rptr_[row] = a_rptr[row];
  }

  //set data and column indices
  for ( UInt i = 0; i < a_nnz; i++) {
    ilu_data_[i] = a_val[i];
    ilu_cidx_[i] = a_cidx[i];
  }

  //set diagonal entries
  UInt istart, iend;
  UInt counter = 0;
  for ( UInt row = 0; row < a_size; row++) {
    istart = ilu_rptr_[row];
    iend   = ilu_rptr_[row+1];
    for ( UInt k = istart; k < iend; k++) {
      if ( a_cidx[k] == row ) {
        diagPos_[row] = counter;
      }
      counter++;
    }
  }


  // now perform the actual factorization

  UInt i,j,k, j1, j2, jj, jw, jrow;
  T_Mtype  tl;

  //get help array and set it to zero;
  Integer *help;
  NEWARRAY(help,Integer,a_size+1);

  for (i=0; i<a_size+1; i++)
    help[i] = 0;

  // ILU: see Saad book: pp. 296
  // some things might have changed while porting to OLAS
  // since I don't have the book

  Integer percentDone = 0;
  Double  actDone;

  // k walks over all rows
  for (k=0; k<a_size; k++) {
    // Keep user informed on progress
    actDone = (double)(k*100) / (double)(size_-1);
    actDone = (UInt)(actDone/10.0)*10;
    if ( actDone > percentDone ) {
      percentDone = (UInt)actDone;
      (*cla) << " .. " << percentDone << "%" << std::flush;
    }

    // set help array for nonzero column indices of row k 
    // the help array contains the positions of nonzero entries 
    // in the modified CRS format
    j1 = ilu_rptr_[k];
    j2 = ilu_rptr_[k+1];
    for (j=j1; j<j2; j++)
      help[ilu_cidx_[j]] = j;

    //we now simultaneously treat L and U:

    // walk over all nonzero values before diagonal
    for (j=j1; j< diagPos_[k]; j++) {

      // jrow keeps track of the column index of the
      // column at position j we are treating.
      // the diagonal element of row jrow is our pivot element
      jrow = ilu_cidx_[j];

      //tl stores Akj/piv
      tl = ilu_data_[j]/ilu_data_[diagPos_[jrow]];

      ilu_data_[j] = tl;

      // jj selects entries in the row of the pivot element in U
      for (jj=diagPos_[jrow]+1; jj < ilu_rptr_[jrow+1]; jj++) {

        // find out if there is an entry in row k under jj
        // and if so, where it is
        jw = help[ilu_cidx_[jj]];

        if (jw!=0){
          ilu_data_[jw] -= tl*ilu_data_[jj];
        }

      }//for jj

      // if we have finished the row we stop the j loop
      // obviously there is no diagonal element and the 
      // algorithm fails
      if (j >= j2) {
        (*cla) << " terminating ILU setup. k=" << k << "jrow: " << jrow << 
            std::endl;
        break;    
      }

    }//for j (lower)

    jrow = ilu_cidx_[j];

    // if the diagonal has not been found or is zero
    // we terminate Setup and state that it failed.
    if (jrow!=k || ilu_data_[j] == 0.0) {
      (*cla) << "Zero pivot in ILU setup: row: " << jrow << " k=" << k
          << " val=" << ilu_data_[j] << std::endl;
      EXCEPTION("Zero pivot in ILU setup!");
    }// 0 pivot

    //set help-array to zero
    for ( j = j1; j < j2; j++ ) {
      help[ilu_cidx_[j]] = 0;
    }
  }//k

  (*cla) << "\n \n " << std::endl;

  delete [] (help); help = NULL;

  // If the user wishes, we can export the LU factorisation to a file
  std::string saveFacFile = "ilu0_fac.out";
  if(pNode->Has("saveFacFile")) {
    pNode->GetValue("saveFacFile", saveFacFile, ParamNode::INSERT);

    this->ExportILUFactorisation( saveFacFile);
    if ( logging_ == true ) {
      (*cla) << " Exported factor matrix to file '" << saveFacFile << "'"
          << std::endl;
    }
  }
}// Setup for ILU

// ============================
//   Export ILU factorisation
// ============================
template <typename T>
void ILU0Precond<T>::ExportILUFactorisation( const std::string& fileName ) {
  // ********************
  //   Open output file
  // ********************
  FILE *fp = fopen( fileName.c_str(), "w" );
  if ( fp == NULL ) {
    EXCEPTION( "CroutLU::ExportILUFactorisation: Cannot open file "
        << fileName << " for writing!" );
  }

  // *********************
  //   Write file header
  // *********************

  // Matrix Market Format Specification
  if ( EntryType<T>::M_EntryType == BaseMatrix::DOUBLE ) {
    fprintf( fp, "%%%%MatrixMarket matrix coordinate real general\n" );
  }
  else if ( EntryType<T>::M_EntryType == BaseMatrix::COMPLEX ) {
    fprintf( fp, "%%%%MatrixMarket matrix coordinate complex general\n" );
  }
  else {
    EXCEPTION( "ILU0Precond::ExportILUFactorisation: Can't identify "
        << "template parameter" );
  }

  // Comment
  fprintf( fp, "%%\n%% (I)LU0 factorisation matrix F = L + U - I \n" );

  // Information on number of rows, columns and entries
  fprintf( fp, "%d\t%d\t%d\n", dim_, dim_,
           ilu_rptr_[dim_-1]);

  // ******************************
  //   Write entries of ILU-array
  // ******************************
  UInt i, j, nblocks;

  // loop over all rows
  for ( i = 0; i < dim_; i++ ) {

    // get number of blocks in i-th row
    nblocks = ilu_rptr_[i+1] - ilu_rptr_[i];

    // loop over all blocks in this row
    for ( j = 0; j < nblocks; j++ ) {

      // store row and column index
      fprintf( fp, "%6d\t%6d\t", i  + 1,
               ( ilu_cidx_[ilu_rptr_[i]+j] )  + 1);

      // store non-zero entry
      OpType<T>::ExportEntry( this->ilu_data_[ilu_rptr_[i]+j],
                              0, 0, fp );
      fprintf( fp, "\n" );
    }
  }


  // *********************
  //   Close output file
  // *********************
  if ( fclose( fp ) == EOF ) {
    WARN("CroutLU::ExportILUFactorisation: Could not close file "
        << fileName << " after writing!");
  }
}

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class ILU0Precond<Double>;
template class ILU0Precond<Complex>;
#endif

}//namespace

