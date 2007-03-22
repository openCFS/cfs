// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "precond/ilu0precond.hh"

namespace OLAS {


  // =====================================================
  //   Constructor (for use in GenerateStdPrecondObject)
  // =====================================================
  template <typename T>
  ILU0Precond<T>::ILU0Precond( const StdMatrix& mat, OLAS_Params *myParams,
			       OLAS_Report *myReport ) {

    ENTER_FCN( "ILU0Precond::ILU0Precond" );
    
    // Set pointers to communication objects
    this->myParams_ = myParams;
    this->myReport_ = myReport;
   
    // Set size information
    size_ = mat.GetNrows();

    // Initialise internal data arrays
    NewArray( diagPos_, UInt, size_ );
    NewArray( ilu_rptr_, UInt, size_ + 1 );
    ilu_cidx_ = NULL;
    ilu_data_ = NULL;

  }


  // ==============
  //   Destructor
  // ==============
  template <typename T>
  ILU0Precond<T>::~ILU0Precond() {

    ENTER_FCN( "ILU0Precond::~ILU0Precond" );

    // De-allocate internal data arrays
    DeleteArray( diagPos_  );
    DeleteArray( ilu_rptr_ );
    DeleteArray( ilu_cidx_ );
    DeleteArray( ilu_data_ );

  }


  // =================================
  //   Application of Preconditioner
  // =================================
  template <typename T>
  void ILU0Precond<T>::Apply( const CRS_Matrix<T> &mat,
			     const Vector<T> & f, Vector<T> & u ) const {

    ENTER_IFCN( "ILU0Precond::Apply" );

    UInt i, j;
    T_Mtype sum;

    // Step 1: Solve L * u~ = f by forward substitution
    for ( i = 1; i <= size_; i++ ) {
      sum = 0.0;
      for ( j = ilu_rptr_[i]; j < diagPos_[i]; j++ ) {
	sum += ilu_data_[j] * u[ilu_cidx_[j]];
      }
      u[i] = f[i] - sum;
    }
  
    // Step 2: Solve U * u = u~ by backward substitution
    for ( i = size_; i > 0; i-- ) {
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

    ENTER_IFCN( "ILU0Precond::Setup" );

    // Assure that matrix is in proper format for current implementation
    // NOTE: By default while using CFS++ the matrix will be in LEX format,
    // so it would be more efficient to re-write the code to take advantage
    // of that format, or simply to use a diagPtr_ for acessing diagonal
    // entries!
    //    mat.ChangeLayout( CRS_Matrix<T>::LEX_DIAG_FIRST );

    // Get the data of the matrix. Note that we use only CRS here
    const Integer *a_rptr = mat.GetRowPointer();
    const Integer *a_cidx = mat.GetColPointer();
    const T_Mtype *a_val  = mat.GetDataPointer();

    // How large is the matrix to be factored?
    UInt a_nnz  = mat.GetNnz();
    UInt a_size = mat.GetNrows();

    // Allocate memory for the ILU entries
    NewArray( ilu_data_, T_Mtype, a_nnz );
    NewArray( ilu_cidx_, UInt, a_nnz );

    //set ilu-matrix structure according to system matrx

    UInt row;
    // first we copy the row pointers since they remain the same
    for ( row = 1; row <= a_size + 1; row++ ) {
      ilu_rptr_[row] = a_rptr[row];
    }

    //set data and column indices
    for ( UInt i = 1; i <= a_nnz; i++) {
      ilu_data_[i] = a_val[i];
      ilu_cidx_[i] = a_cidx[i];
    }

    //set diagonal entries
    UInt istart, iend;
    UInt counter = 1;
    for ( UInt row = 1; row <= a_size; row++) {
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
 
    Integer i,j,k, j1, j2, jj, jw, jrow;
    T_Mtype  tl;

    //get help array and set it to zero;
    Integer *help;
    NewArray(help,Integer,a_size+1);
  
    for (i=1; i<=a_size+1; i++)
      help[i] = 0;

    // ILU: see Saad book: pp. 296
    // some things might have changed while porting to OLAS
    // since I don't have the book
  
    Integer percentDone = 0;
    Double  actDone;

    // k walks over all rows
    for (k=1; k<=a_size; k++) {
      // Keep user informed on progress
      actDone = (double)(k*100) / (double)size_;
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
      for (j=j1; j<diagPos_[k]; j++) {
         
	// jrow keeps track of the column index of the
	// column at position j we are treating.
	// the diagonal element of row jrow is our pivot element
	jrow = ilu_cidx_[j];

	//tl stores Akj/piv
	tl = ilu_data_[j]/ilu_data_[diagPos_[jrow]];
         
	ilu_data_[j] = tl;
	
	// jj selects entries in the row of the pivot element in U
	for (jj=diagPos_[jrow]+1; jj<ilu_rptr_[jrow+1]; jj++) {
            
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
	Error("Zero pivot in ILU setup!",__FILE__,__LINE__);
      }// 0 pivot

      //set help-array to zero
      for ( j = j1; j < j2; j++ ) {
	help[ilu_cidx_[j]] = 0;
      }
    }//k

    (*cla) << "\n \n " << std::endl;

    DeleteArray(help);

  }// Setup for ILU


}//namespace

