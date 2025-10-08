#ifndef OLAS_DIAG_MATRIX_HH
#define OLAS_DIAG_MATRIX_HH

#include <iostream>



#include "SparseOLASMatrix.hh"
#include "Vector.hh"
#include "CoordFormat.hh"

namespace CoupledField {

#if 0
  //forward declaration
  template <typename T>
  class PreMatrix;
#endif
  
  //! Class for representing matrices just storing the diagonal entries

  //! This class allows to store an diagonal matrix
  template<typename T>
  class Diag_Matrix : public SparseOLASMatrix<T> {

  public:
    using SparseOLASMatrix<T>::Mult;
    using SparseOLASMatrix<T>::MultT;
    using SparseOLASMatrix<T>::MultAdd;
    using SparseOLASMatrix<T>::MultTAdd;
    using SparseOLASMatrix<T>::MultSub;
    using SparseOLASMatrix<T>::CompRes;
    using SparseOLASMatrix<T>::Add;
    using SparseOLASMatrix<T>::SetSparsityPattern;
    using SparseOLASMatrix<T>::SetMatrixEntry;
    using SparseOLASMatrix<T>::GetMatrixEntry;
    using SparseOLASMatrix<T>::SetDiagEntry;
    using SparseOLASMatrix<T>::GetDiagEntry;
    using SparseOLASMatrix<T>::AddToMatrixEntry;

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION AND SETUP
    // =======================================================================

    //! \name Construction, Destruction and Setup
    //! The following methods are related to construction and destruction of
    //! Diag_Matrix objects, to initialisation of matrix values and to the
    //! setting and manipulation of the sparsity pattern of the matrix.
    //@{

    //! Default Constructor
    Diag_Matrix() {
      data_             = NULL;
      this->nnz_        = 0;
      this->ncols_      = 0;
      this->nrows_      = 0;
    }

    //! Copy Constructor

    //! This is a deep copy constructor. It is mainly provided to allow for
    //! certain sparse arithmetic operations of the SBM_Matrix class.
    //! \note  The constructed matrix is a copy and thus has the same layout
    //!        as the original matrix it is copied from.
    Diag_Matrix( const Diag_Matrix<T> &origMat ) {
      EXCEPTION("Diag_Matrix do not support a copy form a matrix");
    }

#if USE_MULTIGRID    
    //! Constructor using a PreMatrix
    
    //! If a matrix has been created not from a graph, but by successively
    //! inserting non-zero entries into a PreMatrix object, this method can be
    //! used to convert that PreMatrix object to a matrix in Diag format. The
    //! PreMatrix object itself remains unchanged by this operation.
    //! \param nrows number of rows
    //! \param ncols number of columns
    //! \param premat the pre-matrix with the matrix data
    Diag_Matrix( UInt nrows, UInt ncols, const PreMatrix<T>& premat ) {
      EXCEPTION("Diag_Matrix do not support a construction from a prematrrix");
    }
#endif
    
    //! Constructor with memory allocation

    //! This constructor generates a matrix of given dimensions and allocates
    //! space to store the specified number of non-zero entries.
    //! Initially all matrix entries are set to zero.
    //! \note  Since no column index information is provided the constructor
    //!        sets the currentLayout_ format to UNSORTED.
    //! \param nrows number of rows of the matrix
    //! \param ncols number of columns of the matrix
    //! \param nnz   number of non-zero matrix entries
    Diag_Matrix( UInt nrows, UInt ncols, UInt nnz ) {


      // assign basic properties
      this->nnz_   = nrows;
      this->nrows_ = nrows;
      this->ncols_ = nrows;

      // Allocate memory for matrix entries and initialise to zero
      NEWARRAY( data_, T, this->nnz_ );
      Init();
    }

    //! Constructor employing a CoordFormat matrix (not allowed)
    Diag_Matrix( CoordFormat<T> &sparseMat ) {
      EXCEPTION("Diag_Matrix do not support a copy form a CoordFormat matrix");
    }

    //! Destructor

    //! The default destructor is deep. It frees all memory dynamically
    //! allocated for attributes of this class.
    ~Diag_Matrix() {
      delete [] ( data_ );
    }

    //! Setup the sparsity pattern of the matrix (not necessary)
    void SetSparsityPattern( BaseGraph &graph ) {
      //nothing to do
    };


    //! Transfer ownership of sparsity pattern from matrix to pool
    PatternIdType TransferPatternToPool( PatternPool *patternPool ) {
      // nothing to dp
      return NO_PATTERN_ID;
    }

    //! Copy the sparsity pattern of another matrix to this matrix 
    
    //! This method copies the sparsity pattern of the source matrix to the
    //! current matrix
    void SetSparsityPattern ( const StdMatrix &srcMat ) {
      EXCEPTION("Diag_Matrix do not support a SetSparsityPattern");
    };
    
    //! Re-set internal dimensions

    //! This method can be used to re-size the matrix. If the internal
    //! dimensions differ from the ones previously set, then the internal
    //! data arrays will be thrown away and re-allocated with the new sizes.
    void SetSize( UInt nrows, UInt ncols, UInt nnz );
 
    //! Initialize matrix to zero
    
    //! This matrix sets all matrix entries to zero. Note that all entries
    //! in the non-zero structure of the matrix get a zero value. However,
    //! the structure/data layout remains the same. All positions in the
    //! structure can still be over-written with non-zero values.
    inline void Init() {
      for ( UInt i = 0; i < this->nnz_; i++ ) {
        data_[i] = 0;
      }
    }

    // =======================================================================
    // QUERY / MANIPULATE MATRIX FORMAT
    // =======================================================================

    //! \name Query matrix / manipulate format
    //@{

    //! Return the storage type of the matrix
  
    //! The method returns the storage type of the matrix. The latter is
    //! encoded as a value of the enumeration data type MatrixStorageType.
    //! \return DIAG
    inline BaseMatrix::StorageType GetStorageType() const {
      return BaseMatrix::DIAG;
    }
    
    //! \copydoc BaseMatrix::GetMemoryUsage()
    Double GetMemoryUsage() const;
    //@}


    // =======================================================================
    // OBTAIN / MANIPULATE MATRIX ENTRIES
    // =======================================================================

    //! \name Obtain/manipulate matrix entries
    //@{

    //! Set a specific matrix entry
    
    //! This functions sets the Matrix entry in row i and column j
    //! to the value v. Depending on the entry type of the matrix,
    //! v might be a Double, Complex or a tiny Matrix of either type.
    //! If setCounterPart is set, also the transposed entry is set (only
    //! if the position is non-diagonal).
    void SetMatrixEntry( UInt i, UInt j, const T &v, 
                         bool setCounterPart ) {
      EXCEPTION("Diag_Matrix do not support SetMatrixEntry");
    };
    
    //! Get a specific matrix entry
    
    //! This function returns a reference to the Matrix entry in 
    //!  row i and column j. Depending on the entry type of the matrix,
    //!  v might be a Double, Complex or a tiny Matrix of either type
    void GetMatrixEntry( UInt i, UInt j, T &v ) const {
      EXCEPTION("Diag_Matrix do not support GetMatrixEntry");
    };

    //! Set the diagonal entry of row i to the value of v
    void SetDiagEntry( UInt i, const T &v ) {
      data_[i] = v;
    }

    //! Returns the value of the diagonal entry of row i
    void GetDiagEntry( UInt i, T &v ) const {
      v = data_[i];
    }

    T GetDiagEntry(unsigned int row) const { return data_[row]; }

    //! This routine adds the value of v to the matrix entry at (i,j)

    //! This routine can be used to add the value of the parameter v to the
    //! matrix entry belonging to the index pair (i,j). This method relies on
    //! the possibly slow "dense matrix style" access.  
    void AddToMatrixEntry( UInt i, UInt j, const T& v );

    //! Return the length (i.e. number of non-zero entries) of i-th row
    UInt GetRowSize(UInt i) const {
      return 1;
    }

    //! Set the length (i.e. number of non-zero entries) of i-th row

    //! The method can be used to set the length (i.e. number of non-zero
    //! entries) of row i to the value size.
    void SetRowSize( UInt i, UInt size ) {
      EXCEPTION("Diag_Matrix do not support SetRowSize");
    }

    //! Insert a value/column index at position (row,pos)
    void Insert( UInt row, UInt col, UInt pos ) {
      EXCEPTION("Diag_Matrix do not support Insert");
    }

    //@}


    // =======================================================================
    // MATRIX-VECTOR PRODUCTS
    // =======================================================================

    //! \name Matrix-Vector products 

    //! The following methods perform arithmetic operations of the form
    //! \f$ \vec{y} = \alpha \vec{y} + \beta A \vec{x}\f$ or similar

    //@{

    //! Perform a matrix-vector multiplication rvec = this*mvec

    //! This method performs a matrix-vector multiplication rvec = this*mvec.
    void Mult(const Vector<T> &mvec, Vector<T> &rvec) const;

    //! Perform a matrix-vector multiplication rvec = rvec + this*mvec

    //! This method performs a matrix-vector multiplication
    //! rvec = rvec + this*mvec.
    void MultAdd( const Vector<T> &mvec, Vector<T> &rvec ) const;

    //! Perform a matrix-vector multiplication rvec = rvec - this*mvec

    //! This method performs a matrix-vector multiplication
    //! rvec = rvec - this*mvec.
    void MultSub( const Vector<T> &mvec, Vector<T> &rvec ) const;

    //! Compute the residual \f$r=b-Ax\f$
    
    //! For a given a right-hand side b and an approximate solution vector x
    //! this function computes the residual vector \f$r=b-Ax\f$ of the
    //! corresponding linear system of equations.
    void CompRes( Vector<T>& r, const Vector<T>& x, const Vector<T>& b) const;

    //@}


    // =======================================================================
    // MATRIX-VECTOR-PRODUCT WITH TRANSPOSE
    // =======================================================================

    //! \name Matrix-vector multiplication with transpose of matrix

    //! The following methods perform arithmetic operations of the form
    //! \f$ \vec{y} = \alpha \vec{y} + \beta A^T \vec{x}\f$ or similar

    //@{

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object rvec = transpose(this)*mvec.
    void MultT(const Vector<T> & mvec, Vector<T> & rvec) const;

    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec

    //! This method performs a matrix-vector multiplication with the transpose
    //! of the matrix object followed by an addition:
    //! rvec += transpose(this)*mvec.
    void MultTAdd( const Vector<T> & mvec, Vector<T> & rvec ) const;

    //! Perform a matrix-vector multiplication rvec -= transpose(this)*mvec

     //! This method performs a matrix-vector multiplication with the transpose
     //! of the matrix object followed by an addition:
     //! rvec -= transpose(this)*mvec.
     void MultTSub( const Vector<T> & mvec, Vector<T> & rvec ) const;
    //@}


    // =======================================================================
    // MISCELLANEOUS ARITHMETIC OPERATIONS
    // =======================================================================

    //! \name Miscellaneous arithmetic matrix operations
    //@{

    //! Method to scale all matrix entries by a fixed factor
    void Scale( Double factor );
    
    //! \copydoc StdMatrix::Scale(Double, std::set<UInt>)
    void Scale( Double factor, 
                const std::set<UInt>& rowIndices,
                const std::set<UInt>& colIndices ); 

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix,
    //! \f$A = A + \alpha B\f$. In doing so the sparsity structure of the
    //! matrix mat is assumed to be identical to this matrix' structure.
    void Add( const Double alpha, const StdMatrix& mat );
    void Add( const Complex alpha, const StdMatrix& mat );
    
    //! \copydoc StdMatrix::Add(Double,StdMatrix,std::set<UInt>,std::set<UInt>)
    void Add( const Double a, const StdMatrix& mat,
              const std::set<UInt>& rowIndices,
              const std::set<UInt>& colIndices );
    void Add( const Complex a, const StdMatrix& mat,
              const std::set<UInt>& rowIndices,
              const std::set<UInt>& colIndices );
    //@}


    // =======================================================================
    // I/O operations
    // =======================================================================

    //! \name I/O operations
    //@{

    //! Print the matrix to an output stream
    
    //! Prints the matrix in indexed format to the output stream os. This
    //! method is primarily intended for debugging.
    virtual std::string ToString( char colSeparator = ' ',
                                      char rowSeparator = '\n' ) const;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. Since it is a Diag_Matrix the output
    //! format will be symmetric, co-ordinate, real or complex.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    //! \param comment string to be inserted into file header
    void ExportMatrixMarket( const std::string& fname, const std::string& comment = "" ) const {
      EXCEPTION("Diag_Matrix does not support ExportMatrixMarket().");
    }

    //@}


    // =======================================================================
    // DIRECT ACCESS TO Diag DATA STRUCTURE
    // =======================================================================

    //! \name Direct access to Diag data structure
    //! The following methods allow a direct manipulation of the internal
    //! data structures. Though this does violate the concept of data
    //! encapsulation, these methods are provided for the sake of efficiency.
    //! \note
    //! Since %OLAS is one-based, the internal data arrays have an offset of
    //! one. The pointers returned by the methods below include this offset,
    //! i.e. the actual arrays start at "pointer value" + 1.

    //@{

    //! Get row pointer array
    UInt* GetRowPointer() {
      EXCEPTION("Diag_Matrix do not support GetRowPointer");
      return NULL;
    }

    //! Get row pointer array (read only)
    const UInt* GetRowPointer() const {
      EXCEPTION("Diag_Matrix do not support GetRowPointer");
      return NULL;
    }

    //! Get diagonal pointer array
    UInt* GetDiagPointer() {
      EXCEPTION("Diag_Matrix do not support GetDiagPointer");
      return NULL;
    }

    //! Get diagonal pointer array (read only)
    const UInt* GetDiagPointer() const {
      EXCEPTION("Diag_Matrix do not support GetDiagPointer");
      return NULL;
    }

    //! Get array of column indices
    UInt* GetColPointer() {
      EXCEPTION("Diag_Matrix do not support GetColPointer");
      return NULL;
    }

    //! Get array of column indices (read only)
    const UInt* GetColPointer() const {
      EXCEPTION("Diag_Matrix do not support GetColPointer");
      return NULL;
    }

    //! Get array of matrix entries
    T* GetDataPointer(){
      return data_;
    }

    //! Get array of matrix entries (read only)
    const T* GetDataPointer() const {
      return data_;
    }

    //@}

  private:

    //! Array containing the (potentially) non-zero matrix entries
    T* data_;

  };

}

#endif
