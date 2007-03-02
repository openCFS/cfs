#ifndef OLAS_DIAG_MATRIX_HH
#define OLAS_DIAG_MATRIX_HH

#include <iostream>

#include "utils/utils.hh"
#include "matvec/sparseolasmatrix.hh"
#include "matvec/vector.hh"
#include "matvec/coordformat.hh"

namespace OLAS {

  //forward declaration
  template <typename T>
  class PreMatrix;

  //! Class for representing matrices just storing the diagonal entries

  //! This class allows to store an diagonal matrix
  template<typename T>
  class Diag_Matrix : public SparseOLASMatrix<T> {

  public:

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
      ENTER_FCN( "Diag_Matrix::Diag_Matrix" );
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
      Error("Diag_Matrix do not support a copy form a matrix",
	    __FILE__,__LINE__);
    }

    //! Constructor using a PreMatrix
    
    //! If a matrix has been created not from a graph, but by successively
    //! inserting non-zero entries into a PreMatrix object, this method can be
    //! used to convert that PreMatrix object to a matrix in Diag format. The
    //! PreMatrix object itself remains unchanged by this operation.
    //! \param nrows number of rows
    //! \param ncols number of columns
    //! \param premat the pre-matrix with the matrix data
    Diag_Matrix( Integer nrows, Integer ncols, const PreMatrix<T>& premat ) {
      Error("Diag_Matrix do not support a construction from a prematrrix",
	    __FILE__,__LINE__);
    }

    //! Constructor with memory allocation

    //! This constructor generates a matrix of given dimensions and allocates
    //! space to store the specified number of non-zero entries.
    //! Initially all matrix entries are set to zero.
    //! \note  Since no column index information is provided the constructor
    //!        sets the currentLayout_ format to UNSORTED.
    //! \param nrows number of rows of the matrix
    //! \param ncols number of columns of the matrix
    //! \param nnz   number of non-zero matrix entries
    Diag_Matrix( Integer nrows, Integer ncols, Integer nnz ) {

      ENTER_FCN( "Diag_Matrix::Diag_Matrix" );

      // assign basic properties
      this->nnz_   = nrows;
      this->nrows_ = nrows;
      this->ncols_ = nrows;

      // Allocate memory for matrix entries and initialise to zero
      NewArray( data_, T, this->nnz_ );
      Init();
    }

    //! Constructor employing a CoordFormat matrix (not allowed)
    Diag_Matrix( CoordFormat<T> &sparseMat ) {
      Error("Diag_Matrix do not support a copy form a CoordFormat matrix",
	    __FILE__,__LINE__);
    }

    //! Destructor

    //! The default destructor is deep. It frees all memory dynamically
    //! allocated for attributes of this class.
    ~Diag_Matrix() {
      ENTER_FCN( "Diag_Matrix::~Diag_Matrix" );
      DeleteArray( data_    );
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

    //! Copy the sparsity pattern of this matrix to another 
    
    //! This method copies the sparsity pattern of this matrix and
    //! sets it in another matrix
    void CopySparsityPattern ( StdMatrix &mat ) const {
      Error("Diag_Matrix do not support a CopySparsityPattern",
	    __FILE__,__LINE__);
    };
    
    //! Re-set internal dimensions

    //! This method can be used to re-size the matrix. If the internal
    //! dimensions differ from the ones previously set, then the internal
    //! data arrays will be thrown away and re-allocated with the new sizes.
    void SetSize( Integer nrows, Integer ncols, Integer nnz );
 
    //! Initialize matrix to zero
    
    //! This matrix sets all matrix entries to zero. Note that all entries
    //! in the non-zero structure of the matrix get a zero value. However,
    //! the structure/data layout remains the same. All positions in the
    //! structure can still be over-written with non-zero values.
    inline void Init() {
      ENTER_FCN( "Diag_Matrix::Init" );
      for ( Integer i = 1; i <= this->nnz_; i++ ) {
        data_[i] = 0;
      }
    }

    //! Method to force instantiation of all public member functions
    void InstantiatePublicMethods();

    //@}


    // =======================================================================
    // QUERY / MANIPULATE MATRIX FORMAT
    // =======================================================================

    //! \name Query matrix / manipulate format
    //@{

    //! Return the storage type of the matrix
  
    //! The method returns the storage type of the matrix. The latter is
    //! encoded as a value of the enumeration data type MatrixStorageType.
    //! \return DIAG
    inline MatrixStorageType GetStorageType() const {
      return DIAG;
    }


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
    void SetMatrixEntry( Integer i, Integer j, T &v, 
                         bool setCounterPart ) {
      Error("Diag_Matrix do not support SetMatrixEntry",
	    __FILE__,__LINE__);
    };
    
    //! Get a specific matrix entry
    
    //! This function returns a reference to the Matrix entry in 
    //!  row i and column j. Depending on the entry type of the matrix,
    //!  v might be a Double, Complex or a tiny Matrix of either type
    void GetMatrixEntry( Integer i, Integer j, T &v ) const {
      Error("Diag_Matrix do not support GetMatrixEntry",
	    __FILE__,__LINE__);
    };


    //! Return the diagonal entry of row i
    inline
    T& GetDiag( Integer i ) {
      ENTER_IFCN( "Diag_Matrix::GetDiag" );
      return data_[i];
    }

    //! Return the diagonal entry of row i (read only)
    inline
    const T& GetDiag( Integer i ) const {
      ENTER_IFCN( "Diag_Matrix::GetDiag" );
      return data_[i];
    }

    //! Set the diagonal entry of row i to the value of v
    void SetDiagEntry( Integer i, T &v ) {
      ENTER_IFCN( "Diag_Matrix::SetDiagEntry" );
      data_[i] = v;
    }

    //! Returns the value of the diagonal entry of row i
    void GetDiagEntry( Integer i, T &v ) const {
      ENTER_IFCN( "Diag_Matrix::GetDiagEntry" );
      v = data_[i];
    }

    //! Determine maximum absolute value of diagonal entries

    //! This method determines the the maximal absolute value (on the scalar)
    //! level of the entries on the main diagonal of the matrix.
    Double GetMaxDiag() const;

    //! This routine adds the value of v to the matrix entry at (i,j)

    //! This routine can be used to add the value of the parameter v to the
    //! matrix entry belonging to the index pair (i,j). This method relies on
    //! the possibly slow "dense matrix style" access.  
    void AddToMatrixEntry( Integer i, Integer j, T& v );

    //! Return the length (i.e. number of non-zero entries) of i-th row
    Integer GetRowSize(Integer i) const {
      return 1;
    }

    //! Set the length (i.e. number of non-zero entries) of i-th row

    //! The method can be used to set the length (i.e. number of non-zero
    //! entries) of row i to the value size.
    void SetRowSize( Integer i, Integer size ) {
      Error("Diag_Matrix do not support SetRowSize",
	    __FILE__,__LINE__);
    }

    //! Insert a value/column index at position (row,pos)
    void Insert( Integer row, Integer col, Integer pos ) {
      Error("Diag_Matrix do not support Insert",
	    __FILE__,__LINE__);
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

    //@}


    // =======================================================================
    // MISCELLANEOUS ARITHMETIC OPERATIONS
    // =======================================================================

    //! \name Miscellaneous arithmetic matrix operations
    //@{

    //! Method to scale all matrix entries by a fixed factor
    void Scale( Double factor );

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix,
    //! \f$A = A + \alpha B\f$. In doing so the sparsity structure of the
    //! matrix mat is assumed to be identical to this matrix' structure.
    void Add( const Double alpha, const StdMatrix& mat );
    
    //@}


    // =======================================================================
    // I/O operations
    // =======================================================================

    //! \name I/O operations
    //@{

    //! Print the matrix to an output stream
    
    //! Prints the matrix in indexed format to the output stream os. This
    //! method is primarily intended for debugging.
    void Print(std::ostream& os) const;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. Since it is a Diag_Matrix the output
    //! format will be symmetric, co-ordinate, real or complex.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    //! \param comment string to be inserted into file header
    void Export( const Char *fname, const Char *comment = NULL ) const {
      Error("Diag_Matrix do not support  Export",
	    __FILE__,__LINE__);
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
    Integer* GetRowPointer() {
      Error("Diag_Matrix do not support GetRowPointer",
	    __FILE__,__LINE__);
      return NULL;
    }

    //! Get row pointer array (read only)
    const Integer* GetRowPointer() const {
      Error("Diag_Matrix do not support GetRowPointer",
	    __FILE__,__LINE__);
      return NULL;
    }

    //! Get diagonal pointer array
    Integer* GetDiagPointer() {
      Error("Diag_Matrix do not support GetDiagPointer",
	    __FILE__,__LINE__);
      return NULL;
    }

    //! Get diagonal pointer array (read only)
    const Integer* GetDiagPointer() const {
      Error("Diag_Matrix do not support GetDiagPointer",
	    __FILE__,__LINE__);
      return NULL;
    }

    //! Get array of column indices
    Integer* GetColPointer() {
      Error("Diag_Matrix do not support GetColPointer",
	    __FILE__,__LINE__);
      return NULL;
    }

    //! Get array of column indices (read only)
    const Integer* GetColPointer() const {
      Error("Diag_Matrix do not support GetColPointer",
	    __FILE__,__LINE__);
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
