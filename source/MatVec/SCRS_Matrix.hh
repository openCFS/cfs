#ifndef OLAS_SCRS_MATRIX_HH
#define OLAS_SCRS_MATRIX_HH

#include <iostream>



#include "SparseOLASMatrix.hh"
#include "Vector.hh"
#include "PatternPool.hh"

namespace CoupledField {

  template <typename T> class CoordFormat;
  
  //! Class for representing matrices in the Symmetric Compressed Row Format

  //! This class allows to store a symmetric sparse matrix in the
  //! <em>Symmetric Compressed Row Storage Format</em> (SCRS).
  //! In this format only the upper
  //! triangular part of the matrix is stored. This is accomplished in the
  //! following fashion. For every row of the matrix we take the part of the
  //! row right of the main diagonal including the diagonal entry and store
  //! from it all non-zero entries in a contiguous fashion. The entries are
  //! ordered lexicographically starting from the diagonal entry. We make use
  //! of one array for storing the index of the non-zero entry and one for its
  //! value. A third array is used to store the array indices that mark the
  //! start of a new matrix row in the two previous arrays.
  //! Like all %OLAS matrix classes also this one uses one-based indexing.
  //! \note
  //! - This class implicitely assumes that the diagonal matrix entries
  //!   are all non-zero, see the description of SetSparsityPattern below.
  //! - Debugging for this class can be turned on by defining the macro
  //!   DEBUG_SCRS_MATRIX
  template<typename T>
  class SCRS_Matrix : public SparseOLASMatrix<T> {


  public:
    using SparseOLASMatrix<T>::Mult;
    using SparseOLASMatrix<T>::MultT;
    using SparseOLASMatrix<T>::MultAdd;
    using SparseOLASMatrix<T>::MultTAdd;
    using SparseOLASMatrix<T>::MultSub;
    using SparseOLASMatrix<T>::MultTSub;
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
    //! SCRS_Matrix objects, to initialisation of matrix values and to the
    //! setting and manipulation of the sparsity pattern of the matrix.
    //@{

    //! Default Constructor
    SCRS_Matrix() {
      colInd_      = NULL;
      rowPtr_      = NULL;
      data_        = NULL;
      diagPtr_    = NULL;
      patternPool_ = NULL;
      patternID_   = NO_PATTERN_ID;
      this->nnz_   = 0;
      numEntries_  = 0;
      this->nrows_ = 0;
      this->ncols_ = 0;
    }

    //! Copy Constructor

    //! This is a deep copy constructor. It is mainly provided to allow for
    //! certain sparse arithmetic operations of the SBM_Matrix class.
    SCRS_Matrix( const SCRS_Matrix<T> &origMat );

    //! Constructor with memory allocation

    //! This constructor generates a matrix of given dimensions and allocates
    //! space to store the specified number of non-zero entries. The symmetry
    //! is taken into account, when determining the storage requirements.
    //! Initially all matrix entries are set to zero.
    //! \note We do not allocate memory to store the sparsity pattern here.
    //!       This is taken care of by the SetSparsityPattern() method.
    //! \param nrows number of rows of the matrix
    //! \param ncols number of columns of the matrix
    //! \param nnz   number of non-zero matrix entries. \b Note that this must
    //!              be the total number of non-zero entries, neglecting that
    //!              the matrix is symmetric and only (nnz + nrows) / 2 will
    //!              actually be stored. Therefore here the number of non-zero
    //!              must be passed, as if the matrix was not symmetric.
    SCRS_Matrix( Integer nrows, UInt ncols, UInt nnz ) : diagPtr_( NULL ){


      // assign storage type for dynamic type information
      // this->storagetype_ = SPARSE_SYM;
      this->nnz_ = nnz;
      this->nrows_ = nrows;
      this->ncols_ = ncols;
 
      // NOTE: This is correct iff there are no zeros on the main diagonal!
      numEntries_ = (this->nnz_ - this->nrows_) / 2 + this->nrows_;

      // Allocate memory for matrix entries and initialise to zero
      NEWARRAY( data_, T, numEntries_ );
      Init();

      // No memory allocation for the pattern at this point
      // We either do this in SetSparsityPattern() or use
      // a pattern from the pool
      colInd_ = NULL;
      rowPtr_ = NULL;
      diagPtr_ = NULL;

      // We do not know a pool or pattern id yet
      patternPool_ = NULL;
      patternID_ = NO_PATTERN_ID;
    }

    //! Alternate constructor

    //! This constructor variant is a sort of copy constructor. It expects as
    //! input a sparse matrix in coordinate format. It uses this information
    //! to generate and populate an equivalent SCRS_Matrix.
    //! \note
    //! - will only work for real, scalar entries currently (RMatrixDof1)
    //! - Currently only square matrices are supported.
    //! \param sparseMat sparse matrix in coordinate format
    SCRS_Matrix(CoordFormat<T> &sparseMat);

    SCRS_Matrix(UInt nr, UInt nc, UInt nnz, UInt* row_ptr, UInt* col_ptr, T* data_ptr);

    //! Destructor

    //! The default destructor is deep. It frees all memory dynamically
    //! allocated for attributes of this class.
    ~SCRS_Matrix() {


      // Delete data vector
      delete[] (data_); 

      // Check if we own the pattern or whether it belongs to a pool.
      // In the former case we must free the memory. In the latter we
      // must deregister
      if ( patternPool_ == NULL ) {
        delete[] (rowPtr_);
        delete[] (colInd_);
        delete[] (diagPtr_);
      }
      else {
        patternPool_->DeRegisterUser( patternID_ );
        patternPool_ = NULL;
        patternID_ = NO_PATTERN_ID;
      }

      // It is not our task to destroy the pattern pool object!
      patternPool_ = NULL;
    }

    //! Setup the sparsity pattern of the matrix
  
    //! This method uses the graph to generate the internal information needed
    //! for representing the sparsity pattern of the matrix correctly. Before
    //! this function has not been called the matrix is not ready to be used.
    //! \note
    //! The implementation of this routine expects the graph to be complete,
    //! i.e. it stores all neighbour hoodrelations. It does not expect,
    //! however, that the neighbours for a node are ordered in the array
    //! obtained with GetGraphRow(). But, it is also expected that the array
    //! contains an entry for the node itself.
    void SetSparsityPattern( BaseGraph &graph );

    //! The sparsity-pattern with data, which was already set outside this class

    //! We provide the three CRS-vectors externally as StdVectors and
    //! convert them to arrays in this method
    //! \param rowP  row-pointer
    //! \param colI  column-index
    //! \param data  data-pointer
    void SetSparsityPatternData( const StdVector<UInt>& rowP,
                                 const StdVector<UInt>& colI,
                                 const StdVector<T>& data);

    //! Setup the sparsity pattern of the matrix
  
    //! This method provides an alternative form to set the sparsity pattern
    //! of the matrix. If this method is used a fitting sparsity pattern must
    //! already have been generated, e.g. by another instance of this class
    //! that has the same pattern, and been added to the specified pool. In
    //! this case we can simply obtain copy the address information from the
    //! pool object.
    //! \param patternPool  pointer to a PatternPool object that stores the
    //!                     desired sparsity pattern
    //! \param patternID    identifier required to obtain the sparsity
    //!                     pattern from the pool of patterns
    void SetSparsityPattern( PatternPool *patternPool,
                             PatternIdType patternID );

    //! Transfer ownership of sparsity pattern from matrix to pool

    //! This method can be used to instruct the matrix to transfer the
    //! ownership of its sparsity pattern to a PatternPool object such
    //! that the sparsity pattern can be re-used by other instances of
    //! the class. After a call to the method the matrix object will no
    //! longer attempt to de-allocate memory for the pattern in its
    //! destructor, but instead only de-register itself from the pool.
    //! \param patternPool pointer to the PatternPool object to which the
    //!        ownership of the sparsity pattern is to be transfered.
    //! \return An identifier that can be used to identify the matrix
    //!         pattern when communicating with the PatternPool object.
    PatternIdType TransferPatternToPool( PatternPool *patternPool );

    //! Setup the sparsity pattern of the matrix (manual definition)
  
    //! This method provides an alternative form to set the sparsity pattern
    //! of the matrix. If we already know the pointers usually generated via
    //! SetSparsityPattern(), we can also set it directly.
    //! \param rowPtr       Array containing starting indices of the 
    //!                     different matrix rows
    //! \param colPtr       Array containing the column indices of the 
    //!                     non-zero matrix entries
    //! \param diagPtr      //! Array containing the indices of the 
    //!                     diagonal matrix entries
    void SetSparsityPattern( UInt *rowPtr, UInt *colPtr, UInt *diagPtr );

    //! Copy the sparsity pattern of another matrix to this matrix 
    
    //! This method copies the sparsity pattern of the source matrix to the
    //! current matrix
    void SetSparsityPattern( const StdMatrix &srcMat );
      

    //! Re-set internal dimensions

    //! This method can be used to re-size the matrix. If the internal
    //! dimensions differ from the ones previously set, then the internal
    //! data arrays will be thrown away and re-allocated with the new sizes.
    //! In this case all matrix entries and sparsity pattern information is
    //! lost, of course.
    void SetSize( UInt nrows, UInt ncols, UInt nnz ) {


      this->ncols_ = ncols;
      if ( this->nrows_ != nrows ) {
        this->nrows_ = nrows; 
        delete[] (rowPtr_);
        
        NEWARRAY( rowPtr_, UInt, this->nrows_ + 1 );
        rowPtr_[0] = 0;
      }
      this->nnz_ = nnz;

      //this is correct iff there are no 0s on the diagonal!
      if ( numEntries_!= (this->nnz_ + this->nrows_) / 2 ) {
        delete[] (data_);
        delete[] (colInd_);
        numEntries_ = (this->nnz_ + this->nrows_) / 2;
        NEWARRAY( this->data_, T, numEntries_ );
        NEWARRAY( this->colInd_, UInt, numEntries_ );
      }
    }
 
    //! Initialize matrix to zero
    
    //! This matrix sets all matrix entries to zero. Note that all entries
    //! in the non-zero structure of the matrix get a zero value. However,
    //! the structure/data layout remains the same. All positions in the
    //! structure can still be over-written with non-zero values.
    inline void Init() {
      for ( UInt i = 0; i< numEntries_; i++ ) {
        this->data_[i] = 0;
      }
    }

    //! Return the length of row i omitting all entries before the diagonal
    UInt GetRowSize(UInt i) const {
      return rowPtr_[i+1] - rowPtr_[i];
    }

    /** Return the maximal row size.
     * @see GetRowSize() */
    unsigned int GetMaxRowSize() const {
      unsigned int max = 0;
      for(unsigned int i = 0; i < this->nrows_; i++)
        max = std::max(GetRowSize(i), max);
      return max;
    }


    //! Returns the number of entries
    UInt GetNumEntries() const {
      return numEntries_;
    }

    //! Sets the number of entries
    void SetNumEntries(UInt numEntries) {
      numEntries_ = numEntries;
    }

    //! Set the length of row i.

    //! The method can be used to set the length of row i to the value size.
    //! Size should be chosen such that it includes only entries with a
    //! column index >= i, i.e. right of and including the main diagonal.
    void SetRowSize( UInt i, UInt size ) {
      rowPtr_[i+1] = rowPtr_[i] + size;
    }

    //! Insert a value/column index at position (row,pos)
    void Insert( UInt row, UInt col, UInt pos );

    //@}


    // =======================================================================
    // QUERY / MANIPULATE MATRIX FORMAT
    // =======================================================================

    //! \name Query matrix format and pattern
    //@{

    //! Return the storage type of the matrix
  
    //! The method returns the storage type of the matrix. The latter is
    //! encoded as a value of the enumeration data type MatrixStorageType.
    inline BaseMatrix::StorageType GetStorageType() const {
      return BaseMatrix::SPARSE_SYM;
    }

    //! Get pointer to sparsity pattern pool 

    //! The method returns the pointer to the pattern pool, from which the
    //! matrix got its sparsity pattern. If the matrix was created by using
    //! the graph object, the return pointer will be NULL.
    inline PatternPool* GetPatternPool() const {
      return patternPool_;
    }

    //! Get id of the sparsity pattern of the pattern pool

    //! This method return the id for the sparsity pattern which is used
    //! for communication with the pattern pool. If no pattern pool was used
    //! for the creation of the matrix, NO_PATTERN_ID is returned.
    inline PatternIdType GetPatternId() const {
      return patternID_;
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
    //! v might be a Double, Complex or a tiny Matrix of either type
    //! \note If setCounterPart is false and i != j, an error is thrown, 
    //!       since in this case the symmetry of the matrix can no longer be 
    //!       maintained!
    void SetMatrixEntry( UInt i, UInt j, const T &v,
                         bool setCounterPart );
    
    //! Get a specific matrix entry
    
    //! This function returns a reference to the Matrix entry in 
    //!  row i and column j. Depending on the entry type of the matrix,
    //!  v might be a Double, Complex or a tiny Matrix of either type
    void GetMatrixEntry( UInt i, UInt j, T &v ) const;

    //! Returns the Indices of column, row and the data as a vector.
    //! Should be the same result as in ExportMatrixMarket().

    void GetVectors( StdVector<UInt> *Vec_col, StdVector<UInt> *Vec_row, StdVector<T> *Vec_val ) const;

    //! Set the diagonal entry of row i to the value of v
    virtual void SetDiagEntry( UInt i, const T &v ) {
      data_[rowPtr_[i]] = v;
    }

    //! Returns the value of the diagonal entry of row i
    virtual void GetDiagEntry( UInt i, T &v ) const {
      v = data_[rowPtr_[i]];
    }

    T GetDiagEntry(unsigned int row) const { return data_[rowPtr_[row]]; }

    //! This routine adds the value of v to the matrix entry at (i,j)

    //! This routine can be used to add the value of the parameter v to the
    //! matrix entry belonging to the index pair (i,j). If the index pair
    //! belongs to a matrix entry in the strictly lower triangular part of the
    //! matrix, then the method will ignore the call and return without
    //! altering the matrix object. Thus, the value v will only be added,
    //! if it belongs to the upper triangular part of the matrix. This method
    //! relies on the possibly slow "dense matrix style" access.
    void AddToMatrixEntry( UInt i, UInt j, const T& v );

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
    void Scale( Complex factor );
    
    //! \copydoc StdMatrix::Scale(Double, std::set<UInt>, std::set<UInt>)
    void Scale( Double factor, 
                const std::set<UInt>& rowIndices,
                const std::set<UInt>& colIndices );

    //! Add the multiple of a matrix to this matrix.

    //! The method adds the multiple of a matrix to this matrix,
    //! \f$A = A + \alpha B\f$. In doing so the sparsity structure of the
    //! matrix mat is assumed to be identical to this matrix' structure.
    void Add( const Double alpha, const StdMatrix& mat );

    //! \copydoc StdMatrix::Add(Double,StdMatrix,std::set<UInt>,std::set<UInt>)
    void Add( const Double a, const StdMatrix& mat,
              const std::set<UInt>& rowIndices,
              const std::set<UInt>& colIndices );
    //@}


    // ========================================================================
    // BLOCK STRUCTURE RELATED METHODS
    // ========================================================================
    //@{
    //! \name Block Structure Related Methods

    //! \copydoc StdMatrix::GetNumBlocks
    virtual void GetNumBlocks( UInt& nRowBlocks, UInt& nColBlocks, 
                               UInt& numBlocks ) const;

    //! \copydoc StdMatrix::GetDiagBlock
    virtual void GetDiagBlock( UInt blockRow, DenseMatrix& diagBlock ) const;
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


    std::string Dump() const;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. Since it is a CRS_Matrix the output
    //! format will be symmetric, co-ordinate, real or complex.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    //! \param comment string to be inserted into file header
    void ExportMatrixMarket(const std::string& fname, const std::string& comment = "" ) const;

    //! Copy the matrix into full-storage format 

    //! Trasforms matrix into vector containing all upper
    // triangle elements further usage in openCFS
    void CopySCRSMatrix2Vec(Complex* &dataVec);

    //@}


    // =======================================================================
    // DIRECT ACCESS TO SCRS DATA STRUCTURE
    // =======================================================================

    //! \name Direct access to SCRS data structure
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
      return rowPtr_;
    }

    //! Get diag pointer array
    UInt* GetDiagPointer() {
      return diagPtr_;
    }

    //! Get diag pointer array
    const UInt* GetDiagPointer() const{
      return diagPtr_;
    }

    //! Get row pointer array (read only)
    const UInt* GetRowPointer() const {
      return rowPtr_;
    }

    //! Get array of column indices
    UInt* GetColPointer() {
      return colInd_;
    }

    //! Get array of column indices (read only)
    const UInt* GetColPointer() const {
      return colInd_;
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

    //! Number of stored non-zero entries

    //! This variable keeps track of the number of stored matrix entries.
    //! This value equals the number of non-zero entries in the lower
    //! triangle of the matrix including the main diagonal. In contrast nnz_
    //! stores the number of non-zero matrix entries in the complete matrix.
    UInt numEntries_;

    //! Array containing starting indices of the different matrix rows
    UInt *rowPtr_;

    //! Array containing the column indices of the non-zero matrix entries
    UInt *colInd_;

    //! Array containing the indices of the diagonal matrix entries

    //! This array stores for each matrix row k the index into the #data_
    //! array at which the corresponding diagonal matrix entry is stored.
    //! In the case that the diagonal entry is not part of the sparsity
    //! pattern, or there is no position (k,k) in the matrix, e.g. in case
    //! of a rectangular matrix, an index of 0 is stored.
    UInt* diagPtr_;


    //! Array containing the (potentially) non-zero matrix entries
    T* data_;

    //! Final representation of blocks
    StdVector<std::pair<UInt,UInt> > rowBlocks_;
    StdVector<std::pair<UInt,UInt> > colBlocks_;
    
    //! Flag indicating whether the sparsity pattern belongs to a pattern pool

    //! This is a pointer to a PatternPool object. If the pointer is NULL
    //! in an instance of this matrix class it indicates that the sparsity
    //! pattern belongs to that instance, i.e. the destructor must de-allocate
    //! the associated memory. If the pointer is not NULL this signals that
    //! the sparsity pattern is shared between different matrix instances and
    //! memory is administrated by the respective PatternPool object. Thus,
    //! the destructor must de-register itself with that object.
    PatternPool *patternPool_;

    //! Attribute storing identifier for sparsity pattern if shared

    //! This attribute stores the pattern identifier for the sparsity pattern
    //! of the matrix, if the latter belongs to the PatternPool object to
    //! which patternPool_ points to. If the matrix owns the pattern then the
    //! value is NO_PATTERN_ID.
    PatternIdType patternID_;

  };


} // namespace


#endif // OLAS_SCRS_MATRIX_HH
