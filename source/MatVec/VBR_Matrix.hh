#ifndef OLAS_VBR_MATRIX_HH
#define OLAS_VBR_MATRIX_HH

#include <iostream>



#include "SparseOLASMatrix.hh"
#include "Vector.hh"
#include "CoordFormat.hh"

namespace CoupledField {

  
  //! Class for representing matrices in Variable Block Row (VBR) format

  //! This class allows to store an unsymmetric sparse matrix in the
  //! <em>Variable Block Row</em> (VBR) format.
  //! This format is similar to CRS, but instead of single entries,
  //! it stores dense matrices with variable size. Within each block,
  //! all the zeros are stores explicitly. 
  //! 
  //! \note
  //!   - Insertion takes two searches: one for row and one for column. We
  //!     use already an averageBlockSize to estimate the blockRow of a row.
  //!     Additional ideas could be to re-use last block (assuming subsequent
  //!     entries are in same blocks), use several insertions at one time
  //!     or we could use an average distance of the blockCols to the 
  //!     diagonal to get rough estimate for the starting value.
 
  template<typename T>
  class VBR_Matrix : public SparseOLASMatrix<T> {

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
    //! VBR_Matrix objects, to initialization of matrix values and to the
    //! setting and manipulation of the sparsity pattern of the matrix.
    //@{

    //! Default Constructor
    VBR_Matrix();

    //! Copy Constructor

    //! This is a deep copy constructor. It is mainly provided to allow for
    //! certain sparse arithmetic operations of the SBM_Matrix class.
    //! \note  The constructed matrix is a copy and thus has the same layout
    //!        as the original matrix it is copied from.
    VBR_Matrix( const VBR_Matrix<T> &origMat );

    
    //! Constructor with memory allocation

    //! This constructor generates a matrix of given dimensions and allocates
    //! space to store the specified number of non-zero entries.
    //! Initially all matrix entries are set to zero.
    //! \note  Since no column index information is provided the constructor
    //!        sets the currentLayout_ format to UNSORTED.
    //! \param nrows number of rows of the matrix
    //! \param ncols number of columns of the matrix
    //! \param nnz   number of non-zero matrix entries
    VBR_Matrix( UInt nrows, UInt ncols, UInt nnz );

    //! Destructor

    //! The default destructor is deep. It frees all memory dynamically
    //! allocated for attributes of this class.
    virtual ~VBR_Matrix();

    //! Setup the sparsity pattern of the matrix
  
    //! This method uses the graph to generate the internal information needed
    //! for representing the sparsity pattern of the matrix correctly. Before
    //! this function has not been called the matrix is not ready to be used.
    //! \note
    //! The implementation of this routine expects the GetGraphMethod
    //! of the BaseGraph to return the list of nodes (including the
    //! node itself, in case of a self-reference) in lexicographic
    //! ordering, i.e. in increasing number of node indices. Thus, the
    //! resulting %VBR_Matrix is in LEX sub-format.
    //! \note
    //! If the graph provides information about minor blocks (via 
    //! BaseGraph::GetRowBlockDefinition or BaseGraph::GetColBlockDefinition),
    //! the internal block structure will be set in accordance to that layout.
    //! If not block information is present (either for rows or columns),
    //! a default block size of 1 is assumed. 
    void SetSparsityPattern( BaseGraph &graph );

    //! Copy the sparsity pattern of another matrix to this matrix 
    
    //! This method copies the sparsity pattern of the source matrix to the
    //! current matrix
    void SetSparsityPattern ( const StdMatrix &srcMat );
    
    //! Re-set internal dimensions

    //! This method can be used to re-size the matrix. If the internal
    //! dimensions differ from the ones previously set, then the internal
    //! data arrays will be thrown away and re-allocated with the new sizes.
    //! In this case all matrix entries and sparsity pattern information is
    //! lost, of course an the currentLayout_ attribute is set to UNSORTED.
    void SetSize( UInt nrows, UInt ncols, UInt nnz );
 
    //! Initialize matrix to zero
    
    //! This matrix sets all matrix entries to zero. Note that all entries
    //! in the non-zero structure of the matrix get a zero value. However,
    //! the structure/data layout remains the same. All positions in the
    //! structure can still be over-written with non-zero values.
    void Init();

    /** Transpose the matrix, which is transforming from the row to a
     * column storage format. So this method converts from the internal row
     * storage to a column storage. The space has to be provided!
     * @param col_ptr required space is number of rows/columns +1 for one based +1 for tail
     * @param row_ptr required space is nnz +1 for one based
     * @param data_ptr required space is nnz +1 for one based */
   void Transpose(UInt* col_ptr, UInt* row_ptr, T* data_ptr) const;

    //@}


    // =======================================================================
    // QUERY / MANIPULATE MATRIX FORMAT
    // =======================================================================

    //! \name Query matrix / manipulate format
    //@{

    //! Return the storage type of the matrix
  
    //! The method returns the storage type of the matrix. The latter is
    //! encoded as a value of the enumeration data type MatrixStorageType.
    //! \return VAR_BLOCK_ROW
    inline BaseMatrix::StorageType GetStorageType() const {
      return BaseMatrix::VAR_BLOCK_ROW;
    }
    
    //! \copydoc StdMatrix::GetNumBlocks
    void GetNumBlocks(UInt& nRowBlocks, UInt& nColBlocks, UInt& numBlocks ) const ;
    
    
    //! Get effective number of nonzeros (including padding zeros)
    UInt GetEffNnz() const {
      return nNzEff_;
    }
    
    //! Get average row block size
    Double GetAvgRowBlockSize() const {
      return 1.0 / oneOverBlockSize_;
    }
    
    //! Get number of scalar non-diagonal-block entries
    UInt GetNumOffDiagonalEntries() const {
      return numNonDiagBlockEntires_;
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
                         bool setCounterPart );

    //! Get a specific matrix entry

    //! This function returns a reference to the Matrix entry in 
    //!  row i and column j. Depending on the entry type of the matrix,
    //!  v might be a Double, Complex or a tiny Matrix of either type
    void GetMatrixEntry( UInt i, UInt j, T &v ) const;

    //! Set the diagonal entry of row i to the value of v
    void SetDiagEntry( UInt i, const T &v ) {
      UInt aux = diagPtr_[i];
      data_[aux] = v;
    }

    //! Returns the value of the diagonal entry of row i
    void GetDiagEntry( UInt i, T &v ) const {
      UInt aux = diagPtr_[i];
      v = data_[aux];
    }

    T GetDiagEntry(unsigned int row) const { return data_[diagPtr_[row]]; }

    //! This routine adds the value of v to the matrix entry at (i,j)

    //! This routine can be used to add the value of the parameter v to the
    //! matrix entry belonging to the index pair (i,j). This method relies on
    //! the possibly slow "dense matrix style" access.  
    void AddToMatrixEntry( UInt i, UInt j, const T& v );

    //! Return the length (i.e. number of non-zero entries) of i-th row
    UInt GetRowSize(UInt i) const {
      return rowPtr_[i+1] - rowPtr_[i];
    }

    //! Set the length (i.e. number of non-zero entries) of i-th row

    //! The method can be used to set the length (i.e. number of non-zero
    //! entries) of row i to the value size.
    void SetRowSize( UInt i, UInt size ) {
      rowPtr_[i+1] = rowPtr_[i] + size;
    }

    //! Insert a value/column index at position (row,pos)
    void Insert( UInt row, UInt col, UInt pos );
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

    //! \copydoc StdMatrix::Scale(Double, std::set<UInt>, std::set<UInt>)
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
    void Add( const Complex a, const StdMatrix& mat, const std::set<UInt>& rowIndices, const std::set<UInt>& colIndices );
    //@}


    // =======================================================================
    // I/O operations
    // =======================================================================

    //! \name I/O operations
    //@{

    //! Print the matrix to an output stream
    
    //! Prints the matrix in indexed format to the output stream os. This
    //! method is primarily intended for debugging.
    std::string ToString( char colSeparator = ' ',
                          char rowSeparator = '\n' ) const;
    
    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. Since it is a VBR_Matrix the output
    //! format will be symmetric, co-ordinate, real or complex.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    //! \param comment string to be inserted into file header
    void ExportMatrixMarket(const std::string& fname, const std::string& comment = "") const;
    //@}


    // =======================================================================
    // DIRECT ACCESS TO VBR DATA STRUCTURE
    // =======================================================================

    //! \name Direct access to VBR data structure
    //! The following methods allow a direct manipulation of the internal
    //! data structures. Though this does violate the concept of data
    //! encapsulation, these methods are provided for the sake of efficiency.
    //@{
    
    //! \copydoc StdMatrix::GetDiagBlock
    void GetDiagBlock( UInt blockRow, DenseMatrix& diagBlock ) const; 
    
    //! \copydoc StdMatrix::SetDiagBlock
    void SetDiagBlock( UInt blockRow, const DenseMatrix& diagBlock );
    //@}

  private:
    
    //! Auxiliary method for setting up the diagonal index array

    //! This auxiliary method uses the existing VBR structure information,
    //! i.e. the rowPtr_ and colInd_ arrays, to determine the index positions
    //! of the diagonal matrix entries in the data_ and colInd_ arrays. 
    void FindDiagonalEntries();
    
    //! Find for given row and col index the rowBlock and colBlock
    
    //! This method find for a given row and column the corresponding
    //! rowBlock and columnBlock. Thus, it has to perform two search runs,
    //! one for the row and one for the column.
    //! The former one is accelerated by assuming an average blockSize
    //! and setting the startBlock to row * oneOverBlockSize.
    //! \param row row index of the searched position
    //! \param col column index of the searched position
    //! \param bRow row block index of the searched position
    //! \param bCol column block index of the searched position
    //! \param offset offset to the #data_ array, i.e. the entry
    //!        can be accessed by _data[offset] 
    void FindBlock(UInt row, UInt col, UInt& bRow, UInt& bCol,
                   UInt& offset ) const;
    
    //! Effective number of nonzero elements (including padding zeros)
    
    //! This value contains the total number of entries of the matrix,
    //! including additional zeros introduced due to the block structure.
    //! This value is the effective length of the #data_ array
    UInt nNzEff_;
    
    //! Number of block rows
    UInt nbRows_;
    
    //! Number of block columns
    UInt nbCols_;
    
    //! Total number of non-zero blocks
    
    //! This is total number of non zero blocks.
    //! \note: #nBlocks_ != #nbRows_ x #nbCols_
    UInt nBlocks_;
    
    //! Array containing the starting row positions of each block row
    
    //! This array stores for each blockRow the starting rowNumber.
    //! The i-th block row starts at row #bRow_[i] and ends at row 
    //! #bRow_[i+1]-1. The last entry stores the number of row #nrows_.
    //! The array has a length of #nbRows_ + 1. 
    UInt *bRow_;
    
    //! Array containing the starting column positions of each column block
    
    //! This array stores for each blockCol the starting column number.
    //! The j-th block column starts at at column #bCol_[j] and end at column
    //! #bCol_[j+1]-1. The last entry stores the number of columns #ncols_;
    //! The array has a length of #nbCols_ +1.
    UInt *bCol_;
    
    //! Array containing the starting index for each block
 
    //! This array stores for each block the starting index within #data_.
    //! The b-th block starts at position #valPtr_[b] in the array #data_.
    //! The last element #valPtr_[nBlocks_] contains the total number of 
    //! nonzero scalar values #nnz_. 
    //! The array has a length of #nBlocks_+1.
    UInt *valPtr_;
    
    //! Array containing the column indices of the blocks
    
    //! This array stores fore each block the block column index.
    //! The b-th block begins at column #bCol_[colInd_[b]].
    //! The array has a length of #nBlocks_. 
    UInt *colInd_;
    
    //! Array containing the starting block row indices in #colInd_
    
    //! This array stores for each block row the starting offset within
    //! #colInd_. The i-th block row starts at position rowPtr_[i] in
    //! #colInd_. The last entry contains the total number of blocks
    //! #rowPtr_[nbRows_] = #nBlocks_.
    //! The array has a length of #nbRows_ + 1;
    UInt *rowPtr_;
    
    //! Array containing the indices of the scalar diagonals

    //! This array points directly to the scalar diagonal entries of the
    //! matrix (in case it is square).
    UInt *diagPtr_;
    
    //! Array containing the indices to the diagonal blocks
    
    //! This array contains the indices within the #colInd_ array of the
    //! diagonal blocks (in case the matrix is square).
    UInt *diagBlockPtr_; 

    //! Array containing block-wise all entries. Blocks are in row format.
    
    //! This array contains all values stores in a block-wise fashion.
    //! Each block is dense and can contain zeros. The layout of the blocks
    //! is row-wise.
    T* data_;
    
    //! Helper variable: 1 / average block width
    float oneOverBlockSize_;
    
    //! Number of scalar entries not contained in the diagonal block
    
    //! This value can be used to estimate, how many scalar entries are
    //! off-diagonal, i.e. are not contained in the diagonal blocks.
    UInt numNonDiagBlockEntires_;
    
  };

}
#endif
