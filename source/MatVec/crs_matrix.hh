// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_CRS_MATRIX_HH
#define OLAS_CRS_MATRIX_HH

#include <stddef.h>
#include <complex>
#include <iostream>
#include <string>

#include "General/defs.hh"
#include "General/exception.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/sparseolasmatrix.hh"
#include "def_expl_templ_inst.hh"

namespace CoupledField {

class BaseGraph;
class StdMatrix;
template <typename T> class CoordFormat;
template <typename T> class Vector;
#ifdef USE_MULTIGRID
  //forward declaration
  template <typename T>
  class PreMatrix;
#endif
  
  //! Class for representing matrices in Compressed Row Format

  //! This class allows to store an unsymmetric sparse matrix in the
  //! <em>Compressed Row Storage Format</em> (CRS).
  //! In this format all non-zero entries of the matrix are stored together
  //! with the corresponding index information. This happens in a row-wise
  //! fashion in the following way. In a data array (#data_) the non-zero
  //! entries are stored for one row after the other. The column indices
  //! are stored in an integer array (#colInd_) in the same fashion.
  //! A second integer array (#rowPtr_) stores the starting indices of the
  //! different rows in the two other vectors.\n\n
  //! The CRS format does not specify any order in which the non-zero entries
  //! of a row should be stored. For some operations, like e.g. computing
  //! a matrix-vector-product, such an ordering would be irrelevant. Other
  //! operations, like e.g. a matrix-factorisation on the other hand, would
  //! prefer a lexicographic ordering, while a typical smoother like
  //! Gauss-Seidel requires a fast way to extract the diagonal entry of a
  //! matrix row.
  //! This class, thus, supports the following sub-formats of CRS with
  //! respect to the ordering of row entries:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="2" align="center">
  //!         <b>Supported sub-formats</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>UNSORTED</td>
  //!       <td>row entries are stored in arbitrary order</td>
  //!     </tr>
  //!     <tr>
  //!       <td>LEX</td>
  //!       <td>row entries are stored in lexicographic ordering with
  //!         respect to increasing column index</td>
  //!     </tr>
  //!     <tr>
  //!       <td>LEX_DIAG_FIRST</td>
  //!       <td>same as LEX, but for each row the diagonal entry
  //!           (if it exists!) is stored in the leading position
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \n
  //! The sub-format of the matrix may be altered by the ChangeLayout()
  //! method and the current sub-format may be queried by the
  //! GetCurrentLayout() method.
  //! Fast access to diagonal entries is provided by storing in a
  //! fourth array the array positions of the diagonal entries. In the case
  //! that row k does not have a (non-zero) diagonal entry a zero is stored
  //! as array position. The latter is sensible, since this class, like all
  //! %OLAS matrix classes, uses one-based indexing.
  //! \note
  //! - Debugging for this class can be turned on by defining the macro
  //!   DEBUG_NEWCRS_MATRIX
  //! - Is is possible to generate a %CRS_Matrix object that is stored in a
  //!   certain sub-format, without the object itself being aware of this
  //!   fact. One may, e.g. use a constructor involving a PreMatrix that
  //!   is in LEX format, but the CRS_Matrix will set the currentLayout_
  //!   attribute to UNSORTED in this case. In this case it is possible to
  //!   inform the object on its actual sub-format by using the
  //!   SetCurrentLayout() method. However, care must be taken in doing so,
  //!   in order to avoid serious inconsistencies.
  template<typename T>
  class CRS_Matrix : public SparseOLASMatrix<T> {

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

    //! Enumeration data type describing storage layout of the matrix

    //! This enumeration data type describes the storage layout of the matrix
    //! entries in each row. Currently there are the following values
    //! - UNSORTED
    //! - LEX
    //! - LEX_DIAG_FIRST
    //!
    //! See the class' description for details on the values.
    typedef enum { UNSORTED, LEX, LEX_DIAG_FIRST } subFormat;

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION AND SETUP
    // =======================================================================

    //! \name Construction, Destruction and Setup
    //! The following methods are related to construction and destruction of
    //! CRS_Matrix objects, to initialisation of matrix values and to the
    //! setting and manipulation of the sparsity pattern of the matrix.
    //@{

    //! Default Constructor
    CRS_Matrix() {
      colInd_           = NULL;
      rowPtr_           = NULL;
      diagPtr_          = NULL;
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
    CRS_Matrix( const CRS_Matrix<T> &origMat );

#ifdef USE_MULTIGRID
   
    //! Constructor using a PreMatrix
    
    //! If a matrix has been created not from a graph, but by successively
    //! inserting non-zero entries into a PreMatrix object, this method can be
    //! used to convert that PreMatrix object to a matrix in CRS format. The
    //! PreMatrix object itself remains unchanged by this operation.
    //! \note Querying the PreMatrix object can for its layout status can,
    //!       at the moment, only tell us if it is LEX_DIAG_FIRST, or not.
    //!       Thus, using this constructor results in a matrix with either
    //!       LEX_DIAG_FIRST or UNSORTED as layout format. If you know that
    //!       the PreMatrix has LEX as sub-format use SetCurrentLayout() to
    //!       fix this.
    //! \param nrows number of rows
    //! \param ncols number of columns
    //! \param premat the pre-matrix with the matrix data
    CRS_Matrix( UInt nrows, UInt ncols, const PreMatrix<T>& premat );
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
    CRS_Matrix( UInt nrows, UInt ncols, UInt nnz )
      : diagPtr_( NULL ) {


      // assign basic properties
      this->nnz_   = nnz;
      this->nrows_ = nrows;
      this->ncols_ = ncols;

      // Allocate memory for row and diagonal pointers
      NEWARRAY( rowPtr_ , UInt, this->nrows_ + 1 );
      NEWARRAY( diagPtr_, UInt, this->nrows_     );
      rowPtr_[0]  = 0;
      diagPtr_[0] = 0;

      // Allocate memory for column indices
      NEWARRAY( colInd_, UInt, this->nnz_ );

      // Allocate memory for matrix entries and initialise to zero
      NEWARRAY( data_, T, this->nnz_ );
      Init();

      // No data provided, thus we assume unsortedness
      currentLayout_ = UNSORTED;
    }

    //! Constructor employing a CoordFormat matrix

    //! This constructor variant is a sort of copy constructor. It expects as
    //! input a sparse matrix in coordinate format. It uses this information
    //! to generate and populate an equivalent CRS_Matrix.
    //! \note
    //! - The construct performs e re-ordering to LEX at the end
    //! \param sparseMat sparse matrix in coordinate format
    CRS_Matrix( CoordFormat<T> &sparseMat );

    //! Destructor

    //! The default destructor is deep. It frees all memory dynamically
    //! allocated for attributes of this class.
    virtual ~CRS_Matrix() {
      delete []( data_ );
      delete []( rowPtr_ );
      delete []( colInd_ );
      delete []( diagPtr_ );
    }

    //! Setup the sparsity pattern of the matrix
  
    //! This method uses the graph to generate the internal information needed
    //! for representing the sparsity pattern of the matrix correctly. Before
    //! this function has not been called the matrix is not ready to be used.
    //! \note
    //! The implementation of this routine expects the GetGraphMethod
    //! of the BaseGraph to return the list of nodes (including the
    //! node itself, in case of a self-reference) in lexicographic
    //! ordering, i.e. in increasing number of node indices. Thus, the
    //! resulting %CRS_Matrix is in LEX sub-format.
    void SetSparsityPattern( BaseGraph &graph );

    //! Copy the sparsity pattern of this matrix to another 
    
    //! This method copies the sparsity pattern of this matrix and
    //! sets it in another matrix
    void CopySparsityPattern ( StdMatrix &mat ) const;
    
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
    inline void Init() {
      for ( UInt i = 0; i < this->nnz_; i++ ) {
        data_[i] = 0;
      }
    }

    /** Transpose the matrix, which is transforming from the row to a
     * column storrage format. So this method converts from the internal row
     * storage to a column storrage. The space has to be provided!
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
    //! \return SPARSE_NONSYM
    inline BaseMatrix::StorageType GetStorageType() const {
      return BaseMatrix::SPARSE_NONSYM;
    }

    //! Query the current sub-format of the matrix object
    inline subFormat GetCurrentLayout() const {
      return currentLayout_;
    }

    //! Set the current sub-format of the matrix object
    inline void SetCurrentLayout( const subFormat myLayout ) {
      currentLayout_ = myLayout;
    }

    //! Sort all rows according to the specified sub-format

    //! Calling this method will sort all matrix rows in such a manner that
    //! they conform to the specified sub-format layout. This method allows
    //! e.g. to assemble a CRS_Matrix in an unsorted fashion, using a
    //! PreMatrix, and to sort it afterwards. After the transformation was
    //! performed the method adapts the currentLayout_ attribute accordingly.
    //! \note
    //! - The function does not touch the matrix, if the currentLayout_
    //!   attribute matches the sub-format specified for sorting.
    //! - ChangeLayout will re-built the information in the diagPtr_ array
    //! \param newLayout specifies the new sub-format into which the layout
    //!                  of the matrix should be transformed
    void ChangeLayout( const subFormat newLayout );

    //! old way of sorting, will be removed in future release

    //! \deprecated This method is a remainder from the first implementation
    //! of the %CRS_Matrix class. It assumes that LEX_DIAG_FIRST is the
    //! standard format and will sort the matrix accordingly. The method will
    //! be removed in a future release
    void SortConformingToLayout( const bool forced = true ) {
      ChangeLayout( LEX_DIAG_FIRST );
    }
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

    //! Return the diagonal entry of row i
    inline
    T& GetDiag( UInt i ) {
      UInt aux = diagPtr_[i];
      return data_[aux];
    }

    //! Return the diagonal entry of row i (read only)
    inline
    const T& GetDiag( UInt i ) const {
      UInt aux = diagPtr_[i];
      return data_[aux];
    }

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

    //! Determine maximum absolute value of diagonal entries

    //! This method determines the the maximal absolute value (on the scalar)
    //! level of the entries on the main diagonal of the matrix.
    Double GetMaxDiag() const;

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
    std::string ToString( char colSeparator = ' ',
                          char rowSeparator = '\n' ) const;
    
    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. Since it is a CRS_Matrix the output
    //! format will be symmetric, co-ordinate, real or complex.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    //! \param comment string to be inserted into file header
    void Export( const char *fname, const char *comment = NULL ) const;

    //! Transform layout enumeration value to string format
    static std::string Enum2String( CRS_Matrix<T>::subFormat myEnum ) {
      std::string out;
      if ( myEnum == LEX ) {
        out = "LEX";
      }
      else if ( myEnum == LEX_DIAG_FIRST ) {
        out = "LEX_DIAG_FIRST";
      }
      else if ( myEnum == UNSORTED ) {
        out = "UNSORTED";
      }
      else {
        EXCEPTION("Congratulations! You have found a missing case "
                 << "implementation! Parameter myEnum = " << myEnum);
      }
      return out;
    }
    //@}


    // =======================================================================
    // DIRECT ACCESS TO CRS DATA STRUCTURE
    // =======================================================================

    //! \name Direct access to CRS data structure
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

    //! Get row pointer array (read only)
    const UInt* GetRowPointer() const {
      return rowPtr_;
    }

    //! Get diagonal pointer array
    UInt* GetDiagPointer() {
      return diagPtr_;
    }

    //! Get diagonal pointer array (read only)
    const UInt* GetDiagPointer() const {
      return diagPtr_;
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

    //! Static auxiliary function for ChangeLayout()

    //! This is an auxiliary function for the ChangeLayout() method.
    //! It implements a quick sort algorihm, where the comparison is based on
    //! the column indices passed in cols, and the values in data & cols
    //! are reordered accordingly. The resulting arrays are in ascending
    //! column index ordering.
    //! \param cols    integer array containing the column indices of a row
    //! \param data    array containing the entry data of a row
    //! \param length  number of valid entries starting at cols and data,
    //!                respectively
    static void QuickSort( UInt *const cols, T *const data,
                           const UInt length );

    //! Auxilliary method for setting up the diagonal index array

    //! This auxilliary method uses the existing CRS structure information,
    //! i.e. the rowPtr_ and colInd_ arrays, to determine the index positions
    //! of the diagonal matrix entries in the data_ and colInd_ arrays. It
    //! is used e.g. by SortConformingToLayout() in some cases and the
    //! constructor that employs a PreMatrix.
    void FindDiagonalEntries();

    //! Specialised conversion routine LEX -> LEX_DIAG_FIRST
    void SortLex2LexDiagFirst();

    //! Specialised conversion routine LEX_DIAG_FIRST -> LEX
    void SortLexDiagFirst2Lex();

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
    UInt *diagPtr_;

    //! Array containing the (potentially) non-zero matrix entries
    T* data_;

    //! Attribute storing the current sub-format of the matrix layout
    subFormat currentLayout_;

  };

}

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "crs_matrix.cc"
#endif

#endif
