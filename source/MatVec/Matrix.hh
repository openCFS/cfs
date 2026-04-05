#ifndef FILE_MATRIX_2004
#define FILE_MATRIX_2004

#include <def_build_type_options.hh>
#include <boost/container/small_vector.hpp>

#include "MatVec/promote.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/opdefs.hh"
#include "Utils/ToolsCore.hh"

#ifdef USE_EXPRESSION_TEMPLATES
#include "MatVec/exprt/xpr2.hh"
#endif

#include "DenseMatrix.hh"
#include "General/Exception.hh"

namespace CoupledField
{      

  //! Forward class declaration
  template<class TYPE> class Vector;
  template<class TYPE> class StdVector;

  
  /** Matrix is mathematical matrix for integers, doubles and complex.
   * Similar to Vector it uses a boost::small_vector with stack data for 6x6 entries.
   * If resized above, transparently dynamic memory is allocated on the heap.
   * check with cfs -d to see the actual (remaining) allocations in the .info.xml */
  template<class TYPE>
#ifdef USE_EXPRESSION_TEMPLATES
  class Matrix: public DenseMatrix, public Dim2<TYPE, Matrix<TYPE> >
#else
  class Matrix: public DenseMatrix 
#endif
  {
  public:

    friend class Matrix<Double>;
    friend class Matrix<Complex>;

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION, INITIALIZATION, RESIZING
    // =======================================================================
    
    //! \name Construction, Destruction, Initialization and Resizing

    /** @see buffer_ */
    static constexpr unsigned int CalcBufferSizeBytes(unsigned int nRows, unsigned int nCols) {
       return sizeof(TYPE*) * nRows + sizeof(TYPE) * nCols * nRows;
    }

    //@{ 
    //! Default constructor 
    
    //! Creates an empty matrix of size 0x0
    Matrix( );
    
    //! Constructor for matrix with given size
    
    //! Creates a matrix of size nRows x nCols, initialized
    //! with zeroes
    //! \param nRows (input) Number of rows
    //! \param nCols (input) Number of columns)
    explicit Matrix( const UInt nRows, const UInt nCols );

    //! Constructor from array of column-vectors

    //! Creates a matrix from an array of column vectors.
    //! The number of rows will be the number of entries of one vector
    //! and the number of columns will be \a numVec.
    //! \param numVec number of column vector contained in vecs
    Matrix( const UInt numVec, const Vector<TYPE> * const vecs );

    //! Default copy constructor
    Matrix( const Matrix & );

    //! Templatized copy constructor

    //! Generalized copy constructor. It is only implemented 
    //! to create a complex-valued matrix from a given real-valued one.
    //    template<class T2> 
    //    Matrix( const Matrix<T2> &  );

      //! Destructor
    virtual ~Matrix( );
      
    //! Initialize matrix with a given scalar entry.

    //! Initializes the matrix with a given scalar entry
    //! If no entry given, it gets initialized with zeroes.
    //! \param val (input,opt.) Entry the matrix gets initialized with
    //! \note This method does not change the size of the matrix
    inline void Init()
    {
      InitValue();
    }

    inline void InitValue( const TYPE val = TYPE() )
    {
      if(size_row_*size_col_ > 0)
        std::fill(&data_[0][0],&data_[0][size_row_*size_col_], val);
    }
    
    //! Change the size of the matrix

    //! Change size of general matrix 
    //! \param nRows (input) Number of rows
    //! \param nCols (input) Number of columns
    void Resize(const UInt nRows, const UInt nCols );

    //! Changes the size so that the matrix gets quadratic
    
    //! Changes the size of the matrix according to \a size.
    //! \param size (input) Number of rows / columns
    void Resize(const UInt size);

    /** Resize if necessary to the other matrix */
    void Resize(const Matrix<TYPE>& other);
    
    //@}
    
    // =======================================================================
    // GENERAL INFORMATION
    // =======================================================================
    
    //! \name General Matrix Information
    
    //@{

    //! Get entry type of matrix
    inline BaseMatrix::EntryType GetEntryType()
    {
      return  EntryType<TYPE>::M_EntryType;
    }

    /** Is the Matrix quadratic */
    bool IsQuadratic() const { return GetNumRows() == GetNumCols(); }

    //! Check if the matrix is symmetric

    /** return true, if the matrix is symmetric.
     * Does a binary comparison
     * @see IsSymmetric(bool) */
    bool IsSymmetric() const;

    /** return true, if the matrix is complex (imag = 0 for all entries). */
    bool IsComplex() const;

    /** symmetry check with eps sensitivity. */
    inline bool IsSymmetric(double eps) const
    {
      if(!IsQuadratic())
        return false;

      for(UInt i = 1; i < size_row_; ++i)
        for(UInt j = i+1; j < size_col_; ++j)
          if(!close(data_[i][j], data_[j][i]))
            return false;

      return true;
    }

    /** check if the matrix is a Hermitian matrix. In the non complex case symmetry is checked.
     * @param eps if true use close() to compare the values by an eps  */
    bool IsHermitian(double eps = 1e-6) const;


    //! Get the number of rows
    inline UInt GetNumRows() const
    {       
      return size_row_;
    }
    
    //! Get the number of columns
    inline UInt GetNumCols() const
    {       
      return size_col_;
    }

    /** rows time columns */
    unsigned int GetNumEntries() const { return size_row_ * size_col_; }

    //@}

    // =======================================================================
    // OBTAIN / MANIPULATE MATRIX ENTRIES
    // =======================================================================
    
    //! \name Obtain / Manipulate Matrix Entries

    //@{
    //! General access operator

    //! Access operator of one entries
    //! \param row (input) Row number
    //! \param col (input) Column number
    inline TYPE & operator()( UInt row, UInt col)
    {
      assert(row < size_row_ && col < size_col_);
      return data_[row][col];
    }

    //! General access operator (const)

    //! Access operator of one entries
    //! \param row (input) Row number
    //! \param col (input) Column number
    inline TYPE operator()( UInt row, UInt col ) const
    {
      assert(row < size_row_ && col < size_col_);
      return data_[row][col];
    }
    
    //! Returns pointer to row \a row
    inline TYPE * operator[]( const UInt i ) const
    { 
#ifdef CHECK_INITIALIZED
      if (size_row_ == 0 || size_col_ == 0) 
        EXCEPTION( "undefined Matrix" );
#endif

#ifdef CHECK_INDEX
      if (i >= size_row_) 
        EXCEPTION( "invalid index " << i << " for " << size_row_ << " rows");
#endif

      return data_[i];
    }

    //! Get the entry 'val' at position (row,col) in the matrix
    
    //! Return entry at position (\a row, \a col) in the matrix
    //! \param row (input) row index of entry
    //! \param col (input) column index of entry
    //! \param val (output) on return contains value of entry
    inline void GetEntry( const UInt row, const UInt col, TYPE & val ) const
    {
      assert(row < size_row_ && col < size_col_);
      val = data_[row][col];
    }
     
    //! Set the entry 'val' at position (\a row, \a col) in the matrix
    //! \param row (input) Row of entry
    //! \param col (input) Column of entry
    //! \param val (input) Value to be set
    inline void SetEntry(const UInt row, const UInt col, const TYPE val)
    {
      assert(row < size_row_ && col < size_col_);
      data_[row][col] = val;
    }

    //! Add'val' to the matrix entry at position (\a row, \a col) in the 
    //! matrix
    //! \param row (input) Row of entry
    //! \param col (input) Column of entry
    //! \param val (input) Value to be added
    inline void AddToEntry( const UInt row, const UInt col, const TYPE val )
    {
      assert(row < size_row_ && col < size_col_);
      data_[row][col] += val;
    }
    
    /** give a specific row */
    void GetRow(Vector<TYPE>& vec_out, UInt row) const;

    /** give a specific column */
    void GetCol(Vector<TYPE>& vec_out, UInt col) const;

    //! Set the given column 'vec_in' at position (:, \a col) in the matrix
    //! \param vec_in (input) data column
    //! \param col (input) Column idx to be set
    void SetCol(Vector<TYPE>& vec_in, UInt col) const;

    /** For each row the minimum over all columns. For Complex see Vector::Min() */
    void GetColMin(Vector<TYPE>& vec_out) const;

    /** See GetColMin() */
    void GetColMax(Vector<TYPE>& vec_out) const;

    /* get maximal and minimal matrix entry */
    TYPE GetMax() const;
    TYPE GetMin() const;

    void GetAbsValues(Matrix<TYPE>& AbsMatrix) const;

    //! Gets the diagonal elements of a  matrix in a one column matrix
    void GetDiagInMatrix( Matrix<TYPE>& columnMat ) const;

    //! Get submatrix indirect by row/col indices

    //! This matrix returns a sub-matrix, defined by the row- and column-
    //! indices.
    //! \param ret requested submatrix of size rowInd.Size() x colInd.Size()
    //! \param rowInd indices of the requested row indices
    //! \param colInd indices of the requested column indices
    void GetSubMatrixByInd( Matrix<TYPE>&ret, 
                            const StdVector<UInt>& rowInd,
                            const StdVector<UInt>& colInd ) const;
        
    
    //! Set submatrix by row/col indices
    
    //! This matrix sets a sub-matrix, defined by the row- and column-
    //! indices.
    //! \param ret  submatrix to set of size rowInd.Size x colInd.Size()
    //! \param rowInd indices of the requested row indices
    //! \param colInd indices of the requested column indices
    void SetSubMatrixByInd( const Matrix<TYPE>&ret, 
                            const StdVector<UInt>& rowInd,
                            const StdVector<UInt>& colInd );
    
    //@}

    // =======================================================================
    // NAMED ARITHMETIC OPERATIONS
    // =======================================================================
    
    //! \name Named Arithmetic Operations
    //@{

    /** Add the multiple of another matrix this = this + fac * mat.
     * If you have mixed types use the tools version of Add */
    void Add(const TYPE fac, const Matrix<TYPE> & mat);
    
    /** Add the multiple of the transpose of another matrix this = this + fac * transpose(mat).
     * If you have mixed types use the tools version of Add */
    void AddT(const TYPE fac, const Matrix<TYPE> & mat);

    /** Set this matrix with a multiple of another matrix.
     * This and a mixed variant is also a stand alone method in tools.
     * Anybody knows how to do the mixed form (complex <- double * complex) here? 
     * this = factor * other_mat
     * @param size_tolerant if set this matrix and other_mat may have different size with 0 entries for the unused. */
    void Assign(const Matrix<TYPE>& other_mat, TYPE factor, bool size_tolerant = false);
    
    /** Set the matrix out of a vector.
     * This matrix is resized. rows times cols needs to match vec.GetSize()
     * @param row_major = true assumes a11, a12, a13, a21, a22, ... (C style)
     *        row_major = false = col_major assumes a11, a21, a31, ... (Fortran style)
     * @see https://en.wikipedia.org/wiki/Row-major_order */
    void Assign(const Vector<TYPE>& vec, unsigned int rows, unsigned int cols, bool row_major);

    //! Perform a matrix-matrix multiplication rMat = this*mMat
    void Mult(const DenseMatrix & mMat, DenseMatrix & rMat) const;

    //! Perform generalized matrix-matrix multiplication using BLAS
    
    //! This method calculates the matrix-matrix product
    //! 
    //!     rMat = alpha * this * mMat + beta * rMat
    //! 
    //! using BLAS optimized d/xgemm method.
    //! If \a trans_a or \a trans_b are set to yes, the corresponding matrix is 
    //! transposed.
    //! \param mMat second argument for the matrix-matrix product
    //! \param rMat final matrix result
    //! \param trans_a true, if own matrix should multiplied in transposed 
    //!        state
    //! \param trans_b true, if own \a mMat should multiplied in transposed 
    //!        state
    //! \param alpha additional scalar factor for matrix-matrix product
    //! \param beta scalar factor for re-use of \a rMat
    //! \param conjugate interpret in complex case the transposed by conjugate complex transpose?
    //! \note Currently we assume the \a rMat to have the correct size already
    //! \note If CFS is compiled without BLAS support, we use as fallback the
    //!       internal matrix-matrix multiplication.
    void Mult_Blas(const Matrix& mMat, Matrix& rMat, 
                   bool trans_a, bool trans_b, TYPE alpha,
                   TYPE beta, bool conjugate = false) const;

    //! Perform a matrix-matrix multiplication rMat = Transpose(this)*mMat
    void MultT(const DenseMatrix & mMat, DenseMatrix & rMat) const;

    //! Perform a matrix-vector multiplication rvec = this*mvec
    void Mult( const SingleVector & mvec, SingleVector & rvec ) const;

    /** scale the matrix by the scalar*/
    void Mult(TYPE scale);

//    //! Perform generalized matrix-vector multiplication using BLAS
//    
//    //! This method performs the generalized matrix-vector multiplication of
//    //! the form 
//    //!
//    //!     rvec = alpha* this * mvec + beta * rvec.
//    //! 
//    //! \param mvec source multiplication vector
//    //! \param rvec destination vector
//    //! \param alpha additional scalar factor for multiplication
//    //! \param beta scalar factor for re-use of \a rcvec
//    //! \param transposed if true, the transposed of the matrix is taken
//    //! \note If CFS is compiled without BLAS support, we use as fallback the
//    //!       internal matrix-vector multiplication.
//    void Mult_Blas( const Vector<TYPE> &mvec, Vector<TYPE> &rvec,
//                    Double alpha, Double beta, bool transposed) const;

    /** Perform a matrix-vector multiplication rvec = this*mvec via the Inner product.
     * Hence in the complex case this is the conjugate complex rvec = this*conj(mvec) */
    void MultInner( const SingleVector & mvec, SingleVector & rvec ) const;

    /** This implements the Frobenius inner product of two matrices. This is NOT the Frobenius norm!
     * @return the sum of the element wise product: sum this_ij * other_ij */
    TYPE FrobeniusProduct(const Matrix<TYPE>& other_mat) const;
    
    /** This implements the Frobenius norm of two matrices.
     * @return the sum of the element wise product: sum this_ij * other_ij */
    TYPE ScalarProduct(const Matrix<TYPE>& other_mat) const;

    //! Entry-wise multiplication with another matrix
    Matrix<double> EntryMult(const Matrix<double>& other_mat) const;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    void MultT( const SingleVector & mvec, SingleVector & rvec ) const;
  
    //! Perform a matrix-vector multiplication rvec += this*mvec
    void MultAdd( const SingleVector & mvec, SingleVector & rvec ) const
      { EXCEPTION("!!! Not implemented !!!" ); }
  
    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    void MultTAdd( const SingleVector & mvec, SingleVector& rvec ) const
      { EXCEPTION("!!! Not implemented !!!" ); }
  
    //! Perform a matrix-vector multiplication rvec -= this*mvec
    void MultSub( const SingleVector & mvec, SingleVector & rvec ) const
      { EXCEPTION("!!! Not implemented !!!" ); }

    //! Assign the matrix the dyadic product of a vector with itself
    
    //! Assigns the matrix itself the dyadic product of a vector vec1 
    //! with itself
    //!\param vec1 (input) Vector which gets multiplied with itself
    //!  \f[ \left( \begin{array}{ccc} m_{11} & m_{12} & \cdots \\
    //!  m_{21} & m_{22} & \cdots \\
    //!  \cdots & \cdots & \cdots 
    //!  \end{array} \right) 
    //!  =
    //!  \left( \begin{array}{c} v_1  \\ v_2 \\ \cdots \end{array} \right) 
    //!  \cdot
    //!  \left( \begin{array}{ccc} v_1 & v_2 & \cdots  \end{array} \right)
    //!  \f]
    inline void DyadicMult( const SingleVector & v1 )
    {
      DyadicMult(v1, v1);
    }
  
    //! Assign the matrix the dyadic product of two vectors

    //! Assigns the matrix itself the dyadic product of a vector vec1 
    //! with a vector vec2
    //! \param vec1 (input) Vector which gets multiplied with itself
    //! \f[ \left( \begin{array}{ccc} m_{11} & m_{12} & \cdots \\
    //! m_{21} & m_{22} & \cdots \\
    //! \cdots & \cdots & \cdots 
    //! \end{array} \right) 
    //!  =
    //!  \left( \begin{array}{c} v_1  \\ v_2 \\ \cdots \end{array} \right) 
    //!  \cdot
    //!  \left( \begin{array}{ccc} v_1 & v_2 & \cdots  \end{array} \right)
    //!  \f]
    void DyadicMult( const SingleVector & vec1, const SingleVector & vec2 );

    //! Calculate the Determinant (up to size 3)

    //! Calculates the determinant for a square matrix up to size 3x3.
    //! For larger matrices an error is returned.
    //! \param val (output) Return value of the method
    void Determinant( TYPE & val ) const;
    
    //! Calculates the Trace
    //! works for non-square matrices of any size
    TYPE Trace() const;

    /** Sum up the square of all entries */
    TYPE NormL2() const;

    /** Computes the average of all entries. Weaker than L1 norm */
    TYPE Avg() const;

    /** does something like (this - other).NormL2().
     * @see NormL2() */
    TYPE DiffNormL2(const Matrix<TYPE>& other) const;

    /** @see NormL2() but L1 norm */
    TYPE DiffNormL1(const Matrix<TYPE>& other) const;

    
    //@}

    //! Invert the matrix and store it in 'inv'
    
    //! This method calculates the inverse of the matrix and stores it
    //! into \a inv. The original matrix remains unchanged.
    //! This method is explicitly and efficient coded for matrices up 
    //! to size 3 x 3 and thus intended e.g. for Jacobian matrices.
    //! For larger ones it uses a factorization scheme.
    //! \param inv matrix which will hold the inverse
    void Invert ( Matrix <TYPE> & inv ) const;
    
    //! PseudoInvert the matrix and store it in 'inv'

    //! This method calculates the inverse of the matrix and stores it
    //! into \a inv. The original matrix remains unchanged.
    void PseudoInvert ( Matrix <TYPE> & inv ) const;

    //! Invert the matrix itself with Lapack
    
    //! This methods inverts a general matrix using a LU-factorization of
    //! Lapack. 
    //! \note The matrix itself gets overwritten in this method.
    void Invert_Lapack();
    
    //! Compute the condition number of the matrix with Lapack

    //! This methods estimates the condition number of the matrix
    //! \param k  estimated condition number
    //! \param info returns if the LU factorization was succesfull
    void Invert_Lapack(double & k, int & info);

    //! Transpose the matrix and store the result in \a transposedMat
    //! \note The matrix itself gets not changed.
    //! \note If the transposed of a matrix is needed for a operation
    //! with a vector, the according function like 'MultT' should be used
    void Transpose( Matrix<TYPE> & transposedMat ) const;
    
    /** Check if the matrix contains NAN. To be used by asserts() */
    bool ContainsNaN() const;

    /** Check if the matrix contains +/- INF. To be used by asserts() */
    bool ContainsInf() const;

#ifdef USE_EXPRESSION_TEMPLATES
    // =======================================================================
    // INTERFACE TO EXPRESSION TEMPLATES
    // =======================================================================
        
    //@{ 
    //! \name Interface To Expression Template Headers

    //! Matrix assignment operator using expression templates
    inline Matrix<TYPE>& operator=( const Matrix<TYPE>& rhs ) { 
      return this->assignFrom( rhs ); 
    }
    
    //! Scalar assignment operator using expression templates
    inline Matrix<TYPE>& operator=( TYPE rhs ) { 
      return this->assignFrom( rhs ); 
    }
    
    //! Matrix-Expression assignment operator using expression templates
    template <class X> inline Matrix<TYPE>& 
    operator=( const Xpr2<TYPE,X>& rhs ) {
      return this->assignFrom( rhs );
    }
    
    //! Abstract matrix assignment operator
    template <class M> inline Matrix<TYPE>& 
    operator=( const Dim2<TYPE,M>& rhs ) {
      return this->assignFrom(rhs);
    }
    
    
    //! Return number of rows
    inline unsigned int rows() const { return size_row_; }
    
    //! Return number of columns
    inline unsigned int cols() const { return size_col_; }
    
    //@}
#else
    // =======================================================================
    // MATHEMATICAL OPERATORS
    // =======================================================================

    //! \name Mathematical Operators
    //! \note Due to problems in Doxygen the binary operators +,-,*,/ 
    //!       (which use type promotion) are not shown, although they exist!
    //@{
    
    //! Assignment operator
    Matrix<TYPE> & operator=( const Matrix &y );
    
    //! Unary plus operator (this = +this)
    Matrix<TYPE> operator+() const;

    //! Create new matrix by addition (new = this + y)(type promotion)
    template <class TYPE2>
    Matrix<PROMOTE(TYPE,TYPE2)> operator+( const Matrix<TYPE2> &y ) const;

    //! Add a second matrix to own one (this += y)
    Matrix<TYPE> & operator+=(const Matrix<TYPE> &y );

    //! Unary minus operator (this = -this)
    Matrix<TYPE> operator-() const;

    //! Create new matrix by subtraction (new = this - y) (type promotion)
    template <class TYPE2> Matrix<PROMOTE(TYPE,TYPE2)> 
    operator-( const Matrix<TYPE2> &y ) const;

    //! Subtract a second matrix from own one (this -= y)
    Matrix<TYPE> & operator-=( const Matrix<TYPE> &y );

    //! Create new matrix by multiplication with scalar value 
    //!(type promotion)    
    template <class TYPE2>
    Matrix<PROMOTE(TYPE,TYPE2)> operator* ( const TYPE2 &y ) const;

    //! Create new vector by matrix-vector multiplication (type promotion)
    template <class TYPE2>
    Vector<PROMOTE(TYPE,TYPE2)> operator* ( const Vector<TYPE2> &y ) const;

    //! Create new matrix by matrix-matrix multiplication (type promotion)
    template <class TYPE2>
    Matrix<PROMOTE(TYPE,TYPE2)> operator*( const Matrix<TYPE2> &y ) const;

    //! Divide matrix by a scalar value
    Matrix<TYPE>  & operator/=( const TYPE y );

    //@}

#endif // USE_EXPRESSION_TEMPLATES
    
    //! Multiply matrix by a scalar (this *= y)
    Matrix<TYPE> & operator*=( const TYPE y );

    //! Perform matrix-matrix multiplication (this = this * arg)
    Matrix<TYPE> & operator*=( const Matrix<TYPE> &y );
    
    // =======================================================================
    // BOOLEAN OPERATORS
    // =======================================================================

    //! \name bool operators

    //@{
    
    //! Returns true if \a mat has the same entries as own matrix
    inline bool operator==( const Matrix<TYPE> & mat ) const;

    //! Returns true if \a mat has different entries than own matrix
    inline bool operator!=( const Matrix<TYPE> & mat ) const;
 
    //@}

    // =======================================================================
    // LAPACK INTERFACE
    // =======================================================================

    //! \name LAPACK Interface

    //@{
    //! Solves system of algebraic equation AX = B

    //! Solves system of algebraic equation AX=B
    //! where A is a quadratic matrix, and B a collection of 
    //! right hand side vectors which will be replaced by the 
    //! solution vectors. The enumeration LAPACK_MATRIX_TYPE
    //! describes the qualities of the system matrix A, 
    //! like symmetric, hermitian or general
    void solveWithLapack( Matrix<Complex> & b1, lapackSysMatType & LAPACK_MATRIX_TYPE );

    //! Computes eigenvalues of an hermitian matrix and eigen vectors if necessary
    void eigenvaluesWithLapack(Vector<Double> & b1,Matrix<double> * b2 = NULL);
    //@}
  
    // =======================================================================
    // MISCELLANEOUS METHODS
    // =======================================================================

    //! \name Miscellaneous methods

    //@{

    //! Solves a small system of equations (Ax=b) directly

    //! Solves directly a small system of equations of the form Ax=b
    //! using LU - decomposition (without pivoting!)
    //! \param x (output) solution vector      
    //! \param b (input) right-hand-side vector
    //! \param throw_exception throws an exception if we need to invert zero
    //! \return false only when not throw_exception and same issue
    //! \note The Matrix A=LU contains afterwards the the values of L 
    //! in the lower triangular, and the values of U in the upper part.
    bool DirectSolve(SingleVector& x, const SingleVector& b, bool throw_exception = true);

    /** variant of DirectSolve() which uses Lapack, based on cholesky decomposition, therefore only for S.P.D. matrices.
     * computes this * x = b
     * @param chol  will get the Cholesky decomposition (lower triangle) but the upper part not zeroed (see it as working array)
     * @param b rhs vector, could easily be extended for multiple rhs as the lapack kernel provides this */
    bool CholeskySolveLapack(Matrix<TYPE>& chol, Vector<TYPE>& x, const Vector<TYPE>& b, bool throw_exception = true);

    //! scales the diagonal elements of a  matrix by a factor
    inline void ScaleDiagElems(const TYPE factor);
    
    //! Return a special part ( real, imag, amplitude, phase) of a matrix
    Matrix<Double> GetPart(  Global::ComplexPart part ) const;

    //! Set special part ( real, imag, amplitude, phase) of a matrix
    
    //! This method explicitly set the real/imaginary part of a matrix.
    //! By default, the other part is left unchanged. If zeroOtherPart 
    //! is set to yes, the other part gets initialized to zero.
    void SetPart(Global::ComplexPart part, const Matrix<Double> & partMatrix, bool zeroOtherPart = false);

    //! S explicitly set the real/imaginary part of a matrix times a factor.

    //! By default, the other part is left unchanged. If zeroOtherPart
    //! is set to yes, the other part gets initialized to zero.
    void SetPartMult( Global::ComplexPart part,
                      const Matrix<Double> & partMatrix,
                      Double factor,
                      bool zeroOtherPart = false );

    //! Return a sub-part of the own matrix
    
    //! Copies a sub-matrix at the position (row, col) into subMat, 
    //! the amount of copied elements depends on the size of subMat
    void GetSubMatrix( Matrix<TYPE>& subMat, UInt row, UInt col ) const;

    //! Set a sub-part of the matrix
    
    //! Overwrites the matrix elements at the position (row, col) with subMat
    //! in a rectangular (submatrix) way
    void SetSubMatrix( const Matrix<TYPE>& subMat, UInt row, UInt col );

    //! Adds a subMat to the matrix elements at the position (row, col)
    //! in a rectangular (submatrix) way
    void AddSubMatrix( const Matrix<TYPE>& subMat, UInt row, UInt col );

    //! Converts a matrix into a vector, by appending successively all rows
    void ConvertToVec_AppendRows( SingleVector& vec ) const;

    //! Converts a matrix into a vector, by appending successively all cols
    void ConvertToVec_AppendCols( SingleVector& vec ) const;
 
    /** Converts the upper triangular of a quadratic matrix into a vector the way the stress and strain is defined for Voigt-Notation
     * 2*2 -> 11 22 12
     * 3*3 -> 11 22 33 23 13 12 */
    void ConvertToVec_UpperTriangular( SingleVector& vec ) const;

    /** Material notation. Only for FMO we assume the design to be Hill-Mandel, in LinElastInt we use Voigt. The CFS-B-operator is also Voigt, _NO_DENSITY sets topology variable to 1 in simultaneous material and top. opt. */

    //! Only for testing the switching state of Preisach planes
    void matrix2Bmp(UInt upscale, std::string filename,Matrix<TYPE>* greenChannel = NULL);
    void matrix2Bmp_v2(UInt upscale, std::string filename,Matrix<TYPE>* rotX, Matrix<TYPE>* rotY);
    void matrix2Bmp_v3(UInt upscale, std::string filename,Matrix<TYPE>* rotX, Matrix<TYPE>* rotY);

    /** print content for logging or nice output.
     * See unittests.cc Matrix_ToString for test cases
     * TS_PLAIN is space separated with line breaks for rows
     * TS_MATLAB is with brackets, no comma between elements and line break for rows
     * TS_PYTHON is with brackets, commas and no line break for rows as default
     * TS_INFO gives summary
     * TS_NONZEROS names coordinate as in matrix market
     * @param linesep separator string, when empty meaningful default is used
     * @param digits when given enforces format */
    std::string ToString(ToStringFormat format = TS_PYTHON, const std::string& linesep="", int digits=-1) const override;

    /** Creates a xml string of the following form.
     * <name dim1="6" dim2="6">
         <real>
           1.65682E+11 1.84091E+10 1.84091E+10 0.00000E+00 0.00000E+00 0.00000E+00
           1.84091E+10 1.65682E+11 1.84091E+10 0.00000E+00 0.00000E+00 0.00000E+00
           1.84091E+10 1.84091E+10 1.65682E+11 0.00000E+00 0.00000E+00 0.00000E+00
           0.00000E+00 0.00000E+00 0.00000E+00 7.36364E+10 0.00000E+00 0.00000E+00
           0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 7.36364E+10 0.00000E+00
           0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 7.36364E+10
         </real>
       </name>
       The ident are two spaces + offset.
       Such stuff is read from the material file
      @param name the name for above, the label form the parent ParamNode, in general 'tensor'!
      @param offset spaces in front of line */
    std::string ToXMLFormat(const std::string& name, const int n_offset) const override;

    /** Parses a string generated by ToString(0).
        Reconstructs the dimension in that way but you should know the propert type */
    void Parse(const std::string& data);
    
    //! Rotate this matrix by a given rotation matrix
    
    //! This method generates a copy of this matrix, which contains the rotated
    //! content, as defined by the rotation matrix rotMatrix.
    //! \note This method will only work with matrices of size 2,3, and 6.
    void PerformRotation( const Matrix<Double>& rotMatrix,  Matrix<TYPE>& matMatrix ) const;

    //@}

  private:
    /** Helper method for Parse() */
    unsigned int ParseLineHelper(const std::string& input, StdVector<TYPE>& out);

    /** Helper for MultBLAS() */
    void CallGEMM(char* transa, char* transb, int* m, int* n, int* k, TYPE* alpha, TYPE* a, int* lda, TYPE* b, int* ldb, TYPE* beta, TYPE* c, int* ldc) const;

    //! Calculates the adjunct of the matrix at position (i,j)
    TYPE Adjunct (UInt i, UInt j) const;

    //! Number of rows 
    UInt size_row_ = 0;
  
    //! Number of columns
    UInt size_col_ = 0;

    //! Data of the matrix
    TYPE** data_ = nullptr;

    /** keeps the actual data with single allocation and rows and data as 
     * continuous memory for better cache performance */
    boost::container::small_vector<std::byte, CalcBufferSizeBytes(6,6)> buffer_;
    // one can replace small_vector by std::vector<std::byte> buffer_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================
  //! \class Matrix
  //! 
  //! \purpose This class implements a general, templatized dense matrix with 
  //! the following additional features:
  //! - If the macro USE_EXPR_TEMPLATES is defined, it utilizes an expression
  //! template library, which evaluates mathematical expression at compile 
  //! time.
  //! In this case, the operators(+,-,*,/,-=,+=,...) do not have to be defined
  //! in this class, but are evaluated by the expression template library.
  //! The used library is a modified version of the MET-library 
  //! (<a href="http://met.sourceforge.net">met.sourceforge.net</a>).
  //! - In order to be able to handle mixed Double-Complex valued mathematical
  //! expressions, the concept of Type Promotion / Traits is utilized (see also
  //! <a href="http://osl.iu.edu/~tveldhui/papers/techniques/"> Techniques for
  //! Scientific C++ </a>). This means, expressions like
  //! \verbatim
  //! Matrix<Double> realMat1;
  //! Matrix<Complex> complexMat1, complexMat2;
  //! Vector<Double> realVec1;
  //! Vector<Double> complexVec1, complexVec2;
  //! Double realFactor = 1.0;
  //! Complex complexFactor = Complex(1.0, 1.0);
  //!
  //! complexMat1 = realMat1 * complexFactor;
  //! complexMat2 = complexMat1 * realFactor;
  //! complexMat1 = complexMat2 + realMat1;
  //!
  //! complexVec1 = realMat1 * complexVec2;
  //! complexVec2 = complexMat1 * realVec1;
  //! \endverbatim
  //! can be written and the conversion is done automatically.
  //! \note - Multiple Double <-> Complex conversion in one statement are
  //!       not possible!
  //! 
  //! \note -If expression templates are used, statements like
  //! \verbatim
  //! Matrix<Double> mat = mat1 * 5.0;
  //! \endverbatim 
  //! have to be replaced by
  //! \verbatim
  //! Matrix<Double> mat;
  //! Double factor = 5.0;
  //! mat = mat1 * factor;
  //! \endverbatim 
  //! 
  //! \collab The Matrix class can be used together with the templatized Vector
  //! class.
  //! 
  //! \implement This class uses the concept of type promotion / traits and
  //! can additionally utilize expression templates.
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve 
  //! - Check 'const'-correctness of class!
  //! - Add safety checks for initialization
  //! 

#endif
  


  // =======================================================================
  // RELATED FUNCTIONS
  // =======================================================================
  //! \relates Matrix
  //! Output operator for std::ostream
  template<class TYPE>  std::ostream& operator << ( std::ostream & , 
                                                    const Matrix<TYPE> &);
  
#ifndef USE_EXPRESSION_TEMPLATES
  //! Explicit Transpose function
  template<class TYPE>
  Matrix<TYPE> Transpose( const Matrix<TYPE>& m ) {
    UInt numRows = m.GetNumRows();
    UInt numCols = m.GetNumCols();
    Matrix<TYPE> trans(numCols, numRows);
    for( UInt i = 0; i < numCols; i++ ) {
      for (UInt j = 0; j < numRows; j++ ) {
        trans[i][j] = m[j][i];
      }
    }
    return trans;
  }
  
  //! Explicit conjugation of matrix
   template<class TYPE>
   Matrix<TYPE> Conj( const Matrix<TYPE>&m);
  
  //! Explicit Hermitian of matrix
  template<class TYPE>
  Matrix<TYPE> Herm( const Matrix<TYPE>&m);

  #endif

  //! Explicit Transpose function
  template<class TYPE>
  Matrix<TYPE> TransposeConjugate( const Matrix<TYPE>& m )
  {
    Matrix<TYPE> trans(m.GetNumCols(), m.GetNumRows());

    for( UInt i = 0, in = m.GetNumCols(); i < in; i++ )
      for (UInt j = 0, jn = m.GetNumRows(); j < jn; j++ )
        trans[i][j] = conj(m[j][i]);

    return trans;
  }



  // =======================================================================
  // INLINE MEMBER IMPLEMENTATION
  // =======================================================================

  template<class TYPE>
  inline void Matrix<TYPE>::Determinant (TYPE & ret) const
  {       
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0|| size_col_ == 0)
      EXCEPTION( "Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_col_ ) 
      EXCEPTION( "No quadratic matrix!" );
#endif
  
    switch (size_row_)
      {
      case 1: ret =  data_[0][0];
        break;
      case 2: ret =   data_[0][0]*data_[1][1]-data_[0][1]*data_[1][0];
        break;
      case 3: ret = data_[0][0]*data_[1][1]*data_[2][2] +
          data_[0][1]*data_[1][2]*data_[2][0] +
          data_[0][2]*data_[1][0]*data_[2][1] -
          data_[0][2]*data_[1][1]*data_[2][0] -
          data_[0][1]*data_[1][0]*data_[2][2] -
          data_[0][0]*data_[1][2]*data_[2][1];
        break;
      default: 
        EXCEPTION( "Dimension larger than 3!" );
      }
  }

  template<class TYPE>
  inline TYPE Matrix<TYPE>::Trace() const {
    assert(!(size_row_ == 0|| size_col_ == 0));
    UInt smallersize = size_row_ < size_col_ ? size_row_ : size_col_;
    TYPE ret = data_[0][0];
    for(UInt i = 1; i < smallersize; i++)
      ret += data_[i][i];
    return ret;
  }

  // Perform a matrix-matrix multiplication rMat = this*mMat
  template<class TYPE>
  inline void Matrix<TYPE>::Mult(const DenseMatrix & mMat, DenseMatrix & rMat) const
  {
    Matrix<TYPE> const & mMat1 = dynamic_cast<const Matrix<TYPE>& >(mMat);
    Matrix<TYPE> & rMat1 = dynamic_cast<Matrix<TYPE>& >(rMat);
  
    UInt size_mMatRow = mMat1.GetNumRows();
    UInt size_mMatCol = mMat1.GetNumCols();

#ifdef CHECK_INITIALIZED
    UInt size_rMatRow = rMat1.GetNumRows();
    UInt size_rMatCol = rMat1.GetNumCols();

    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mMatRow == 0 || size_mMatCol==0) 
      EXCEPTION("undefined Matrix");
    if (size_rMatRow == 0||size_rMatCol==0) 
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mMatRow) {
      EXCEPTION("incompatible dimension while matrix-matrix multiplication" );
    }
    if (size_row_ != size_rMatRow) {
      EXCEPTION("incompatible dimension while matrix-matrix multiplication" );
    }
    if (size_mMatCol != size_rMatCol) {
      EXCEPTION("incompatibel dimension while matrix-matrix multiplication" );
    }
#endif
   
    for (UInt i = 0; i < size_row_; i++ ) {
      for (UInt j = 0; j < size_mMatCol; j++ ) {
        rMat1[i][j] = data_[i][0] * mMat1[0][j];
        for ( UInt k = 1; k < size_mMatRow; k++ ) {
          rMat1[i][j] += data_[i][k] * mMat1[k][j];
        }
      }
    }
  }

  // Perform a matrix-matrix multiplication rMat = Transpose(this)*mMat
  template<class TYPE>
  inline void Matrix<TYPE>::MultT(const DenseMatrix & mMat, 
                                 DenseMatrix & rMat) const {


    Matrix<TYPE> const & mMat1 = dynamic_cast<const Matrix<TYPE>& >(mMat);
    Matrix<TYPE> & rMat1 = dynamic_cast<Matrix<TYPE>& >(rMat);
  
    UInt size_mMatRow = mMat1.GetNumRows();
    UInt size_mMatCol = mMat1.GetNumCols();

#ifdef CHECK_INITIALIZED
    UInt size_rMatRow = rMat1.GetNumRows();
    UInt size_rMatCol = rMat1.GetNumCols();

    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mMatRow == 0 || size_mMatCol==0) 
      EXCEPTION("undefined Matrix");
    if (size_rMatRow == 0||size_rMatCol==0) 
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_mMatRow) {
      EXCEPTION("incompatible dimension while matrix-matrix multiplication" );
    }
    if (size_col_ != size_rMatRow) {
      EXCEPTION("incompatible dimension while matrix-matrix multiplication" );
    }
    if (size_mMatCol != size_rMatCol) {
      EXCEPTION("incompatibel dimension while matrix-matrix multiplication" );
    }
#endif
   
    for (UInt i = 0; i < size_col_; i++ ) {
      for (UInt j = 0; j < size_mMatCol; j++ ) {
        rMat1[i][j] = data_[0][i] * mMat1[0][j];
        for ( UInt k = 1; k < size_mMatRow; k++ ) {
          rMat1[i][j] += data_[k][i] * mMat1[k][j];
        }
      }
    }
  }

  // =======================================================================
  //  Inline part for all operators using type promotion
  //  rr being defined only in non-template-expression case
  // =======================================================================

#ifndef USE_EXPRESSION_TEMPLATES

  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::operator+(const Matrix<TYPE2> &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif 
    
#ifdef CHECK_INDEX
    if (size_row_ != x.GetNumRows() || size_col_ != x.GetNumCols())
      EXCEPTION("incompatible dimension");
#endif
  
    Matrix<PROMOTE(TYPE,TYPE2)> z(size_row_,size_col_);
    
    for(UInt k = 0; k < size_row_*size_col_; ++k)                                                                                
      z [0][k] = x[0][k]+data_[0][k];
    
    return z;
  }

  
  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::operator-(const Matrix<TYPE2> &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.GetNumRows() == 0 || x.GetNumCols() == 0)
      EXCEPTION("undefined Matrix");
#endif
  
#ifdef CHECK_INDEX  
    if (size_row_ != x.GetNumRows() || size_col_ != x.GetNumCols())
      EXCEPTION("incompatible dimension for +"); 
#endif
  
    Matrix<PROMOTE(TYPE,TYPE2)> z(size_row_,size_col_);
    
    for(UInt k = 0; k < size_row_*size_col_; ++k)
      z[0][k] = -x[0][k]+data_[0][k];

    return z;
  }

  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::operator*(const TYPE2 &x) const 
  { 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif
  
    Matrix<PROMOTE(TYPE,TYPE2)> z(size_row_,size_col_);
  
    for(UInt k = 0; k < size_row_*size_col_; ++k)
      z [0][k] = data_[0][k]*x;
    
    return z;
  }

  template<class TYPE>  template<class TYPE2>
  Vector<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::operator*(const Vector<TYPE2> &x) const
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (x.GetSize() == 0) EXCEPTION("undefined Vector");
#endif

#ifdef CHECK_INDEX
    if (size_col_ != x.GetSize()) 
    {
	    std::cout << "Matrix r x c " << size_row_ << " x " << size_col_ << std::endl;
	    std::cout << size_col_ << " vs " << x.GetSize() << std::endl;
      EXCEPTION("incompatible dimension");
    }
#endif
  
    Vector<PROMOTE(TYPE,TYPE2)> z(size_row_);
  
    UInt k,kk;
    for ( k = 0; k < size_row_; k++)
      for ( kk = 0; kk < size_col_; kk++)
        z[k] += data_[k][kk] * x[kk];

    return z;
  }

  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::
  operator*(const Matrix<TYPE2> &x) const
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.GetNumRows() == 0 || x.GetNumCols()== 0)
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX  
    if (size_col_ != x.GetNumRows())
      EXCEPTION("incompatible dimension");
#endif
 
    PROMOTE(TYPE,TYPE2) a;
    Matrix<PROMOTE(TYPE,TYPE2)>  z (size_row_, x.GetNumCols());
  
    for (UInt i = 0; i < size_row_; i++)
    {
      for (UInt j = 0, cols = x.GetNumCols(); j < cols; j++)
        {       
          a = data_ [i] [0] * x[0][j];
          for (UInt k = 1; k < size_col_; k++)
            a += data_ [i] [k] * x[k][j];
          z(i,j) = a;
        }
    }
    return z;
  }
#endif //USE_EXPRESSION_TEMPLATES

} //end of namespace


#endif  // FILE_MATRIX

