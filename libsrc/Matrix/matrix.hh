#ifndef FILE_MATRIX_2004
#define FILE_MATRIX_2004

#include "cfsmatrix.hh"
#ifdef USE_LAPACK
#include "matrixLapackSupport.hh"
#endif

namespace CoupledField
{      

  //! Overloading << for class Matrix
  template<class TYPE> class Vector;

  //! Concrete implementation of a dense matrix
  template<class TYPE>
  class Matrix: public CFSMatrix
  {
  public:

    // Friend declarations
    friend class Vector<TYPE>;
  
    //! Constructor 
    //! creates an empty matrix of size 0x0
    Matrix();

    //! Constructor
    //! creates a matrix of size nRows x nCols, initalized
    //! with zeroes
    /*! 
      \param nRows (input) Number of rows
      \param nCols (input) Number of columns)
    */
    Matrix (const UInt nRows, const UInt nCols);

    //! Constructor (number of vectors, array of vectors(colomn))
    Matrix (const UInt, const Vector<TYPE> * const);

    //! Default Copy Construcctor
    Matrix(const Matrix &);

    //! Special copy constructor (e.g. convert double to complex)
    template<class T2> 
    Matrix(const Matrix<T2> &);

    //! Destructor
    ~Matrix();

    //! Hard coded query if values are complex
    Boolean IsComplex() const;

    //! Initialize matrix with a given entry.
    //! If no entry given, it gets initalized with zeroes
    /*!
      \param val (input,opt.) Entry the matrix gets initalized with
    */
    //! \note This method does not change the size of the matrix
    void Init(const TYPE val = TYPE());

    //! Change size of quadratic matrix
    /*!
      \param size (input) Number of rows / columns
    */
    //! \note the matrix contains afterwards only zeroes
    void Resize(const UInt size);
  
    //! Change size of general matrix 
    /*!
      \param nRows (input) Number of rows
      \param nCols (input) Number of columns)
    */
    //! \note the matrix contains afterwards only zeroes
    void Resize(const UInt nRows, const UInt nCols);

    //! Get the number of rows
    UInt GetSizeRow() const;
  
    //! Get the number of columns
    UInt GetSizeCol() const;
  
    //! Set the entry 'val' at position (row,col) in the matrix
    /*!
      \param row (input) Row of entry
      \param col (input) Column of entry
      \param val (input) Value to be set
    */
    void SetEntry( const UInt row, const UInt col, const TYPE val ) {
      data_[row][col] = val;
    }
  
    //! Add'val' to the matrix entry at position (row,col) in the matrix
    /*!
      \param row (input) Row of entry
      \param col (input) Column of entry
      \param val (input) Value to be added
    */
    void AddToEntry(const UInt row, const UInt col, const TYPE val);
   
    //! Get the entry 'val' at position (row,col) in the matrix

    //! \param row (input) row index of entry
    //! \param col (input) column index of entry
    //! \param val (output) on return contains value of entry
    void GetEntry(const UInt row, const UInt col, TYPE & val) const {
      val = *( data_[0] + row * size_col_ + col ); 
    }

    //! Calculates the determinant (up to size 3)
    /*!
      \param val (output) Return value of the method
    */
    void Determinant(TYPE & val) const;

    //! Invert the matrix and store it in 'inv' (up to size 3)
    /*!
      \param inv (output) Inverse of the matrix
    */
    void Invert(CFSMatrix & inv) const {
      Error( "!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    };
  
    //! Transpose the matrix and store it in 'trans'
    /*!
      \param trans (output) Transposed matrix
    */
    //! \note The matrix itself gets not changed.
    //! \note If the transposed of a matrix is needed for a operation
    //! with a vector, the according function like 'MultT' should be used
    void Transpose(CFSMatrix & trans) const {
      Error("!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    };


    //! Solves a small system of equations (Ax=b) directly
    /*!
      \param x (output) solution vector
      \param b (input) right-hand-side vector
    */
    //! Solves directly a small system of equations of the form Ax=b
    //! using LU - decomposition (without pivoting!)
    //! \note The Matrix A=LU contains afterwards the the values of L 
    //! in the lower triangular, and the values of U in the upper part.
    void DirectSolve(CFSVector & x, CFSVector & b);

#ifdef USE_LAPACK
    //! Solves system of algebraic equation AX=B
    //! where A is a quadratic matrix, and B a collection of 
    //! right hans side vectors which will be replaced by the 
    //! solution vectors. The enumeration LAPACK_MATRIX_TYPE
    //! describes the qualities of the system matrix A, 
    //! like symmetric, hermitian or general
    //! Compile with LAPACK - Support (USE_LAPACK = yes)
    void solveWithLapack(Matrix<Complex> & b1,
                         lapackSysMatType & LAPACK_MATRIX_TYPE);

    //! Computes eigenvalues of an hermitian matrix
    void eigenvaluesWithLapack(Vector<Double> & b1);

    //! Converts a fortran 77 matrix to C++ complex
    void F772CC( const F77complex16 &v, std::complex<double> &val ) {
      std::complex<double> aux(v.real,v.imag);
      val = aux;
    }

    void F772CC( const F77real8 &v, double &val ) {
      val = (double)v;
    }

    //! Converts cfs data to fortran 77 format
    void CC2F77( const std::complex<double> &v, F77complex16 &val ) {
      val.real = (F77real8)v.real();
      val.imag = (F77real8)v.imag();
    }
    void CC2F77( const double &v, F77real8 &val ) {
      val = (F77real8)v;
    }
#endif
  
    //! Assignes the matrix itself the dyadic product of a vector vec1 
    //! with itself
    /*!
      \param vec1 (input) Vector which gets multiplied with itself
      \f[ \left( \begin{array}{ccc} m_{11} & m_{12} & \cdots \\ 
      m_{21} & m_{22} & \cdots \\
      \cdots & \cdots & \cdots 
      \end{array} \right) 
      =
      \left( \begin{array}{c} v_1  \\ v_2 \\ \cdots \end{array} \right) 
      \cdot
      \left( \begin{array}{ccc} v_1 & v_2 & \cdots  \end{array} \right)
      \f]
    */                    
    void DyadicMult(const CFSVector & vec1);  
  
    //! Assignes the matrix itself the dyadic product of a vector vec1 
    //! with a vector vec2
    /*!
      \param vec1 (input) Vector which gets multiplied with itself
      \f[ \left( \begin{array}{ccc} m_{11} & m_{12} & \cdots \\ 
      m_{21} & m_{22} & \cdots \\
      \cdots & \cdots & \cdots 
      \end{array} \right) 
      =
      \left( \begin{array}{c} v_1  \\ v_2 \\ \cdots \end{array} \right) 
      \cdot
      \left( \begin{array}{ccc} v_1 & v_2 & \cdots  \end{array} \right)
      \f]
    */          
    void DyadicMult(const CFSVector & vec1, const CFSVector & vec2); 
  
    //! copies a submatrix at the position (row, col) into subMat, 
    //! the amount of copied elements depends on the size of subMat
    void GetSubMatrix( CFSMatrix &subMat, const UInt nRows,
                       const UInt nCols ) const {
      Error( "!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    };
  
    //! overwrites the matrix elements at the position (row, col) with subMat
    //! in a rectangular (submatrix) way
    void SetSubMatrix(const CFSMatrix & subMat, const UInt nRows,
                      const UInt nCols) {
      Error("!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    };
  
    //! scales the diagonal elements of a  matrix by a factor
    void ScaleDiagElems(const TYPE factor);

    //!
    void Add(const TYPE fac, const CFSMatrix & mat){};

    //! Perform a matrix-matrix multiplication rMat = this*mMat
    void Mult(const CFSMatrix & mMat, CFSMatrix & rMat) const ;
  
    //! Perform a matrix-vector multiplication rvec = this*mvec
    void Mult(const CFSVector & mvec, CFSVector & rvec) const;


    //! Perform a matrix(Double)-vector(Complex) multiplication 
    //! rvec = this*mvec where the matrix is supposed to be of
    //! type Double, rvec and mvec are complex valued
    void MatVecMult_DC(const Vector<Complex> & mvec, 
                       Vector<Complex> & rvec) const;

    //! Perform a matrix(Complex)-vector(Double) multiplication
    //! rvec = this*mvec where the matrix is supposed to be of
    //! type Complex as well as rvec; mvec is of type Double
    void MatVecMult_CD(const Vector<Double> & mvec, 
                       Vector<Complex> & rvec) const;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    void MultT(const CFSVector & mvec, CFSVector & rvec) const {;};
  
    //! Perform a matrix-vector multiplication rvec += this*mvec
    void MultAdd(const CFSVector & mvec, CFSVector & rvec)const {;};
  
    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    void MultTAdd(const CFSVector & mvec, CFSVector& rvec) const {;};
  
    //! Perform a matrix-vector multiplication rvec -= this*mvec
    void MultSub(const CFSVector & mvec, CFSVector & rvec) const {;};

    //! Check if the matrix is symmetric
    bool IsSymmetric() const;

    //////////////////////////////////////
    // Functions for working with other //
    // Matrix<Type> and Vector<TYPE>    //
    //////////////////////////////////////

    //! Assignment operator
    Matrix<TYPE> & operator= (const Matrix &);

    //! Return pointer to row number []
    inline TYPE * operator[] (const UInt) const;

    //! fast inversion for matrices smaller than size 3
    void Invert (Matrix <TYPE> & inv) const;

    //! returns pointer to continuos chunck of data
    TYPE ** GetRowPointer() const;

    //! return pointer to array of elements in matrix, row by row
    inline TYPE * GetDataPointer() const { return data_[0];}

    //! Overloading of operations

    //! Access operator
    /*!
      \param row (Input) Row number
      \param col (Input) Column number
    */
    TYPE & operator()(const UInt row , const UInt col);

    //! 
    Matrix<TYPE> operator+() const;

    //!
    Matrix<TYPE> operator+(const Matrix<TYPE> &) const;

    //!
    Matrix<TYPE> & operator+=(const Matrix<TYPE> &);

    //!
    Matrix<TYPE> operator-() const;

    //!
    Matrix<TYPE> operator-(const Matrix<TYPE> &) const;

    //!
    Matrix<TYPE> & operator-=(const Matrix<TYPE> &);

    //! multiplication with scalar value
    Matrix<TYPE> operator* (const TYPE &) const;

    //!
    Vector<TYPE> operator* (const Vector<TYPE> &) const;

    //!
    Matrix<TYPE> operator*(const Matrix<TYPE> &) const;

    //!
    Matrix<TYPE> & operator*=(const TYPE &);


    //!
    Matrix<TYPE> & operator*=(const Matrix<TYPE> &);

    //!
    Matrix<TYPE> & operator/=(const TYPE &);
    //!
    Boolean operator ==(const Matrix<TYPE> &) const;

    //!
    Boolean operator!=(const Matrix<TYPE> &) const;
 

    //  //! Cut part of matrix (left index row, right, upper index col, low )
    //   //   Matrix    part    (const UInt, const UInt,
    //   //                          const UInt, const UInt) const;
  
    //   //! Cut row number i, colomn number j from matrix
    //   //void cut(const UInt i, const UInt j);

    //! Add a row to Matrix at position i
    void AddRow(const Vector<TYPE> & x, const UInt pos );

    //! Add a colomn to Matrix at position i
    void AddColumn(const Vector<TYPE> & x, const UInt pos ); 

    /// Transpose actual matrix
    void Transpose (Matrix<TYPE> &transposedMat) const;  

    /// copies a submatrix at the position (row, col) into subMat, 
    /// the amount of copied elements depends on the size of subMat
    void GetSubMatrix(Matrix<TYPE>& subMat, UInt row, UInt col) const;

    /// overwrites the matrix elements at the position (row, col) with subMat
    /// in a rectangular (submatrix) way
    void SetSubMatrix(const Matrix<TYPE>& subMat, UInt row, UInt col);

    /// converts a matrix into a vector, by appending successively all rows
    void ConvertToVec_AppendRows(CFSVector& vec) const;

    /// converts a matrix into a vector, by appending successively all cols
    void ConvertToVec_AppendCols(CFSVector& vec) const;

    /// gets the diagonal elements of a  matrix in a one column matrix
    void GetDiagInMatrix(Matrix<TYPE>& columnMat);
    //Matrix<TYPE> GetDiagInMatrix();

  private:

    //! calculates the adjunct of the matrix at position (i,j)
    TYPE Adjunct (UInt i, UInt j) const;

    //! number of rows 
    UInt size_row_;
  
    //! number of columns
    UInt size_col_;

    //! data of the matrix
    TYPE** data_;


  };

  /////////////////////////////
  // Inline member functions //
  /////////////////////////////


  template<class TYPE>
  inline void Matrix<TYPE>::Init(const TYPE val)
  {
    UInt i;
    for (i=0; i<size_row_*size_col_; i++) 
      data_[0][i]=val;
  }

  template<class TYPE>
  inline void Matrix<TYPE>::AddToEntry ( const UInt i, const UInt j,
                                         const TYPE value ) {
    data_[i][j]+=value;
  }

  template<class TYPE>
  inline TYPE & Matrix<TYPE>::operator()(const UInt row, const UInt col) 
  {
    return data_[row][col];
  }

  template<class TYPE>
  inline TYPE * Matrix<TYPE>::operator[] (const UInt i) const
  { 
    ENTER_IFCN("Matrix::operator[]");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) Error("undefined Matrix",
                                                __FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (i < 0 || i >= size_row_) Error("invalid index",__FILE__,__LINE__);
#endif
  
    return data_[i];
  }

  template<class TYPE>
  inline TYPE ** Matrix<TYPE>::GetRowPointer () const
  {
    return data_;
  }


  template<class TYPE>
  inline UInt Matrix<TYPE>::GetSizeRow() const
  {       
    return size_row_;
  }
 
  template<class TYPE>
  inline UInt Matrix<TYPE>::GetSizeCol () const
  {       
    return size_col_;
  }

  template<class TYPE>
  inline void Matrix<TYPE>::Determinant (TYPE & ret) const
  {       
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0|| size_col_ == 0) 
      Error("Undefined Matrix!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_col_ ) 
      Error("No quadratic matrix!",__FILE__,__LINE__);
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
      default: Error("Dimension larger than 3!",__FILE__,__LINE__);
      }
  }

  template<class TYPE>  std::ostream& operator << ( std::ostream & , 
                                                    const Matrix<TYPE> &);



  // Perform a matrix-matrix multiplication rMat = this*mMat
  template<class TYPE>
  inline void Matrix<TYPE>::Mult(const CFSMatrix & mMat, 
                                 CFSMatrix & rMat) const {

    ENTER_IFCN( "Matrix::Mult" );

    Matrix<TYPE> const & mMat1 = dynamic_cast<const Matrix<TYPE>& >(mMat);
    Matrix<TYPE> & rMat1 = dynamic_cast<Matrix<TYPE>& >(rMat);
  
    UInt size_mMatRow = mMat1.GetSizeRow();
    UInt size_mMatCol = mMat1.GetSizeCol();

#ifdef CHECK_INITIALIZED
    UInt size_rMatRow = rMat1.GetSizeRow();
    UInt size_rMatCol = rMat1.GetSizeCol();

    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
    if (size_mMatRow == 0 || size_mMatCol==0) 
      Error("undefined Matrix",__FILE__,__LINE__);
    if (size_rMatRow == 0||size_rMatCol==0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mMatRow) {
      Error("incompatible dimension while matrix-matrix multiplication",
            __FILE__,__LINE__);
    }
    if (size_row_ != size_rMatRow) {
      Error("incompatible dimension while matrix-matrix multiplication",
            __FILE__,__LINE__);
    }
    if (size_mMatCol != size_rMatCol) {
      Error("incompatibel dimension while matrix-matrix multiplication",
            __FILE__,__LINE__);
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

#if defined(__GNUC__) 
  template class Matrix<Double>;
  template class Matrix<Integer>;
  template class Matrix<UInt>;
  template class Matrix<Complex>;
#endif

} //end of namespace


#endif  // FILE_MATRIX

