#ifndef FILE_MATRIX_2004
#define FILE_MATRIX_2004

#include "Utils/boost-serialization.hh"
#include "cfsmatrix.hh"
#include "Utils/promote.hh"

#ifdef USE_LAPACK
#include "matrixLapackSupport.hh"
#endif

#ifdef EXPR_TEMPLATES
#include "Utils/exprt/xpr2.hh"
#endif

namespace CoupledField
{      

  //! Forward class declaration
  template<class TYPE> class Vector;
  

  //! Concrete implementation of a dense matrix
  template<class TYPE>
#ifdef EXPR_TEMPLATES
  class Matrix: public CFSMatrix, public Dim2<TYPE, Matrix<TYPE> >
#else
  class Matrix: public CFSMatrix 
#endif
  {
  public:
    
    //! Friend declaration for vector
    friend class Vector<TYPE>;

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION, INITIALIZATION, RESIZING
    // =======================================================================
    
    //! \name Construction, Destruction, Initialization and Resizing

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
    void Init( const TYPE val = TYPE() );

    //! Change the size of the matrix

    //! Change size of general matrix 
    //! \param nRows (input) Number of rows
    //! \param nCols (input) Number of columns
    void Resize(const UInt nRows, const UInt nCols );

    //! Changes the size so that the matrix gets quadratic
    
    //! Changes the size of the matrix according to \a size.
    //! \param size (input) Number of rows / columns
    void Resize( const UInt size );

    //@}
    
    // =======================================================================
    // GENERAL INFORMATION
    // =======================================================================
    
    //! \name General Matrix Information
    
    //@{

    //! Get entry type of matrix
    EntryType::ScalarType GetEntryType() {
      return  EntryTypeMap<TYPE>::S_TYPE;
    }

    //! Check if the matrix is symmetric

    //! Return true, if the matrix is symmetric
    //! \note The results might be incorrect due to numeric rounding errors
    bool IsSymmetric() const;

    //! Get the number of rows
    UInt GetSizeRow() const;
    
    //! Get the number of columns
    UInt GetSizeCol() const;

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
    inline TYPE & operator()( UInt row, UInt col) {
      return data_[row][col];
    }

    //! General access operator (const)

    //! Access operator of one entries
    //! \param row (input) Row number
    //! \param col (input) Column number
    inline TYPE operator()( UInt row, UInt col ) const {
      return data_[row][col];
    }
    
    //! Returns pointer to row \a row
    inline TYPE * operator[]( const UInt row ) const;

    //! Returns pointer to raw data

    //! Returns pointer to continuous chunk of data
    inline TYPE ** GetRowPointer() const;

    //! Returns pointer to array of elements in matrix, row by row
    inline TYPE * GetDataPointer() const { return data_[0];}
    
    //! Get the entry 'val' at position (row,col) in the matrix
    
    //! Return entry at position (\a row, \a col) in the matrix
    //! \param row (input) row index of entry
    //! \param col (input) column index of entry
    //! \param val (output) on return contains value of entry
    inline void GetEntry( const UInt row, const UInt col, 
                          TYPE & val ) const {
      val = *( data_[0] + row * size_col_ + col ); 
    }
     
    //! Set the entry 'val' at position (\a row, \a col) in the matrix
    
    //! Set the entry 'val' at position (\a row, \a col) in the matrix
    //! \param row (input) Row of entry
    //! \param col (input) Column of entry
    //! \param val (input) Value to be set
    inline void SetEntry( const UInt row, const UInt col, const TYPE val ) {
      data_[row][col] = val;
    }

    //! Add'val' to the matrix entry at position (row,col) in the matrix
    
    //! Add'val' to the matrix entry at position (\a row, \a col) in the 
    //! matrix
    //! \param row (input) Row of entry
    //! \param col (input) Column of entry
    //! \param val (input) Value to be added
    void AddToEntry( const UInt row, const UInt col, const TYPE val );
    
    //! Gets the diagonal elements of a  matrix in a one column matrix
    void GetDiagInMatrix( Matrix<TYPE>& columnMat ) const;

    //@}

    // =======================================================================
    // NAMED ARITHMETIC OPERATIONS
    // =======================================================================
    
    //! \name Named Arithmetic Operations
    //@{

    //! Add the multiple of another matrix this = fac * mat
    void Add( const TYPE fac, const CFSMatrix & mat) {
      Error( "Not implemented yet!", __FILE__, __LINE__ );
    }
    
    //! Perform a matrix-matrix multiplication rMat = this*mMat
    void Mult(const CFSMatrix & mMat, CFSMatrix & rMat) const;

    //! Perform a matrix-vector multiplication rvec = this*mvec
    void Mult( const CFSVector & mvec, CFSVector & rvec ) const;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    void MultT( const CFSVector & mvec, CFSVector & rvec ) const {;};
  
    //! Perform a matrix-vector multiplication rvec += this*mvec
    void MultAdd( const CFSVector & mvec, CFSVector & rvec ) const {;};
  
    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    void MultTAdd( const CFSVector & mvec, CFSVector& rvec ) const {;};
  
    //! Perform a matrix-vector multiplication rvec -= this*mvec
    void MultSub( const CFSVector & mvec, CFSVector & rvec ) const {;};

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
    void DyadicMult( const CFSVector & vec1 );  
  
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
    void DyadicMult( const CFSVector & vec1, const CFSVector & vec2 ); 

    //! Calculate the Determinant (up to size 3)

    //! Calculates the determinant for a square matrix up to size 3x3.
    //! For larger matrices an error is returned.
    //! \param val (output) Return value of the method
    void Determinant( TYPE & val ) const;

    //@}

    //@{
    //! Invert the matrix and store it in 'inv' (up to size 3)
    void Invert( CFSMatrix & inv ) const {
      Error( "!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    }
    void Invert ( Matrix <TYPE> & inv ) const;
    //@}
    
    //@{
    //! Transpose the matrix and store the result in \a transposedMat
    //! \note The matrix itself gets not changed.
    //! \note If the transposed of a matrix is needed for a operation
    //! with a vector, the according function like 'MultT' should be used
    void Transpose( Matrix<TYPE> & transposedMat ) const;  
    void Transpose( CFSMatrix & transposedMat ) const {
      Error("!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    }
    //@}

    
#ifdef EXPR_TEMPLATES
    // =======================================================================
    // INTERFACE TO EXPRESSION TEMPLATES
    // =======================================================================
        
    //@{ 
    //! \name Interface To Expression Template Headers

    //! Matrix assignment operator using expression templates
    inline Matrix<TYPE>& operator=( const Matrix<TYPE>& rhs ) { 
      return assignFrom( rhs ); 
    }
    
    //! Scalar assignment operator using expression templates
    inline Matrix<TYPE>& operator=( TYPE rhs ) { 
      return assignFrom( rhs ); 
    }
    
    //! Matrix-Expression assignment operator using expression templates
    template <class X> inline Matrix<TYPE>& 
    operator=( const Xpr2<TYPE,X>& rhs ) {
      return assignFrom( rhs );
    }
    
    //! Abstract matrix assignment operator
    template <class M> inline Matrix<TYPE>& 
    operator=( const Dim2<TYPE,M>& rhs ) {
      return assignFrom(rhs);
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

    //! Create new matrix by Matrix- matrix multiplication (type promotion)
    template <class TYPE2>
    Matrix<PROMOTE(TYPE,TYPE2)> operator*( const Matrix<TYPE2> &y ) const;

    //! Multiply matrix by a scalar (this *= y)
    Matrix<TYPE> & operator*=( const TYPE &y );

    //! Perform matrix-matrix multiplication (this = this * arg)
    Matrix<TYPE> & operator*=( const Matrix<TYPE> &y );

    //! Divide matrix  by a scalar value
    Matrix<TYPE>  & operator/=( const TYPE &y );

    //@}

#endif // EXPR_TEMPLATES
    
    // =======================================================================
    // BOOLEAN OPERATORS
    // =======================================================================

    //! \name bool operators

    //@{
    
    //! Returns true if \a mat has the same entries as own matrix
    bool operator ==( const Matrix<TYPE> & mat ) const;

    //! Returns true if \a mat has different entries than own matrix
    bool operator!=( const Matrix<TYPE> & mat ) const;
 
    //@}

#ifdef USE_LAPACK
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
    //! Compile with LAPACK - Support (USE_LAPACK = yes)
    void solveWithLapack( Matrix<Complex> & b1,
                          lapackSysMatType & LAPACK_MATRIX_TYPE );
    
    //! Computes eigenvalues of an hermitian matrix
    void eigenvaluesWithLapack(Vector<Double> & b1);

    //! Converts a fortran 77 complex to a C++ complex
    void F772CC( const F77complex16 &v, std::complex<double> &val ) {
      std::complex<double> aux(v.real,v.imag);
      val = aux;
    }

    //! Converts a C++ complex to a fortan 77 one
    void CC2F77( const std::complex<double> &v, F77complex16 &val ) {
      val.real = (F77real8)v.real();
      val.imag = (F77real8)v.imag();
    }
    
    //! Converts a fortran 77 double to a C++ double
    void F772CC( const F77real8 &v, double &val ) {
      val = (double)v;
    }
    //! Converts a C++ double to a fortran 77 double
    void CC2F77( const double &v, F77real8 &val ) {
      val = (F77real8)v;
    }
    //@}
#endif
  
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
    //! \note The Matrix A=LU contains afterwards the the values of L 
    //! in the lower triangular, and the values of U in the upper part.
    void DirectSolve( CFSVector & x, CFSVector & b );


    //! scales the diagonal elements of a  matrix by a factor
    void ScaleDiagElems( const TYPE factor );
    
    //! Add a row to Matrix at position i
    void AddRow( const Vector<TYPE> & x, const UInt pos );
    
    //! Add a column to Matrix at position i
    void AddColumn( const Vector<TYPE> & x, const UInt pos ); 
    
    //! Return a special part ( real, imag, amplitude, phase) of a matrix
    Matrix<Double> GetPart(  DataType part ) const;

    //! Set special part ( real, imag, amplitude, phase) of a matrix
    void SetPart( DataType part, const Matrix<Double> & partMatrix );

    //! Return a sub-part of the own matrix

    //! Copies a sub-matrix at the position (row, col) into subMat. 
    //! The amount of copied elements depends on the size of subMat.
    void GetSubMatrix( CFSMatrix &subMat, const UInt nRows,
                       const UInt nCols ) const {
      Error( "!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    };

    //! Return a sub-part of the own matrix
    
    //! Copies a sub-matrix at the position (row, col) into subMat, 
    //! the amount of copied elements depends on the size of subMat
    void GetSubMatrix( Matrix<TYPE>& subMat, UInt row, UInt col ) const;
  
    //! Set a sub-part of the matrix
    
    //! Overwrites the matrix elements at the position (row, col) with subMat
    //! in a rectangular (submatrix) way.
    void SetSubMatrix( const CFSMatrix & subMat, const UInt nRows,
                       const UInt nCols ) {
      Error("!!! IMPLEMENT !!!", __FILE__, __LINE__ );
    };


    //! Set a sub-part of the matrix
    
    //! Overwrites the matrix elements at the position (row, col) with subMat
    //! in a rectangular (submatrix) way
    void SetSubMatrix( const Matrix<TYPE>& subMat, UInt row, UInt col );

    //! Adds a subMat to the matrix elements at the position (row, col)
    //! in a rectangular (submatrix) way
    void AddSubMatrix( const Matrix<TYPE>& subMat, UInt row, UInt col );

    //! Converts a matrix into a vector, by appending successively all rows
    void ConvertToVec_AppendRows( CFSVector& vec ) const;

    //! Converts a matrix into a vector, by appending successively all cols
    void ConvertToVec_AppendCols( CFSVector& vec ) const;

    //@}

  private:

    //! Calculates the adjunct of the matrix at position (i,j)
    TYPE Adjunct (UInt i, UInt j) const;

    //! Number of rows 
    UInt size_row_;
  
    //! Number of columns
    UInt size_col_;

    //! Data of the matrix
    TYPE** data_;

    // =======================================================================
    // SERIALIZATION FUNCTIONS
    // =======================================================================
    // These functions allow us to write a vector directly
    // into an boost::archive, for saving on a disk or in a 
    // iostream object

    //! allow serialization class to access vector entries
    friend class boost::serialization::access;
    
    //! Saving internal state into a boost::archive
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;
    
    //! Reading internal state from a boost::archive
    template<class Archive>
    void load(Archive & ar, const unsigned int version);
    
    //! The following statement is needed for boost
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    

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



  // =======================================================================
  // INLINE MEMBER IMPLEMENTATION
  // =======================================================================

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
  inline TYPE *  Matrix<TYPE>::operator[] (const UInt i) const
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

  // =======================================================================
  //  Inline part for all operators using type promotion
  //  rr being defined only in non-template-expression case
  // =======================================================================

#ifndef EXPR_TEMPLATES

  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::
  operator+(const Matrix<TYPE2> &x) const
  {
    ENTER_IFCN("Matrix::operator+");
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif 
    
#ifdef CHECK_INDEX
    if (size_row_ != x.GetSizeRow() || size_col_ != x.GetSizeCol())
      Error("incompatible dimension",__FILE__,__LINE__);
#endif
  
    Matrix<PROMOTE(TYPE,TYPE2)> z(size_row_,size_col_);
  
    UInt k;
    for ( k = 0; k < size_row_*size_col_; k++)
      z [0][k] = x[0][k]+data_[0][k];
  
    return z;
  }

  
  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::
  operator-(const Matrix<TYPE2> &x) const
  {
    ENTER_IFCN("Matrix::operator-");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.GetSizeRow() == 0 || x.GetSizeCol() == 0)
      Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX  
    if (size_row_ != x.GetSizeRow() || size_col_ != x.GetSizeCol())
      Error("incompatible dimension for +",__FILE__,__LINE__); 
#endif
  
    Matrix<PROMOTE(TYPE,TYPE2)> z(size_row_,size_col_);
  
    UInt k;
    for ( k = 0; k < size_row_*size_col_; k++)
      z[0][k] = -x[0][k]+data_[0][k];
  
    return z;
  }

  template<class TYPE> template<class TYPE2>
  Matrix<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::
  operator* (const TYPE2 &x) const 
  { 
    ENTER_IFCN("Matrix::operator*");
  
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
    UInt k;
  
    Matrix<PROMOTE(TYPE,TYPE2)> z(size_row_,size_col_);
  
    for ( k = 0; k < size_row_*size_col_; k++)
      z [0][k] = data_[0][k]*x;
  
    return z;
  }

  template<class TYPE>  template<class TYPE2>
  Vector<PROMOTE(TYPE,TYPE2)> Matrix<TYPE>::
  operator*(const Vector<TYPE2> &x) const
  {
    ENTER_IFCN("Matrix::operator*");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
    if (x.GetSize() == 0) Error("undefined Vector",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (size_col_ != x.GetSize()) Error("incompatible dimension",
                                        __FILE__,__LINE__);
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
    ENTER_IFCN("Matrix::operator*");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.GetSizeRow() == 0 || x.GetSizeCol()== 0)
      Error("undefined Matrix",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX  
    if (size_col_ != x.GetSizeRow())
      Error("incompatible dimension",__FILE__,__LINE__);
#endif
 
    PROMOTE(TYPE,TYPE2) a;
    Matrix<PROMOTE(TYPE,TYPE2)>  z (size_row_, x.GetSizeCol());
  
    UInt i,j; 
    for (i = 0; i < size_row_; i++)
      for (j = 0; j < x.GetSizeCol(); j++)
        {       
          a = data_ [i] [0] * x[0][j];
          for (UInt k = 1; k < size_col_; k++)
            a += data_ [i] [k] * x[k][j];
          z(i,j) = a;
        }
  
    return z;
  }
#endif //EXPR_TEMPLATES

  template<class TYPE> template< class Archive>
  void Matrix<TYPE>::save(Archive & ar, const unsigned int version) const {
    
    // invoke serialization of the base class 
    ar & boost::serialization::base_object<CFSMatrix>(*this);
    
    // save own members
    ar & size_row_;
    ar & size_col_;
    
    for( UInt i = 0 ; i < size_row_; i++ ) {
      for( UInt j = 0; j < size_col_; j++ ) {
        ar & data_[i][j];
      }
    }
  }
  
  template<class TYPE> template <class Archive>
  void Matrix<TYPE>::load(Archive & ar, const unsigned int version) {

    // invoke serialization of the base class 
    ar & boost::serialization::base_object<CFSMatrix>(*this);

    // check if data is already present
    if (data_ != NULL)
      {
        delete[] data_[0];
        delete[] data_;
      }

    // invoke serialization of the base class 
    ar & boost::serialization::base_object<CFSMatrix>(*this);

    ar & size_row_;
    ar & size_col_;

    // create storage for data to read in
    data_ = new TYPE* [size_row_];
    data_[0]=new TYPE[size_col_*size_row_];
    for (UInt k=1; k < size_row_; k++) 
      data_[k]=data_[k-1]+size_col_;

    // copy data itself from archive
    for( UInt i = 0 ; i < size_row_; i++ ) {
      for( UInt j = 0; j < size_col_; j++ ) {
        ar & data_[i][j];
      }
    }

  }
} //end of namespace


#endif  // FILE_MATRIX

