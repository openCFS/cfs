#ifndef FILE_MATRIX_2004
#define FILE_MATRIX_2004

#include "cfsmatrix.hh"

namespace CoupledField
{      

//! Overloading << for class Matrix
template<class TYPE> class Vector;
template<class TYPE> class Array;

//! Concrete implementation of a dense matrix
template<class TYPE>
class Matrix: public CFSMatrix
{
public:

  // Friend declarations
  friend class Vector<TYPE>;
  friend class Array<TYPE>;
  
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
  Matrix (const Integer nRows, const Integer nCols);

  //! Constructor (number of vectors, array of vectors(colomn))
  Matrix (const Integer, const Vector<TYPE> * const);

  //! Default Copy Construcctor
  Matrix(const Matrix &);

  //! Destructor
  ~Matrix();
  
  //! Initialize matrix with a given entry.
  //! If no entry given, it gets initalized with zeroes
  /*!
    \param val (input,opt.) Entry the matrix gets initalized with
  */
  //! \note This method does not change the size of the matrix
  void Init(const TYPE val = 0.0);

  //! Change size of quadratic matrix
  /*!
    \param size (input) Number of rows / columns
  */
  //! \note the matrix contains afterwards only zeroes
  void Resize(const Integer size);
  
  //! Change size of general matrix 
  /*!
    \param nRows (input) Number of rows
    \param nCols (input) Number of columns)
  */
  //! \note the matrix contains afterwards only zeroes
  void Resize(const Integer nRows, const Integer nCols);

  //! Get the number of rows
  Integer GetSizeRow() const;
  
  //! Get the number of columns
  Integer GetSizeCol() const;
  
  //! Set the entry 'val' at position (row,col) in the matrix
  /*!
    \param row (input) Row of entry
    \param col (input) Column of entry
    \param val (input) Value to be set
  */
  void SetEntry(const Integer row, const Integer col, const TYPE val)
  {Error("!!! IMPLEMENT !!!");};
  
  //! Add'val' to the matrix entry at position (row,col) in the matrix
  /*!
    \param row (input) Row of entry
    \param col (input) Column of entry
    \param val (input) Value to be added
  */
  void AddToEntry(const Integer row, const Integer col, const TYPE val);
   
  //! Get the entry 'val' at position (row,col) in the matrix
  /*!
    \param row (input) Row of entry
    \param col (input) Column of entry
    \param val (output) Variable the value is written into
  */
  void GetEntry(const Integer row, const Integer col, TYPE & val) const
  {Error("!!! IMPLEMENT !!!");}; 
  
  //! Calculates the determinant (up to size 3)
  /*!
    \param val (output) Return value of the method
  */
   void Determinant(TYPE & val) const;

  //! Invert the matrix and store it in 'inv' (up to size 3)
  /*!
    \param inv (output) Inverse of the matrix
  */
  void Invert(CFSMatrix & inv) const
  {Error("!!! IMPLEMENT !!!");};
  
  //! Transpose the matrix and store it in 'trans'
  /*!
    \param trans (output) Transposed matrix
  */
  //! \note The matrix itself gets not changed.
  //! \note If the transposed of a matrix is needed for a operation
  //! with a vector, the according function like 'MultT' should be used
  void Transpose(CFSMatrix & trans) const
    {Error("!!! IMPLEMENT !!!");};
  
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
  //! \deprecated Do not use std::vectors for 
  //! numerical values anymore
  void DyadicMult(std::vector<TYPE> vec1);
  
  
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
  //! \deprecated Do not use std::vectors for 
  //! numerical values anymore
  void DyadicMult(std::vector<TYPE> vec1, std::vector<TYPE> vec2); 
  
  //! copies a submatrix at the position (row, col) into subMat, 
  //! the amount of copied elements depends on the size of subMat
  void GetSubMatrix(CFSMatrix & subMat, const Integer nRows, const Integer nCols) const 
  {Error("!!! IMPLEMENT !!!");};
  
  //! overwrites the matrix elements at the position (row, col) with subMat
  //! in a rectangular (submatrix) way
  void SetSubMatrix(CFSMatrix & subMat, const Integer nRows, const Integer nCols)
  {Error("!!! IMPLEMENT !!!");};
  
  //! scales the diagonal elements of a  matrix by a factor
  void ScaleDiagElems(const TYPE factor);

  //!
  void Add(const TYPE fac, const CFSMatrix & mat){};
  
  //! Perform a matrix-vector multiplication rvec = this*mvec
   void Mult(const CFSVector & mvec, CFSVector & rvec)const {;}

  //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
   void MultT(const CFSVector & mvec, CFSVector & rvec) const {;};

  //! Perform a matrix-vector multiplication rvec += this*mvec
   void MultAdd(const CFSVector & mvec, CFSVector & rvec)const {;};
  
  //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
   void MultTAdd(const CFSVector & mvec, CFSVector& rvec) const {;};
  
  //! Perform a matrix-vector multiplication rvec -= this*mvec
   void MultSub(const CFSVector & mvec, CFSVector & rvec) const {;};

  //////////////////////////////////////
  // Functions for working with other //
  // Matrix<Type> and Vector<TYPE>    //
  //////////////////////////////////////

  //! Assignment operator
  Matrix & operator= (const Matrix &);

  //! Return pointer to row number []
  inline TYPE * operator[] (const Integer) const;

  //! fast inversion for matrices smaller than size 3
  void Invert (Matrix <TYPE> & inv) const;

  //! returns pointer to continuos chunck of data
  TYPE ** GetRowPointer() const;

  //! return pointer to array of elements in matrix, row by row
  inline TYPE * GetDataPointer() const { return data_[0];}

  //! get the pointer to the continuous chunk of data as double *.
  //! If the matrix consists of Complex, it will be decomposed
  //! in a continous array of double (real_1, imag_1, real_2, ...)
  Double * GetDoublePointer();

  //! Overloading of operations

  //! Access operator
  /*!
    \param row (Input) Row number
    \param col (Input) Column number
  */
  TYPE & operator()(const Integer row , const Integer col);

  //! 
  Matrix operator+() const;

  //!
  Matrix operator+(const Matrix &) const;

  //!
  Matrix & operator+=(const Matrix &);

  //!
  Matrix operator-() const;

  //!
  Matrix operator-(const Matrix &) const;

  //!
  Matrix & operator-=(const Matrix &);

  //! multiplication with scalar value
  Matrix operator* (const TYPE &) const;

  //!
  Vector<TYPE> operator* (const Vector<TYPE> &) const;

  //!
  std::vector<TYPE> operator*(const std::vector<TYPE> &x) const;
  
  //! 
  //Array<TYPE> operator*(const Array<TYPE> &) const;

  //!
  Matrix operator*(const Matrix &) const;

  //!
  Matrix & operator*=(const TYPE &);


  //!
  Matrix & operator*=(const Matrix &);

  //!
  Matrix & operator/=(const TYPE &);
  //!
  Integer operator ==(const Matrix &) const;

  //!
  Integer operator!=(const Matrix &) const;
 
  //////////////////////////////////////////
  // Folgende Funktionen müssen grundlegend
  // auf 'Vector' statt 'std::vector' 
  // umgestellt werden
  ///////////////////////////////////////////

 //  //! Cut part of matrix (left index row, right, upper index col, low )
//   //   Matrix	part	(const Integer, const Integer,
//   //                          const Integer, const Integer) const;
  
//   //! Cut row number i, colomn number j from matrix
//   //void cut(const Integer i, const Integer j);

   //! Add a row to Matrix at position i
   void AddRow(const Vector<TYPE> & x, const Integer pos );

   //! Add a colomn to Matrix at position i
   void AddColumn(const std::vector<TYPE> & x, const Integer pos ); 

//   //! get a column out of the matrix
//   std::vector<TYPE> get_col(const Integer acol);

   /// Transpose actual matrix
  void Transpose (Matrix<TYPE> &transposedMat);  

  /// copies a submatrix at the position (row, col) into subMat, 
  /// the amount of copied elements depends on the size of subMat
  void GetSubMatrix(Matrix<TYPE>& subMat, Integer row, Integer col) const;

  /// overwrites the matrix elements at the position (row, col) with subMat
  /// in a rectangular (submatrix) way
  void SetSubMatrix(Matrix<TYPE>& subMat, Integer row, Integer col);

  /// converts a matrix into a vector, by appending successively all rows
  void ConvertToVec_AppendRows(std::vector<TYPE>& vec) const;

  /// converts a matrix into a vector, by appending successively all columns
  void ConvertToVec_AppendCols(std::vector<TYPE>& vec) const;

   /// gets the diagonal elements of a  matrix in a one column matrix
  void GetDiagInMatrix(Matrix<TYPE>& columnMat);
  //Matrix<TYPE> GetDiagInMatrix();

private:

  //! calculates the adjunct of the matrix at position (i,j)
  TYPE Adjunct (Integer i, Integer j) const;

  //! number of rows 
  Integer size_row_;
  
  //! number of columns
  Integer size_col_;

  //! data of the matrix
  TYPE** data_;


};

/////////////////////////////
// Inline member functions //
/////////////////////////////


template<class TYPE>
inline void Matrix<TYPE>::Init(const TYPE val)
{
  Integer i;
  for (i=0; i<size_row_*size_col_; i++) 
    data_[0][i]=val;
}

template<class TYPE>
inline void Matrix<TYPE>::AddToEntry ( const Integer i, const Integer j, const TYPE value)
{
  data_[i][j]+=value;
}

template<class TYPE>
inline TYPE & Matrix<TYPE>::operator()(const Integer row, const Integer col) 
{
  return data_[row][col];
}

template<class TYPE>
inline TYPE * Matrix<TYPE>::operator[] (const Integer i) const
{ 
  ENTER_IFCN("Matrix::operator[]");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) Error("undefined Matrix",__FILE__,__LINE__);
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
inline Integer Matrix<TYPE>::GetSizeRow() const
{       
  return size_row_;
}
 
template<class TYPE>
inline Integer Matrix<TYPE>::GetSizeCol () const
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

template<class TYPE>  std::ostream& operator << ( std::ostream & , const Matrix<TYPE> &);

Double operator* (std::vector<Double> & vec1, Vector<Double> & vec2);

Double operator* (std::vector<Double> & vec1, std::vector<Double> & vec2);

Double L2Norm(std::vector<Double> & vec);

//// EXPERIMENTELL
std::vector<Double> operator*= (std::vector<Double> & vec, Double val);

std::vector<Double> operator/= (std::vector<Double> & vec, Double val);

std::vector<Double> operator+ (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator- (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator+= (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator-= (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator* (Double val, std::vector<Double> & vec);


//#ifdef __GNUC__
template class Matrix<Double>;
template class Matrix<Integer>;
template class Matrix<Complex>;
//#endif

} //end of namespace


#endif	// FILE_MATRIX

