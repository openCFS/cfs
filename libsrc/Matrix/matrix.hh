#ifndef FILE_MATRIX_2001
#define FILE_MATRIX_2001

#include "Utils/vector.hh"
#include "Utils/array.hh"
#include "absmatrix.hh"
 
namespace CoupledField
{      

  //! Overloading << for class Matrix
  template<class TYPE>  std::ostream& operator << ( std::ostream & , const Matrix<TYPE> &);
  
  //! Calculate Spur of Matrix
  template <class T> T Spur(const Matrix<T> &);

//   //! Transporting Matrix 
//   template <class T> Matrix<T> Trans(const Matrix<T> &);


//! Own class for working with full matrixes
template<class TYPE>
class Matrix: public AbsMatrix<Matrix<TYPE> , TYPE>
{
  //! number of row :w
  Integer	row;

  //! number of col
  Integer col;

  //! pointer to array of pointers to row
  TYPE	** p;

public:
       
  friend class SparseMatrix<TYPE>;
  friend class SymSparseMatrix<TYPE>;
  friend class BandMatrix<TYPE>;
  friend class SymBandMatrix<TYPE>;
  friend class SymMatrix<TYPE>;
  friend class Vector<TYPE>;
  friend class Array<TYPE>;
  

  //! Constructor row=0,col=0,p=NULL
  Matrix	();

  //! Constructor(row, col)
  Matrix	(const Integer, const Integer);

  //! Constructor (number of vectors, array of vectors(colomn))
  Matrix	(const Integer, const Vector<TYPE> * const);

  //! Default Copy Construcctor
  Matrix	(const Matrix &);

  //! Change size of quadratic matrix
  void Resize(const Integer);

  //! Change size of matrix (row, col)
  void Resize(const Integer, const Integer);


  //! Initialized by zero matrix 
  void Init();

  //! Deconstructor
  ~Matrix	();

  //! Copying
  Matrix	&operator=	(const Matrix &);

  //! Return pointer to row number []
  TYPE * operator[] (const Integer) const;

  //! Return size of matrix , only for quadratic matrixes
  Integer getSize() const;

  //! Return number of rows
  Integer	size_row() const;

  //! Return number of colomnes
  Integer	size_col() const;

  /// calculates the determinante of matrices up to size 3
  Double Det () const;      


  /// fast inversion for matrices smaller than size 3
  void Invert (Matrix <Double> & inv) const;

  //! return p
  TYPE ** get() const;

  //! return pointer to array of elements in matrix, row by row
  TYPE * getinarray() const { return p[0];}

  //! Add element on position i,j
  void Add(const Integer, const Integer, const TYPE value);

  //! Overloading of operations
  //!
  TYPE & operator()(const Integer, const Integer);

  //! 
  Matrix	operator+	() const;

  //!
  Matrix	operator+	(const Matrix &) const;

  //!
  Matrix	&operator+=	(const Matrix &);

  //!
  Matrix	operator-	() const;

  //!
  Matrix	operator-	(const Matrix &) const;

  //!
  Matrix	&operator-=	(const Matrix &);

  //! multiplication with scalar value
  Matrix	operator*	(const TYPE &) const;

  //!
  Vector<TYPE>	operator*	(const Vector<TYPE> &) const;

  //!
  std::vector<TYPE> operator*(const std::vector<TYPE> &x) const;
  

  //! 
  Array<TYPE>   operator*       (const Array<TYPE> &) const;

  //!
  Matrix	operator*	(const Matrix &) const;

  //!
  Matrix	&operator*=	(const TYPE &);

  //!
  Matrix	&operator*=	(const Matrix &);

  //!
  Matrix	&operator/=	(const TYPE &);
  //!
  Integer	operator ==	(const Matrix &) const;

  //!
  Integer	operator!=	(const Matrix &) const;
 

  //! Cut part of matrix (left index row, right, upper index col, low )
//   Matrix	part	(const Integer, const Integer,
//                          const Integer, const Integer) const;

  //! Cut row number i, colomn number j from matrix
  void cut(const Integer i, const Integer j);

  //! Add a row to Matrix at position i
  void add_row(const Vector<TYPE> & x, const Integer pos );

  //! Add a colomn to Matrix at position i
  void add_col(const std::vector<TYPE> & x, const Integer pos ); 

  //! get a column out of the matrix
  std::vector<TYPE> get_col(const Integer acol);

  //! Precondition: Implemented Jacobi, SSOR(omega=1.2), LU
  void precond(Vector<TYPE> &, const Vector<TYPE>, enum precond type);

  /// Transpose actual matrix
  void Transpose (Matrix<TYPE> &transposedMat);  

  /// Dyadic multiplication of two vectors
  void DyadicMult(std::vector<TYPE> v1);

  /// Dyadic multiplication of two vectors
  void DyadicMult(std::vector<TYPE> v1, std::vector<TYPE> v2);

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

  //
  void CholerskyDecomposition();
  
  void ConvertToVec_RowsFirst(std::vector<TYPE>& vec) const;

  /// scales the diagonal elements of a  matrix by a factor
  void ScaleDiagElems(TYPE factor);

  /// gets the diagonal elements of a  matrix in a one column matrix
  void GetDiagInMatrix(Matrix<TYPE>& columnMat);
  //Matrix<TYPE> GetDiagInMatrix();

private:
  /// calculates the adjunct of the matrix at position (i,j)
  Double Adjunct (Integer i, Integer j) const;
  

  //!
  //       TYPE At( const Integer i, const Integer j) { return p[i][j]; }
};

// inline member functions

template<class TYPE> 
inline Matrix<TYPE>::~Matrix ()
{	
  if (p) { delete [] p[0]; 
         delete  [] p;}
}

template<class TYPE>
inline void Matrix<TYPE>::Init ()
{
Integer i;
for (i=0; i<row*col; i++) p[0][i]=0;
}

template<class TYPE>
inline void Matrix<TYPE>::Add ( const Integer i, const Integer j, const TYPE value)
{
 p[i][j]+=value;
}

template<class TYPE>
inline TYPE & Matrix<TYPE>::operator()(const Integer i, const Integer j) 
{
  return p[i][j];
}

template<class TYPE>
inline TYPE ** Matrix<TYPE>::get () const
{
  return p;
}

template<class TYPE>
inline Integer Matrix<TYPE>::getSize() const
{
 if (row!=col) Error("Function .size() is valid only for quadratic matrix");
 return row;
}

template<class TYPE>
inline Integer Matrix<TYPE>::size_row () const
{       if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
        return row;
}
 
template<class TYPE>
inline Integer Matrix<TYPE>::size_col () const
{       if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
        return col;
}

template<class TYPE>
inline Double Matrix<TYPE>::Det () const
{       
  if (!row || !col) Error("Undefined Matrix!",__FILE__,__LINE__);
  if (row != col ) Error("No quadratic matrix!",__FILE__,__LINE__);
  switch (row)
    {
    case 1: return p[0][0];
      break;
    case 2: return  p[0][0]*p[1][1]-p[0][1]*p[1][0];
      break;
    case 3: return p[0][0]*p[1][1]*p[2][2] +
	      p[0][1]*p[1][2]*p[2][0] +
	      p[0][2]*p[1][0]*p[2][1] -
	      p[0][2]*p[1][1]*p[2][0] -
	      p[0][1]*p[1][0]*p[2][2] -
	      p[0][0]*p[1][2]*p[2][1];
      break;
    default: Error("Dimension larger than 3!",__FILE__,__LINE__);
    }
}


//   template<TYPE>
//   friend  std::vector<TYPE> operator* (std::vector<TYPE> & vec, const Matrix<TYPE> & mat);

// std::vector<Double> operator* (std::vector<Double> & vec, const Matrix<Double> & mat);


Double operator* (std::vector<Double> & vec1, Vector<Double> & vec2);

Double operator* (std::vector<Double> & vec1, std::vector<Double> & vec2);

Double L2Norm(std::vector<Double> & vec);

std::vector<Double> operator*= (std::vector<Double> & vec, Double val);

std::vector<Double> operator/= (std::vector<Double> & vec, Double val);

std::vector<Double> operator+ (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator- (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator+= (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator-= (std::vector<Double> & vec, std::vector<Double> & vec2);

std::vector<Double> operator* (Double val, std::vector<Double> & vec);

//std::vector<Double> operator= (std::vector<Double> & vec, Double val);




#ifdef __GNUC__
template class Matrix<Double>;
template class Matrix<Integer>;
#endif

} //end of namespace

//typedef double Double;



#endif	// FILE_MATRIX

