#ifndef FILE_MATRIX_2001
#define FILE_MATRIX_2001

#include "Utils/vector.hh"
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

  //! Constructor row=0,col=0,p=NULL
  Matrix	();

  //! Constructor(row, col)
  Matrix	(const Integer, const Integer);

  //! Constructor (number of vectors, array of vectors(colomn))
  Matrix	(const Integer, const Vector<TYPE> * const);

  //! Default Copy Construcctor
  Matrix	(const Matrix &);

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

  //!
  Matrix	operator*	(const TYPE &) const;

  //!
  Vector<TYPE>	operator*	(const Vector<TYPE> &) const;

  //!
  Matrix	operator*	(const Matrix &) const;

  //!
  Matrix	&operator*=	(const TYPE &);

  //!
  Matrix	&operator*=	(const Matrix &);

  //!
  Integer	operator ==	(const Matrix &) const;

  //!
  Integer	operator!=	(const Matrix &) const;

  //! Cut part of matrix (left index row, right, upper index col, low )
  Matrix	part	(const Integer, const Integer,
                         const Integer, const Integer) const;

  //! Cut row number i, colomn number j from matrix
  void cut(const Integer i, const Integer j);

  //! Add a row to Matrix at position i
  void add_row(const Vector<TYPE> & x, const Integer pos );

  //! Add a colomn to Matrix at position i
  void add_col(const Vector<TYPE> & x, const Integer pos ); 

  //! Precondition: Implemented Jacobi, SSOR(omega=1.2), LU
  void precond(Vector<TYPE> &, const Vector<TYPE>, enum precond type);

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
{       if (!row || !col) Error("undefined Matrix");
        return row;
}
 
template<class TYPE>
inline Integer Matrix<TYPE>::size_col () const
{       if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
        return col;
}

#ifdef __GNUC__
template class Matrix<Double>;
template class Matrix<Integer>;
#endif

} //end of namespace
#endif	// FILE_MATRIX

