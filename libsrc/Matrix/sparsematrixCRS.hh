#ifndef FILE_SPARSEMATRIXCRS_2001
#define FILE_SPARSEMATRIXCRS_2001

#include "matrix.hh"
#include "absmatrix.hh"
 
namespace CoupledField
{
 
template<class TYPE> class SparseMatrix;

template<class TYPE>  std::ostream& operator<<( std::ostream & , const SparseMatrix<TYPE> &); 
/// Calculate Spur of SparseMatrix
template <class T> T Spur(const SparseMatrix<T> &);

template<class TYPE> class SparseMatrix
{
  /// number of entries
  Integer numentry;
  /// number of rows of matrix
  Integer row;
  /// number of colomnes of matrix
  Integer col;
  /// pointer to array of nonzero entries
  TYPE * p;
  /// pointer to array of column indexes
  Integer * pc;
  /// pointer to array of first position nonzero entry at row n 
  Integer * pf;  

public:
  /// Constructor  
  SparseMatrix();
  /// Constructor(number of entries, number of rows, number of colomnes)
  SparseMatrix(const Integer, const Integer, const Integer);
  ///  Default Copy Constructor
  SparseMatrix(const SparseMatrix & x);
  /// Constructor(numentry, size)
  SparseMatrix(const Integer, const Integer);
  /// Constructor (Matrix)
  SparseMatrix(const Matrix<TYPE> &x);
  /// Change size of matrix (numentry,row,col)
  void Resize(const Integer, const Integer, const Integer); 
  /// Deconstructor
  ~SparseMatrix();
  /// Return number of entries
  Integer GetNumentry() const;
  /// Return number of rows
  Integer GetRow() const;
  /// Return number of colomnes
  Integer GetCol() const;
  /// Return p
  TYPE * get() const;
  /// Return pointer to row number []
  TYPE * operator[] (const Integer) const;
  /// Overloading for operations
  /// Copying
  SparseMatrix &operator=  (const SparseMatrix &);  
  /// Overloading of operations
  SparseMatrix  &     operator+       () ;
  ///
  SparseMatrix  & operator+= (const SparseMatrix<TYPE> &) ; 
  ///
  SparseMatrix operator+ (const SparseMatrix &) const;
  ///
  SparseMatrix       operator-       () const;
  ///
  SparseMatrix       operator-       (const SparseMatrix &) const;
  ///
  SparseMatrix       &operator-=     (const SparseMatrix &);
  ///
  SparseMatrix       operator*       (const TYPE &) const;
  ///
  Vector<TYPE>    operator*       (const Vector<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const Matrix<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const SparseMatrix<TYPE> &) const;
  ///
  SparseMatrix       &operator*=     (const TYPE &);
  ///
  Integer operator==      (const SparseMatrix &) const;
  ///
  Integer operator!=      (const SparseMatrix &) const;
  /// Cut part of matrix(index) /numeration from 1, not from 0 as in C++/
  SparseMatrix       part    (const Integer) const;
  /// void initialize matrix by zero
  void Init();
  /// Return FALSE
  Boolean IsSymmetric() const;
  /// precondition
  void precond(Vector<TYPE> & e, const Vector<TYPE> r, enum precond type);

  template<class S>
friend std::ostream & operator<< (std::ostream & out, const SparseMatrix<S> &mat);

  template<class S>
  friend S Spur(const SparseMatrix<S> &x);

};

template<class TYPE>
inline SparseMatrix<TYPE>::~SparseMatrix()
{
#ifdef TRACE
 (*trace) << "entering SparseMatrix::~SparseMatrix" << std::endl;
#endif

 if (p) delete [] p;
 if (pc) delete [] pc;
 if (pf) delete [] pf;
 
}

template<class TYPE>
inline Integer SparseMatrix<TYPE>::GetRow () const
{
  return row;
}

template<class TYPE>
inline Integer SparseMatrix<TYPE>::GetCol () const
{
  return col;
}

template<class TYPE>
inline Integer SparseMatrix<TYPE>::GetNumentry () const
{
  return numentry;
}

template<class TYPE>
inline TYPE * SparseMatrix<TYPE>::get () const
{
  return p;
}

template<class TYPE>
inline Boolean SparseMatrix<TYPE>::IsSymmetric () const
{
  return FALSE;
}

template class SparseMatrix<Double>;
//template class SparseMatrix<Integer>;

} // end of namespace
#endif //









