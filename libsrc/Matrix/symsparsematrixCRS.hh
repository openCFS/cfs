#ifndef FILE_SYMSPARSEMATRIXCRS_2001
#define FILE_SYMSPARSEMATRIXCRS_2001

#include "symmatrix.hh"
#include "absmatrix.hh"
 
namespace CoupledField
{
 
template<class TYPE> class SymSparseMatrix;

template<class TYPE>  std::ostream& operator<<( std::ostream & , const SymSparseMatrix<TYPE> &); 
/// Calculate Spur of SymSparseMatrix
template <class T> T Spur(const SymSparseMatrix<T> &);

template<class TYPE> class SymSparseMatrix
{
public:
  /// number of entries
  Integer numentry;
  /// size of matrix
  Integer row;
  /// pointer to array of nonzero entries
  TYPE * p;
  /// pointer to array of column indexes
  Integer * pc;
  /// pointer to array of first position nonzero entry at row n 
  Integer * pf;  

public:
  /// Constructor  
  SymSparseMatrix();
  /// Constructor(number of entries, number of rows)
  SymSparseMatrix(const Integer, const Integer);
  ///  Default Copy Constructor
  SymSparseMatrix(const SymSparseMatrix & x);
  /// Constructor (Matrix)
  SymSparseMatrix(const SymMatrix<TYPE> &x);
  /// Change size of matrix (numentry,size)
  void Resize(const Integer, const Integer); 
  /// Deconstructor
  ~SymSparseMatrix();
  /// Return number of entries
  Integer GetNumentry() const;
  /// Return number of rows
  Integer GetRow() const;
  /// Return p
  TYPE * get() const;
  /// Return pointer to row number []
  TYPE * operator[] (const Integer) const;
  /// Overloading for operations
  /// Copying
  SymSparseMatrix &operator=  (const SymSparseMatrix &);  
  /// Overloading of operations
  SymSparseMatrix  &     operator+       () ;
  ///
  SymSparseMatrix  & operator+= (const SymSparseMatrix<TYPE> &) ; 
  ///
  SymSparseMatrix operator+ (const SymSparseMatrix &) const;
  ///
  SymSparseMatrix       operator-       () const;
  ///
  SymSparseMatrix       operator-       (const SymSparseMatrix &) const;
  ///
  SymSparseMatrix       &operator-=     (const SymSparseMatrix &);
  ///
  SymSparseMatrix       operator*       (const TYPE &) const;
  ///
  Vector<TYPE>    operator*       (const Vector<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const Matrix<TYPE> &) const;
  ///
  SymSparseMatrix       &operator*=     (const TYPE &);
  ///
  Integer operator==      (const SymSparseMatrix &) const;
  ///
  Integer operator!=      (const SymSparseMatrix &) const;
  /// Cut part of matrix(index) /numeration from 1, not from 0 as in C++/
  SymSparseMatrix       part    (const Integer) const;
  /// void initialize matrix by zero
  void Init();
  /// Return FALSE
  Boolean IsSymmetric() const;

  template<class S>
friend std::ostream & operator<< (std::ostream & out, const SymSparseMatrix<S> &mat);

  template<class S>
  friend S Spur(const SymSparseMatrix<S> &x);

};

template<class TYPE>
inline SymSparseMatrix<TYPE>::~SymSparseMatrix()
{
#ifdef TRACE
 (*trace) << "entering SymSparseMatrix::~SymSparseMatrix" << std::endl;
#endif

 if (p) delete [] p;
 if (pc) delete [] pc;
 if (pf) delete [] pf;
 
}

template<class TYPE>
inline Integer SymSparseMatrix<TYPE>::GetRow () const
{
  return row;
}

template<class TYPE>
inline Integer SymSparseMatrix<TYPE>::GetNumentry () const
{
  return numentry;
}

template<class TYPE>
inline TYPE * SymSparseMatrix<TYPE>::get () const
{
  return p;
}

template<class TYPE>
inline Boolean SymSparseMatrix<TYPE>::IsSymmetric () const
{
  return FALSE;
}

template class SymSparseMatrix<Double>;
template class SymSparseMatrix<Integer>;

} // end of namespace
#endif //









