#ifndef FILE_SPARSEMATRIX_2002
#define FILE_SPARSEMATRIX_2002

#include "matrix.hh"
#include "absmatrix.hh"

namespace CoupledField
{

template<class TYPE>
struct ElemSparseMatrix
{
  //! number of column
  Integer column;

  //! value
  TYPE value;

  //! pointer to next element
  ElemSparseMatrix<TYPE> * next;
};

template<class TYPE>
class SparseMatrix
{
public:
  //! Constructor
  SparseMatrix();
  //! Constructor(number of column, number of rows)
  SparseMatrix(const Integer anumrows, const Integer anumcol);
  //! Constructor(const Matrix<TYPE> & x)
  SparseMatrix(const Matrix<TYPE> & x);
  /// Default Copy Constructor
  SparseMatrix(const SparseMatrix & x);
  /// Change size of Matrix (numrow, numcol)
  void Resize(const Integer, const Integer);
  /// Deconstructor
  ~SparseMatrix();
  /// Return number of rows
  Integer GetRow() const { return numrows;}
  /// Return number of colomns
  Integer GetCol() const {return numcol;}
  //// access to elements
  TYPE & operator()(Integer row,Integer col);

  //// overloading operators
  /// Copying
  SparseMatrix<TYPE> &operator=  (const SparseMatrix<TYPE> &);
  /// Overloading of operations
  SparseMatrix<TYPE>  &     operator+       () ;
  ///
  SparseMatrix<TYPE>  & operator+= (const SparseMatrix<TYPE> &) ;
  ///
  SparseMatrix<TYPE> operator+ (const SparseMatrix<TYPE> &) const;
  ///
  SparseMatrix<TYPE>       operator-       () const;
  ///
  SparseMatrix<TYPE>       operator-       (const SparseMatrix<TYPE> &) const;
  ///
  SparseMatrix<TYPE>       &operator-=     (const SparseMatrix<TYPE> &);
  ///
  SparseMatrix<TYPE>       operator*       (const TYPE &) const;
  ///
  Vector<TYPE>    operator*       (const Vector<TYPE> &) const;
  ///
  SparseMatrix<TYPE>       &operator*=     (const TYPE &);
  ///
  Integer operator==      (const SparseMatrix<TYPE> &) const;
  ///
  Integer operator!=      (const SparseMatrix<TYPE> &) const;

  template<class S>
  friend std::ostream & operator<< (std::ostream& , const SparseMatrix<S>&);

private:
  //! number of rows
  Integer numrows;
  //! number of column
  Integer numcol;
  //! a pointer to a pointer to Elements
  ElemSparseMatrix<TYPE> ** ptRow; 

  //! Auxialary function for copy matrix 
  void Copy(const SparseMatrix<TYPE> & x);
  //! Auxialary function for inserting element
  TYPE & InsertElement(ElemSparseMatrix<TYPE> * ptElem,const Integer row, const Integer col, const Integer first);
  //! Auxialary function for deleting element
  void DeleteElement(const Integer i, const Integer j);

public:
  void test();

};

//! Overloading << for class Matrix
template<class TYPE> 
std::ostream & operator << ( std::ostream & , const SparseMatrix<TYPE> &);

template<class TYPE>
inline SparseMatrix<TYPE>::~SparseMatrix()
{
  if (ptRow) delete [] ptRow;
}

#ifdef __GNUC__
template class SparseMatrix<Double>;
#endif

} // end of namespace

#endif
