#ifndef FILE_BANDMATRIX_2001
#define FILE_BANDMATRIX_2001
 
namespace CoupledField
{
 
 template<class TYPE> class BandMatrix;

template<class TYPE>  std::ostream& operator<<( std::ostream & , const BandMatrix<TYPE> &); /// Calculate Spur of BandMatrix
template <class T> T Spur(const BandMatrix<T> &);

template<class TYPE> class BandMatrix
{
  /// number of rows
  Integer row;
  /// number of elements on the right of diagonal - right index
  Integer ri;
  ///  number of elements on the left of diagonal - left index
  Integer li;
  /// pointer to array of pointers to entries
  TYPE ** p;
  /// number of entries in matrix
  Integer numentry;

public:
  template<class S>
  friend S Spur(const BandMatrix<S> &);
  /// Constructor  
  BandMatrix();
  ///  Default Copy Constructor
  BandMatrix(const BandMatrix & x);
  /// Constructor(row, right col, left col)
  BandMatrix(Integer,Integer,Integer);
  /// Change size of matrix (row,left index, right index)
  void Resize(const Integer, const Integer, const Integer); 
  /// Deconstructor
  ~BandMatrix();
  /// Return number of rows
  Integer GetRow() const;
  /// Return left index
  Integer GetLi() const;
  /// Return right index
  Integer GetRi() const;
   /// Return number of non-zero entries in matrix 
  Integer GetNumentry() const;
  /// Return p
  TYPE ** get() const;
  /// Return pointer to row number []
  TYPE * operator[] (const Integer) const;
  /// Overloading for operations
  /// Copying
  BandMatrix &operator=  (const BandMatrix &);  
  /// Overloading of operations
  BandMatrix  &     operator+       () ;
  ///
  BandMatrix  & operator+= (const BandMatrix<TYPE> &) ; 
  ///
  BandMatrix operator+ (const BandMatrix &) const;
  ///
  BandMatrix       operator-       () const;
  ///
  BandMatrix       operator-       (const BandMatrix &) const;
  ///
  BandMatrix       &operator-=     (const BandMatrix &);
  ///
  BandMatrix       operator*       (const TYPE &) const;
  ///
  Vector<TYPE>    operator*       (const Vector<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const Matrix<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const BandMatrix<TYPE> &) const;
  ///
  BandMatrix       &operator*=     (const TYPE &);
  ///
  Integer operator==      (const BandMatrix &) const;
  ///
  Integer operator!=      (const BandMatrix &) const;
  /// Cut part of matrix(index)
  BandMatrix       part    (const Integer) const;
  /// Return FALSE
  Boolean IsSymmetric() const;
};

template<class TYPE>
inline BandMatrix<TYPE>::~BandMatrix()
{
#ifdef TRACE
 (*trace) << "entering BandMatrix::~BandMatrix" << endl;
#endif
 if (p) {  
          delete [] p[0];
          delete [] p; }
 
}

template<class TYPE>
inline Integer BandMatrix<TYPE>::GetRow () const
{
  return row;
}

template<class TYPE>
inline Integer BandMatrix<TYPE>::GetRi () const
{
  return ri;
}

template<class TYPE>
inline Integer BandMatrix<TYPE>::GetLi () const
{
  return li;
}

template<class TYPE>
inline Integer BandMatrix<TYPE>::GetNumentry () const
{
  return numentry;
}

template<class TYPE>
inline TYPE ** BandMatrix<TYPE>::get () const
{
  return p;
}

template<class TYPE>
inline Boolean BandMatrix<TYPE>::IsSymmetric () const
{
  return FALSE;
}

#ifdef __GNUC__
template class BandMatrix<Double>;
template class BandMatrix<Integer>;
#endif

} // end of namespace
#endif //
