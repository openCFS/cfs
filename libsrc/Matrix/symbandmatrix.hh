#ifndef FILE_SYMBANDMATRIX_2001
#define FILE_SYMBANDMATRIX_2001
 
namespace CoupledField
{
 
template<class TYPE> class SymBandMatrix;

template<class TYPE>  ostream& operator<<( ostream & , const SymBandMatrix<TYPE> &); /// Calculate Spur of BandMatrix
template <class T> T Spur(const SymBandMatrix<T> &);

template<class TYPE> class SymBandMatrix:public AbsMatrix<SymBandMatrix<TYPE> , TYPE> 
{
public:
  /// number of rows
  /// suppose that matrix is quadratic
  Integer row;
  /// number of elements on the left of diagonal - left index
  Integer li;
  /// pointer to array of pointers to entries
  TYPE ** p;
  /// number of entries in matrix
  Integer numentry;

public:
  template<class S>
  friend S Spur(const SymBandMatrix<S> &);
  /// Constructor  
  SymBandMatrix();
  ///  Default Copy Constructor
  SymBandMatrix(const SymBandMatrix & x);
  /// Constructor(row, left col)
  SymBandMatrix(Integer,Integer);
  /// Change size of matrix (row,left index)
  void Resize(const Integer, const Integer); 
  /// Initialize by zero
  void Init();
  /// Deconstructor
  ~SymBandMatrix();
  /// Return number of rows
  Integer GetRow() const;
  /// Return left index
  Integer GetLi() const;
  /// Return number of entries in matrix 
  Integer GetNumentry() const;
  /// Cut row with number ai, col with number aj 
  void cut(const Integer ai, const Integer aj);
  /// Return p
  TYPE ** get() const;
  /// Return pointer to row number []
  TYPE * operator[] (const Integer) const;
  /// Overloading for operations
  /// Copying
  SymBandMatrix &operator=  (const SymBandMatrix &);  
  ///
  TYPE & At(const Integer, const Integer) ;
  /// Overloading of operations
  SymBandMatrix  &     operator+       () ;
  ///
  SymBandMatrix  & operator+= (const SymBandMatrix<TYPE> &) ; 
  ///
  SymBandMatrix operator+ (const SymBandMatrix &) const;
  ///
  SymBandMatrix       operator-       () const;
  ///
  SymBandMatrix       operator-       (const SymBandMatrix &) const;
  ///
  SymBandMatrix       &operator-=     (const SymBandMatrix &);
  ///
  SymBandMatrix       operator*       (const TYPE &) const;
  ///
  Vector<TYPE>    operator*       (const Vector<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const Matrix<TYPE> &) const;
  ///
  Matrix<TYPE> operator* (const SymBandMatrix<TYPE> &) const;
  ///
  SymBandMatrix       &operator*=     (const TYPE &);
  ///
  Integer operator==      (const SymBandMatrix &) const;
  ///
  Integer operator!=      (const SymBandMatrix &) const;
  /// Cut part of matrix(index)
  SymBandMatrix       part    (const Integer) const;
  ///  Precondition: Implemented Jacobi, SSOR(omega=1.2), LU
   void precond(Vector<TYPE> &, const Vector<TYPE>, const Integer type);
  /// Return TRUE 
  Boolean IsSymmetric() const;
  /// Return size (number of rows)
  Integer getSize(){ return row; }  
  template<class S>
  friend ostream & operator << (ostream & out, const SymBandMatrix<S> &);

};

template<class TYPE>
inline SymBandMatrix<TYPE>::~SymBandMatrix()
{
#ifdef TRACE
 cout << "entering SymBandMatrix::~SymBandMatrix" << endl;
#endif
 if (p) {  
          delete [] p[0];
          delete [] p; }
 
}

template<class TYPE>
inline Integer SymBandMatrix<TYPE>::GetRow () const
{
  return row;
}

template<class TYPE>
inline Integer SymBandMatrix<TYPE>::GetLi () const
{
  return li;
}

template<class TYPE>
inline void SymBandMatrix<TYPE>::Init() 
{
Integer i;
for (i=0; i<numentry; i++) p[0][i]=0;
}

template<class TYPE>
inline Integer SymBandMatrix<TYPE>::GetNumentry () const
{
  return numentry;
}

template<class TYPE>
inline TYPE ** SymBandMatrix<TYPE>::get () const
{
  return p;
}

template<class TYPE>
inline Boolean SymBandMatrix<TYPE>::IsSymmetric () const
{
  return TRUE;
}

template class SymBandMatrix<Double>;
//template class SymBandMatrix<Integer>;

} // end of namespace
#endif //
