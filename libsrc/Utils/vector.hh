#ifndef FILE_VECTOR_2001
#define FILE_VECTOR_2001
 
#include <iostream>

#include "tools.hh"
//#include "matrix_head.hh"


namespace CoupledField
{          

  template<class TYPE> class Vector;

  //! Function for swap  number a and b (Integer, Double)
  template <class S> void swap(S &, S &);

  //! Function for swap Vector<T>;
  template<class T>
  void swap(Vector<T> & a, Vector<T> & b);

  //! Sort of vector with size n
  template <class S> void sort(S* v, Integer n);

  template<class TYPE> class SparseMatrix;
  template<class TYPE> class SymSparseMatrix;
  template<class TYPE> class SymBandMatrix;
  template<class TYPE> class BandMatrix;
  template<class TYPE> class SymMatrix;
  template<class TYPE> class Matrix;

  //    template<class Type> istream &operator>> (istream &,Vector<TYPE> &);
  //! Overloading << for class vector
  template<class TYPE>  std::ostream& operator << ( std::ostream & , const Vector<TYPE> &);

//! class Vector for working with arrays of Integer or Double numbers 
template<class TYPE> class Vector
{
  //! size
  Integer	n;

  //! pointer to array
  TYPE	*p;

  void	help	(const Integer, std::istream &);
        
public:

  static Vector<TYPE> null;

  friend class SymSparseMatrix<TYPE>;
  friend class SparseMatrix<TYPE>;
  friend class BandMatrix<TYPE>;
  friend class SymBandMatrix<TYPE>;
  friend class SymMatrix<TYPE>;
  friend class Matrix<TYPE>;

  template<class T>
  friend void swap(Vector<T> &, Vector<T> &);

  //! Constructor
  Vector	();

  //! Constructor with parameter dimension of Vector
  Vector	(Integer adim);

  //! Constructor of Vector with initialization with Vector x
  Vector	(Integer, const TYPE *const x);

  //! Copy of Vector
  Vector	(const Vector &);

  //! Deconstructor
  ~Vector	();

  //! Change size of vector
  void Resize(const Integer i);

  //! Allocate vector of size i
  void Allocate(const Integer i);       

  //! Initialize vector by zero
  void Init(const Integer l=0);

  //! Overloading of operation =
  Vector	&operator=	(const Vector &);

  //! Element can be referred to as v[i]
  TYPE	&operator[]	(const Integer) const;

  //! Return dimension of Vector
  Integer	size () const;

  //! Return pointer p to array 
  TYPE*  get() const;

  //! Overloading of operations +,+=
  Vector	operator+	() const;
  Vector	operator+	(const Vector &) const;
  Vector	&operator+=	(const Vector &);

  //! Overloading of operations -,-=
  Vector	operator-	() const;
  Vector	operator-	(const Vector &) const;
  Vector	&operator-=	(const Vector &);

  //! Overloading of operations / for number and Vector
  Vector  operator/       (const TYPE &) const;

  //! Overloading of operations * for number and Vector   
  Vector	operator*	(const TYPE &) const;

  //! Overloading of operations * for Vector and Vector 
  TYPE	operator*	(const Vector &) const;

  //! Overloading of operations * for Vector and Matrix  
  Vector	operator*	(const Matrix<TYPE> &) const;

  //! Overloading of operation *= for Vector and number
  Vector	&operator*=	(const TYPE &);

  //! Overloading of operations * for Vector and Matrix 
  Vector	&operator*=	(const Matrix<TYPE> &);

  //! Overloading of operation equal for Vector
  Integer	operator==	(const Vector &) const;

  //! Overloading of operation not equal for Vector
  Integer	operator!=	(const Vector &) const;

  //! Return part of Vector from index i to ii
  Vector	part	(const Integer i, const Integer ii) const;

  //! 
  static Vector	unit	(const Integer n, const Integer i);

  //! Calculate norm L^2 for Vector
  Double norm_2 ();

  //! Add element of the same type at position pos, by default to the beginning Beware of numeration in C++
  void add(const TYPE & y, Integer pos=0);

  //! Add vector to vector at position pos
  void  add (const Vector<TYPE> & y, Integer pos=0);

  //! Delete element from vector on position pos
  void  cut (const Integer pos);

  //! Delete elements from position pos1 to pos2, on pos1, pos2 too
  void  cut (const Integer pos1, const Integer pos2);

  //! Return size of space memory of this vector
  Integer  memory() const;

  //!Friends
  //! Sort of vector: v - vec.p, n - vec.size
  template <class S> void sort(S* v, Integer n);
  //! Swap 2 elements in vector Ex swap(v[i],v[j])
  template<class T> void swap(T & a, T & b);

};

// inline member functios

template<class TYPE>
inline Vector<TYPE>::~Vector ()
{
	if (p)
		delete [] p;
}

template<class TYPE>
inline Integer Vector<TYPE>::size () const
{
	return n;
}

template <class TYPE>
inline TYPE* Vector<TYPE>::get ()  const
{	if (!p)
	Error("Vector: undefined Vector",__FILE__,__LINE__);
  return p;
}

#ifdef __GNUC__
template class Vector<Integer>;
template class Vector<Double>;
#endif
//template class Vector<Integer>;
//template class Vector<Point2D>;

}
#endif	// FILE_VECTOR


