#ifndef FILE_SYMMATRIX_2001
#define FILE_SYMMATRIX_2001

#include "matrix.hh"
#include "absmatrix.hh"
 
namespace CoupledField
{      

 template<class TYPE> class SymMatrix;
/// Overloading << for class SymMatrix
template<class TYPE>  std::ostream& operator << ( std::ostream & , const SymMatrix<TYPE> &);
 /// Calculate Spur of SymMatrix
template <class T> T Spur(const SymMatrix<T> &);

/// Own class for working with symmetric atrixes
template<class TYPE> class SymMatrix: public AbsMatrix<SymMatrix<TYPE>, TYPE >
{
        /// size
	Integer	size;

        /// pointer to array of pointers to row
	TYPE	** p;

public:
   	friend class Vector<TYPE>;
        friend class SymSparseMatrix<TYPE>;

        /// Constructor size=0,p=NULL
	SymMatrix	();
        /// Constructor(size)
	SymMatrix	(const Integer);
        /// Default Copy Construcctor
	SymMatrix	(const SymMatrix &);
        /// Change size of matrix (size)
        void Resize(const Integer, const Integer);
        /// Deconstructor
	~SymMatrix	();
        /// Add value to matrix on position i,j
        void Add(const Integer, const Integer, const TYPE);
        /// Copying
	SymMatrix	&operator=	(const SymMatrix &);
        /// Return pointer to row number []
	TYPE * operator[] (const Integer) const;
        /// Return size
	Integer	getSize() const;
        /// return p
        TYPE ** get() const;
        /// cut row and col with number ai and aj 
        void cut(const Integer ai, const Integer aj);
        /// Overloading of operations
        ///
       TYPE & operator()(const Integer i, const Integer j);
        ///
	SymMatrix	operator+	() const;
        ///
	SymMatrix	operator+	(const SymMatrix &) const;
        ///
	SymMatrix	&operator+=	(const SymMatrix &);
        ///
	SymMatrix	operator-	() const;
        ///
	SymMatrix	operator-	(const SymMatrix &) const;
        ///
	SymMatrix	&operator-=	(const SymMatrix &);

        ///
	SymMatrix	operator*	(const TYPE &) const;
        ///
	Vector<TYPE>	operator*	(const Vector<TYPE> &) const;
        ///
	Matrix<TYPE>    operator* (const Matrix<TYPE> &) const;
        ///
        Matrix<TYPE>    operator* (const SymMatrix<TYPE> &) const;
        ///
	SymMatrix	&operator*=	(const TYPE &);
        ///

	Integer	operator==	(const SymMatrix &) const;
        ///
	Integer	operator!=	(const SymMatrix &) const;

        /// Cut part of matrix (index row, up to which matrix is cut )
	SymMatrix	part	(const Integer) const;

         /// Initialize matrix by zero
        void Init();

        /// Return TRUE
        Boolean IsSymmetric() const;

        ///  Precondition: Implemented Jacobi, SSOR(omega=1.2), LU
          void precond(Vector<TYPE> &, const Vector<TYPE>, const Integer type);

        //
        void CholerskyDecomposition();
};

// inline member functions

template<class TYPE> 
inline SymMatrix<TYPE>::~SymMatrix ()
{	
  if (p) { delete  p[0]; 
         delete  [] p;}
}

template<class TYPE>
inline TYPE & SymMatrix<TYPE>::operator ()(const Integer i, const Integer j)
{
  return  p[i][j];
}

template<class TYPE>
inline void SymMatrix<TYPE>::Add(const Integer i, const Integer j, const TYPE value)
{
  p[i][j]+=value;
}
  
template<class TYPE>
inline TYPE ** SymMatrix<TYPE>::get () const
{
  return p;
}

template<class TYPE>
inline void SymMatrix<TYPE>::Init ()
{
Integer i;
for (i=0; i< size ; i++) p[0][i]=0;
}

template<class TYPE>
inline Boolean SymMatrix<TYPE>::IsSymmetric () const
{
  return TRUE;
}

template<class TYPE>
inline Integer SymMatrix<TYPE>::getSize () const
{
  return size;
}

#ifdef __GNUC__
template class SymMatrix<Double>;
#endif

} //end of namespace
#endif	// FILE_MATRIX

