#ifndef FILE_VECTOR_2004
#define FILE_VECTOR_2004

#include "cfsvector.hh"
#include "DataInOut/WriteInfo.hh"
namespace CoupledField {


// Forward class declarations
template<class TYPE> class Matrix;
template<class TYPE> class Vector;
template<class TYPE> class NodeStoreSol;
template<class TYPE> class ElemStoreSol;

//! Concrete Template class for a general dense vector
template<class TYPE>
class Vector : public CFSVector {
public:

  // Friend declarations
  friend class Matrix<TYPE>;
  friend class NodeStoreSol<TYPE>;
  friend class ElemStoreSol<TYPE>;
  template<class S>
  friend void Swap(Vector<S> &, Vector<S> &);

  //! Constructor
  Vector();

  //! Constructor with inital size.
  //! All entries are filled with zeroes
  Vector(const Integer size, const TYPE entry = TYPE());

  //! Copy constructor
  Vector(const Vector<TYPE> & vec);

  //! Destructor
  ~Vector();  

  //! Hard coded query if values are complex
  Boolean IsComplex();
  
  //! Return a Double pointer to the data of the vector.
  //! If the Vector is complex, a new array is created,
  //! where the entries are sequentilly ordered in real
  //! and imaginary parts (real_1, imag_1, real_2, ...)
  Double* GetDoublePointer();

  //! Initalizes the vector with a given entry
  /*!
    \param entry (input) Entry vector gets inialized with
  */
  //! \note this method does not change the size of the vector!
  void Init(const TYPE entry = TYPE());

  //! Get the length of the vector
  inline Integer GetSize() const {return size_;}

  //! Set the lenght of the vector
  /*!
    \param size (input) Lengh of vector
  */
  //! \note the entries are set to zero afterwards!
  void Resize(const Integer size);
 
  //! Set the entry i of the vector to the given value (x[i] = s)
  /*!
    \param i (input) Index of entry s
    \param s (input) Entry to be set on position i
  */
  void SetEntry(const Integer i, const TYPE &s);
  
  //! Get the entry i of the vector on the given value (ret = x[i])
  /*!
    \param i (input) Index of entry s
    \param ret (output) Entry on position i
  */
  void GetEntry(const Integer i, TYPE &ret) const;
  
  //! Add s to i-th vector entry (x[i] += s)
  /*!
    \param i (input) Index of entry s
    \param s (input) Value to be added to x[i]
  */
  void AddEntry(const Integer i, const TYPE &s);
  
  //! Mult the i-th vector entry with s (x[i] *= s)
  /*!
    \param i (input) Index of entry s
    \param s (input) Factor, which i-the entry gets multiplied with
  */
  void MultEntry(const Integer i, const TYPE &s);
  
  //! Multiply the i-th vector entry with a and add s on it (x[i] = a*x[i]+s)
  /*!
    \param i (input) Index of entry s
    \param a (input) Factor the i-the entry gets multiplied with
    \param s (input) Value to be added to a*x[i]
  */
  void MultAddEntry(const Integer i, const TYPE &a, const TYPE &s);

  //! Adds another vector to itself (x = x+y)
  /*! 
    \param y (input) Addend to the vector
  */
  void Add(const CFSVector& y);

  //! Adds the multiple of another vector to itself (x = x +a*y)
  /*!
    \param a (input) Factor for scaling vector y
    \param y (input) Addend to the vector
  */
  void Add(const TYPE a, const CFSVector &y);

  //! Replaces the vector by the sum of two scaled vectors (x = a*y+b*z)
  /*!
    \param a (input) Factor for vector y
    \param y (input) Vector scaled by factor a
    \param b (input) Factor for vector z
    \param z (input) Vector scaled by factor b
  */
  void Add( const TYPE a, const CFSVector& y,
		    const TYPE b, const CFSVector& z );

  
  //! Scales the vector itself and adds the multiple of another one y (x = a*x + y)
  /*!
    \param a (input) Factor for scaling the vector itself
    \param y (input) Addend for the vector (gets not scaled)
  */
  void Axpy(const TYPE a, const CFSVector &y);
 
  
  //! Performes the dot/inner product of the vector itself with y (=x^T*y)
  /*! 
    \param y (input) Vector to perform innter product with
    \param result (output) Result of x^T * y
  */
  void Inner(const CFSVector &y, TYPE &result) const;
  
  //! Calculates the Eucluidean L2-norm
  Double NormL2() const;
  

  //*************************************************
  //* old interface which is compatible to previous *
  //* version of Vector<TYPE>                       *
  //*************************************************

    //! Overloading of operation =
  Vector	&operator=	(const Vector &);

    //! Assignment operator for base class
  CFSVector & operator= (const CFSVector & vec);


  //! Element can be referred to as v[i]
  inline TYPE	&operator[]	(const Integer i) const
  {	
#ifdef CHECK_INDEX
     if (i >= size_){
      std::string errorMsg;
      errorMsg =  "Vector: invalid access to element ";
      errorMsg += Info->GenStr(i+1);
      errorMsg += "! \n Length of vector: ";
      errorMsg += Info->GenStr(size_);
      Error(errorMsg.c_str(),__FILE__,__LINE__);
    }
#endif
    return  data_[i];
  }

  //! Return pointer p to array 
  TYPE*  GetPointer() const
  {
#ifdef CHECK_MEMORY
    if (!data_)
      Error("Vector: undefined Vector",__FILE__,__LINE__);
#endif
    return data_;
  }

  //! initialize vector with defined data
  /*!
    \param nsize size of vector
    \param ptdata pointer to array with values of vector
  */
  void TransformInVector(const Integer nsize, TYPE* ptdata);
  
  //! Overloading of operations +,+=
  Vector	operator+	(const Vector &) const;
  Vector	&operator+=	(const Vector &);

  //! Overloading of operations -,-=
  Vector        operator-       () const;
  Vector	operator-	(const Vector &) const;
  Vector	&operator-=	(const Vector &);

  //! Overloading of operations / for number and Vector
  Vector  operator/       (const TYPE &) const;

  //! Overloading of operation not equal for Vector
  Vector &operator/= (const TYPE &x);

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

  //! Overloading of assignement operatior for Vector and Matrix  
  Vector	operator=	(const Matrix<TYPE> &) const;

  //! Return part of Vector from index i to ii
  Vector	Part	(const Integer i, const Integer ii) const;

  //! Constructs the unit vector of length n, which only non-zero
  //! entry is a 1 at the i-th position
  /*! 
    \param n (input) length of vector
    \param i (input) component, which is 1
    \f[ \left( \begin{array}{c} 0  \\ \cdots \\ 0 \\ 1 \\ 0 \\ \cdots \\ 0 
    \end{array} \right) \f]
  */
  static Vector	Unit	(const Integer n, const Integer i);

  //! Add element of the same type at position pos, by default to the beginning Beware of numeration in C++
  void AddElement(const TYPE & y, Integer pos=0);

  //! Add element of the same type at the end of the vector
  void Push_back(const TYPE & y);

  //! Insert a vector to this vector at position pos
  void InsertVector (const Vector<TYPE> & y, Integer pos=0);

  //! Delete element from vector on position pos
  void  Cut (const Integer pos);

  //! Delete elements from position pos1 to pos2, on pos1, pos2 too
  void  Cut (const Integer pos1, const Integer pos2);

  //! Return size of space memory of this vector
  Integer  Memory() const;

  //!Friends
  //! Sort of vector: v - vec.p, n - vec.size
  template <class S> void Sort(S* v, Integer n);
 
 //! Swap 2 elements in vector Ex Swap(v[i],v[j])
  template<class T2> void Swap(T2& a, T2 & b);

 protected:

  //! Length of the vector
  Integer size_;

  //! Data of the vector
  TYPE* data_;
  
  //! Capacity of the vector
  Integer capacity_;

};

// ******************************************************
// * Additional functions related with handling vectors *
// ******************************************************
  
//! Function for swap  number a and b (Integer, Double)
template <class S> void Swap(S &, S &);

//! Function for swap Vector<T>;
template<class TYPE>
void Swap(Vector<TYPE> & a, Vector<TYPE> & b);

//! Overloading << for class vector
template<class TYPE>  std::ostream& operator << ( std::ostream & , const Vector<TYPE> &);


#ifdef __GNUC__
  // Template instantiation for used vectors
  template class Vector<Integer>;
  template class Vector<Double>;
  template class Vector<Complex>;
#endif
 
} // end of namespace

#endif
