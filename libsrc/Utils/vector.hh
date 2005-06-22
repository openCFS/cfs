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
  template <class TYPE> class VectorListInitializer;

  //! Concrete Template class for a general dense vector
  template<class TYPE>
  class Vector : public CFSVector {
  public:

    // Friend declarations
    friend class Matrix<TYPE>;
    friend class NodeStoreSol<TYPE>;
    friend class ElemStoreSol<TYPE>;

    //! Constructor
    Vector();

    //! Constructor with inital size.
    //! All entries are filled with zeroes
    Vector(const UInt size, const TYPE entry = TYPE());

    //! Copy constructor
    Vector(const Vector<TYPE> & vec);

    //@{
    //! Converting a point to a vector
    Vector(const Point<2> & p);
    Vector(const Point<3> & p);
    //@}

    //! Destructor
    ~Vector();  

    //! Hard coded query if values are complex
    Boolean IsComplex();
  
    //! Deletes the content of the vector
    void Clear();

    //! Initalizes the vector with a given entry
    /*!
      \param entry (input) Entry vector gets inialized with
    */
    //! \note this method does not change the size of the vector!
    void Init(const TYPE entry = TYPE());

    //! Add functionality of vector class to a data array
  
    //! This method allows to add the functionality of the Vector class,
    //! especially its arithmetic and access methods, to a plain data array.
    //! Calling this method will replace the internal data_ array by the
    //! entries vector and re-set the internal attributes, like e.g. size_.
    //! Note that, the length may
    //! not be larger than the actual length of the allocated memory block.
    //! Calling Replace is memory safe in the sense that the old data_ array
    //! will be de-allocated, if this is the responsibility of the
    //! corresponding vector object. Responsibility for de-allocating the
    //! new data_ array is transferred to the vector object by setting the
    //! transferMem parameter to true.
    //! \param length      the length of the new vector
    //! \param entries     pointer (one-based) to the data array containing
    //!                    the new vector entries
    //! \param transferMem flag signalling transfer of responsibility for
    //!                    memory management
    void Replace( UInt length, TYPE* entries, Boolean transferMem );

    //! Get the length of the vector
    inline UInt GetSize() const {return size_;}

    //! Set the lenght of the vector
    /*!
      \param size (input) Lengh of vector
    */
    //! \note the entries are set to zero afterwards!
    void Resize(const UInt size);
 
    //! Set the entry i of the vector to the given value (x[i] = s)
    /*!
      \param i (input) Index of entry s
      \param s (input) Entry to be set on position i
    */
    void SetEntry(const UInt i, const TYPE &s);
  
    //! Get the entry i of the vector on the given value (ret = x[i])
    /*!
      \param i (input) Index of entry s
      \param ret (output) Entry on position i
    */
    void GetEntry(const UInt i, TYPE &ret) const;
  
    //! Add s to i-th vector entry (x[i] += s)
    /*!
      \param i (input) Index of entry s
      \param s (input) Value to be added to x[i]
    */
    void AddEntry(const UInt i, const TYPE &s);
  
    //! Mult the i-th vector entry with s (x[i] *= s)
    /*!
      \param i (input) Index of entry s
      \param s (input) Factor, which i-the entry gets multiplied with
    */
    void MultEntry(const UInt i, const TYPE &s);
  
    //! Multiply the i-th vector entry with a and add s on it (x[i] = a*x[i]+s)
    /*!
      \param i (input) Index of entry s
      \param a (input) Factor the i-the entry gets multiplied with
      \param s (input) Value to be added to a*x[i]
    */
    void MultAddEntry(const UInt i, const TYPE &a, const TYPE &s);

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
    Vector        &operator=      (const Vector &);

    //! Assignment operator for base class
    CFSVector & operator= (const CFSVector & vec);


    //! Element can be referred to as v[i]
    inline TYPE   &operator[]     (const UInt i) const
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
    void TransformInVector(const UInt nsize, TYPE* ptdata);
  
    //! Overloading of operations +,+=
    Vector        operator+       (const Vector &) const;
    Vector        &operator+=     (const Vector &);

    //! Overloading of operations -,-=
    Vector        operator-       () const;
    Vector        operator-       (const Vector &) const;
    Vector        &operator-=     (const Vector &);

    //! Overloading of operations / for number and Vector
    Vector  operator/       (const TYPE &) const;

    //! Overloading of operation not equal for Vector
    Vector &operator/= (const TYPE &x);

    //! Overloading of operations * for number and Vector   
    Vector        operator*       (const TYPE &) const;

    //! Overloading of operations * for Vector and Vector 
    TYPE  operator*       (const Vector &) const;

    //! Overloading of operations * for Vector and Matrix  
    Vector        operator*       (const Matrix<TYPE> &) const;

    //! Overloading of operation *= for Vector and number
    Vector        &operator*=     (const TYPE &);

    //! Overloading of operations * for Vector and Matrix 
    Vector        &operator*=     (const Matrix<TYPE> &);

    //! Overloading of operation equal for Vector
    Boolean       operator==      (const Vector &) const;

    //! Overloading of operation not equal for Vector
    Boolean       operator!=      (const Vector &) const;

    //! Overloading of assignement operatior for Vector and Matrix  
    Vector        operator=       (const Matrix<TYPE> &) const;

    //! This method is needed to intialize a StdVector like this
    //! StdVector<Integer> A;<br>
    //! A = 1,2,5,10
    inline VectorListInitializer<TYPE> operator=(const TYPE x);

    //! Return part of Vector from index i to ii
    Vector        Part    (const UInt i, const UInt ii) const;

    //! Constructs the unit vector of length n, which only non-zero
    //! entry is a 1 at the i-th position
    /*! 
      \param n (input) length of vector
      \param i (input) component, which is 1
      \f[ \left( \begin{array}{c} 0  \\ \cdots \\ 0 \\ 1 \\ 0 \\ \cdots \\ 0 
      \end{array} \right) \f]
    */
    static Vector Unit    (const UInt n, const UInt i);

    //! Add element of the same type at position pos, by default to the beginning Beware of numeration in C++
    void AddElement(const TYPE & y, UInt pos=0);

    //! Add element of the same type at the end of the vector
    void Push_back(const TYPE & y);

    //! Insert a vector to this vector at position pos
    void InsertVector (const Vector<TYPE> & y, UInt pos=0);

    //! Delete element from vector on position pos
    void  Cut (const UInt pos);

    //! Delete elements from position pos1 to pos2, on pos1, pos2 too
    void  Cut (const UInt pos1, const UInt pos2);

    //! Return size of space memory of this vector
    UInt  Memory() const;

    //! Sort vector in ascending order
    void Sort();
 
  protected:

    //! Length of the vector
    UInt size_;

    //! Data of the vector
    TYPE* data_;
  
    //! Capacity of the vector
    UInt capacity_;

    //! Flag signaling whether management of data array is done by this object
  
    //! This attribute is used to keep track on the fact whether the object
    //! is responsible for managing the memory of the data_ array, especially
    //! its deallocation.
    Boolean memBelongsToMe_;
  
  };

  // ******************************************************
  // HELPER CLASS FOR INITALIZING Vector
  // (ref. 'Techniques for Scientific C++' 
  // by Todd Veldhuizen, page. 43ff
  // ******************************************************
  template <class TYPE>
  class VectorListInitializer
  {
  public:

    //! Constructor
    VectorListInitializer(Vector<TYPE> * vec)
      :vec_(vec) {};
    
    //! Overloading of comma operator
    VectorListInitializer<TYPE> operator, (const TYPE x)
    {
      vec_->Push_back(x);
      return VectorListInitializer<TYPE>(vec_);
    }

  private:
    //! pointer to vector
    Vector<TYPE> * vec_;
  };

  //! 
  template<class TYPE>
  VectorListInitializer<TYPE> Vector<TYPE>::operator= (const TYPE x)
  {
    Clear();
    Push_back(x);
    return VectorListInitializer<TYPE>(this);
  }

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


#if defined(__GNUC__)
  // Template instantiation for used vectors
  template class Vector<Integer>;
  template class Vector<Double>;
  template class Vector<UInt>;
  template class Vector<Complex>;
#endif
 
} // end of namespace

#endif
