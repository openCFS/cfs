#ifndef OLAS_VECTOR_HH
#define OLAS_VECTOR_HH

#include "TypeDefs.hh"
#include "SingleVector.hh"


#include <def_build_type_options.hh>
#include <def_use_embedded_python.hh>

#include "Utils/tools.hh"
#include "Utils/boost-serialization.hh"

#ifdef USE_EXPRESSION_TEMPLATES
#include "exprt/xpr1.hh"
#else
#include "promote.hh"
#endif

#ifdef USE_EMBEDDED_PYTHON
  #define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
  #include <Python.h>
#else
  struct PyObject;
#endif

namespace CoupledField {

// !Forward class declarations
template<typename T> class Matrix;
template<typename T> class Vector;
template<typename T> class NodeStoreSol;
template<typename T> class ElemStoreSol;

  //! Class for dense array-based algebraic vector
  template <typename T>
  
#ifdef USE_EXPRESSION_TEMPLATES
  class Vector : public SingleVector, public Dim1<T, Vector<T> >
#else 
  class Vector : public SingleVector 
#endif
  {

  public:
  
    //! Friend declaration for Matrix
    friend class Matrix<T>;
    //! Friend declaration for NodeStoreSol
    friend class NodeStoreSol<T>;
    //! Friend declaration for ElemStoreSol
    friend class ElemStoreSol<T>;

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION AND REPLACEMENT
    // =======================================================================

    //! \name Methods for construction, destruction and replacement
    //! The vector class can also be used as a form of container. We
    //! provide methods that allow to add to a simple data array the
    //! functionality of this vector class, especially its arithmetic and
    //! access operations.

    //@{

    //! Default Constructor
    // We set this to true by default, to be able to do a re-size
    // afterwards. If a Replace() is done, it will adapt it accordingly.
    // size_ is set to 0 in ctor of SingleVector
    Vector() : SingleVector(), data_(NULL), capacity_(0), memBelongsToMe_(true)
    {  
    }

    //! Constructor with initial size

    //! Creates an vector with size \a size and gets initialized
    //! with \a entry. If no \a entry is provided, the vector
    //! is initialized with zeroes.
    Vector(const UInt size, const T entry = 0);

    //! Copy Constructor

    //! This is a deep copy constructor. It will allocate its own data_
    //! array and generate an actual copy of the entries of the original
    //! vector. Consequently it sets the memBelongsToMe_ attribute to true.
    //! This constructor is mainly provided to allow for certain sparse
    //! arithmetic operations of the SBM_Vector class.
    Vector(const Vector<T> &origVec)
    {
      // Obtain size info and allocate memory
      size_     = origVec.size_;
      capacity_ = size_;
      memBelongsToMe_ = true;

      data_ = size_ > 0 ? new T[size_] : NULL;

      std::copy_n(origVec.data_, size_, data_);
    }

    /** constructs a vector out of two vectors with the size of both vectors */
    Vector(const Vector<T>& lower, const Vector<T>& upper)
    {
      // Obtain size info and allocate memory
      size_     = lower.size_ + upper.size_;
      capacity_ = size_;
      memBelongsToMe_ = true;

      data_ = size_ > 0 ? new T[size_] : NULL;

      for(unsigned int i = 0; i < lower.size_; ++i)
        data_[i] = lower.data_[i];
      for(unsigned int i = 0; i < upper.size_; ++i)
        data_[lower.size_ + i] = upper.data_[i];
    }

    /** construct a data with external data. Either copy the date or be a frontend for it.
     * @param entries when transferMem note, that the const is casted away
     * @param copy if false point to external data and refuse resizing and deletion */
    Vector(UInt length, const T* entries, bool copy)
    {
      memBelongsToMe_ = copy;
      // if copy we start with 0 and copy afterwards via Fill()
      size_ = !copy ? length : 0;
      capacity_ = size_;
      data_ = !copy ? const_cast<T*>(entries) : NULL;

      if(copy)
        Fill(entries, length);
    }

     Vector(UInt length, T* entries, bool copy)
     {
       memBelongsToMe_ = copy;
       // if copy we start with 0 and copy afterwards via Fill()
       size_ = !copy ? length : 0;
       capacity_ = size_;
       data_ = !copy ? entries : NULL;

       if(copy)
         Fill(entries, length);
     }

    /** construct a vector from a numpy array vector (dim=1).
     * Creates own data and copies the content.
     * @param obj numpy array which is validated
     * @param decref shall the object referecence counter be decremented */
    Vector(PyObject* obj, bool decref);

    //! Destructor
    
    //! The default destructor must be deep, i.e. it must free all dynamically
    //! allocated memory.
    virtual ~Vector();

    //! Return the Entry type of the vector

    //! The method returns the entry type of the vector (i.e. Double, Complex,
    //! or whatever it is). This is encoded as a value of the enumeration data
    //! type MatrixEntryType.
    BaseMatrix::EntryType GetEntryType() const
    {
      return  CoupledField::EntryType<T>::M_EntryType;
    }

    //! Re-size the vector

    //! This method can be used to change the length of the vector. Data might be lost!
    //! \param newSize the new length of the vector
    //! \note
    //! - When we grow, all data will be lost - no copy of data as for StdVector()!
    //! - A re-allocation of memory will only be triggered by Resize(),
    //!   if the new length of the vector exceeds the length of the
    //!   internal data array as given by dataSize_.
    //! - Re-size will currently refuse to perform a re-size operation,
    //!   if it is not responsible for the memory management of the data_
    //!   array.
    void Resize(const unsigned int newSize);
    
    //! Resize the vector to new size and initialize entries with val
    void Resize(const unsigned int newSize, const T val);

    unsigned int GetSize() const { return size_; }

    //! Add functionality of vector class to a data array

    //! This method allows to add the functionality of the Vector class,
    //! especially its arithmetic and access methods, to a plain data array.
    //! Calling this method will replace the internal data_ array by the
    //! entries vector and re-set the internal attributes, like e.g. size_.
    //! Note that, since the vector is one-based, the pointer to the new
    //! memory array must have an appropriate off-set and that length may
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
    void Replace(UInt length, T* entries, bool transferMem);

    //! Withdraw responsibility for memory management from vector object

    //! Calling this method will relieve the vector object from the task of
    //! managing the memory block containing the entries of the vector.
    //! It will also return a (one-based) pointer to that memory block.
    //! \note Responsibility for de-allocating the memory should belong to
    //!       a unique object. Therefore, DecoupleMem will issue a warning
    //!       when an attempt is made to obtain that responsibilty from the
    //!       vector object and the latter does not own it.
    T* DecoupleMem()
    {
      if(memBelongsToMe_ == false)
      {
        WARN("DecoupleMem was called on a vector object not "
             << "responsible for managing its memory block! Memory "
             << "problems may arise!");
      }
      memBelongsToMe_ = false;
      return data_;
    }

    //! Clear the vector

    //! Calling this method will clear the vector, i.e. the internal data_
    //! array will be de-allocated, if it belongs to the vector object, and
    //! the internal attributes will be re-set to the state we also obtain
    //! from the default constructor.
    void Clear(bool keepCapacity = true) {
      size_ = 0;
      if( keepCapacity == false ) {
        if(memBelongsToMe_)
          delete[] data_;

        capacity_ = 0;
        data_     = NULL;
      }
      // We set this to true be default, to be able to do a re-size
      // afterwards. If a Replace() is done, it will adapt it accordingly.
      memBelongsToMe_ = true;
    }

    /** Fills the vector with data. Does a resize
     * The data is copied */
    void Fill(const T* source, unsigned int size);

    /** Fills the whole vector by increments
     * @param start this value will be at the first position
     * @param increment second = start + increment, ... */
    void Fill(const T& start, const T& increment);

    /** Resize and fill the object with a numpy array (dim=1) object
     * @see Export() */
    void Fill(PyObject* obj, bool decref = true);

    //@}

    // =======================================================================
    // ARITHMETIC OPERATIONS
    // =======================================================================

    //! \name Arithmetic operations

    //@{
    //! Add vec to this vector object
    void Add(const SingleVector &vec);

    //! Add a scaled version of a vector to this vector object

    //! The method takes this vector object \f$x\f$ and replaces it with
    //! the result of \f$x + \alpha v\f$.
    void Add(T a, const SingleVector &vec);

    //! Override SingleVector functions
    //    virtual void Add(Double a,const SingleVector &vec);
    //    virtual void Add(Complex a,const SingleVector &vec);

    //! Replace this vector object by the sum of two scaled vectors

    //! This method replaces this vector object by the sum
    //! \f$\alpha x +\beta y\f$.
    void Add(T a, const SingleVector &vec1,
             T b, const SingleVector &vec2);

    //! Override SingleVector functions
    //    virtual void Add(Double a, const SingleVector& vec1,
    //                     Double b, const SingleVector& vec2);
    //    virtual void Add(Complex a, const SingleVector& vec1,
    //                     Complex b, const SingleVector& vec2);


    /** this = a * vec. Like Add() with initial 0 */
    void Set(T a, const SingleVector &vec);


    /** Hadamard product. Pointwise this[i] = v1[i] * v2[i] */
    void Hadamard(const Vector<T>& v1, const Vector<T>& v2);

    //! Compute inner product

    //! The method computes the value of the inner product between this vector
    //! and the input vector vec. The value is returned in sum.
    void Inner(const SingleVector& vec, T& sum) const;
    
    /** computes the inner product only for a defined range
     * @param start first index
     * @param end not included index */
    T Inner(const SingleVector& vec, unsigned int start, unsigned int end) const;

    T Inner(const SingleVector& vec) const;

    T Inner() const;

    Double InnerAngle(const Vector<T>& other, bool preferPositiveZ=false) const;
    
    Vector<T> Conj() const;

    //! Override SingleVector functions
    //    virtual void Inner(const SingleVector& vec,Double& s) const;
    //    virtual void Inner(const SingleVector& vec,Complex& s) const;

    //! Compute Euclidean norm of this vector object
    virtual double NormL2() const;
    virtual double NormL2_squared() const;

    /** diff norm */
    double NormL2(const Vector<T>& other) const;

    /**  this functions localized the maximal component (absolute value) and returns it with its original sign
		example: SignedMax([1,0,0]) = 1; SignedMax([-1,0,0]) = -1 */ 
    virtual Double SignedMax() const;

    //** evaluate ModalAssuranceCriterion */
    Double MAC(const Vector<T>& other) const;

    /** Sum of all the vector's elements */
    T Sum() const;

    T Avg() const;

    /** Product of elements.
     * @see Inner()  */
    T Product() const;

    /** Extremal element. For Complex separate for real and imaginary part */
    T Min() const;
 
    /** some index of the smallest number. Not implemented yet for complex */
    unsigned int ArgMin() const;

    T Max() const;

    /** returns the entry with the maximum absolute value, the location is reported in loc */
    T MaxAbs(int& loc) const;

    /** return the minimal and maximal element concurrently.
     * @see Min() for complex */
    void MinMax(T& min, T& max) const;

    /** Calculates the max-norm (of the real part) */ 
    Double NormMax() const; 

    /** Calculates the max-norm (of the real part) of the vector difference */
    Double NormMax(const SingleVector& other) const;

    //! Normalize vector, i.e. scale it to unit length w.r.t. L2 norm.
    //! Returns the vectors L2 norm before normalization.
    Double Normalize();

    /** Count the number of non-zero entries */
    Integer CountNonZero() const;

    //! Create scalar product of vector and second one
    template <typename T2>
    PROMOTE(T,T2)  operator*( const Vector<T2> & x ) const;

    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the
    //! factor alpha and adds the vector y to it. The result will over-write
    //! the vector x.
    void Axpy( const T alpha, const SingleVector &y );

    //! Overload assignment operator
    SingleVector &operator= ( const SingleVector &stdvec );

    //! Divide each vector entry by specified real-valued scalar
    void ScalarDiv( const Double factor );

    //! Multiply each vector entry by specified real-valued scalar
    void ScalarMult( const Double factor );

    //! Divide each vector entry by specified complex-valued scalar
    void ScalarDiv( const Complex factor );

    //! Multiply each vector entry by specified complex-valued scalar
    void ScalarMult( const Complex factor );

    //! Calc cross product of two vectors v = this x b
    void CrossProduct( const Vector<T>& b, Vector<T>& v ) const;
    
    //! Is this vector collinear with another vector?
    bool Collinear( const Vector<T>& vec);
    
//@}


#ifdef USE_EXPRESSION_TEMPLATES
    
    // =======================================================================
    // INTERFACE TO EXPRESSION TEMPLATES
    // =======================================================================
    //! \name Interface To Expression Template Headers
    //@{ 
    
    //! Vector assignment operator using expression templates
    template <class X> Vector<T>&  operator=( const Xpr1<T,X>& rhs ) {
      return this->assignFrom(rhs);
    }
    //! Abstract vector assignment operator using expression templates
    template <class V> Vector<T>&  operator=( const Dim1<T,V>& rhs ) {
      return this->assignFrom(rhs);
    }

    //! Vector assignment operator
    Vector<T>&  operator=( const Vector<T>& rhs ) { 
      return this->assignFrom(rhs); 
    }
    
    //! vector product for general vector expression as second argument.
    //! Not that in the complex case not the complex conjugate is used! Use Inner() for that purpose!
    template <class V> 
    T operator*( const Xpr1<T,V>& rhs ) {
      T ret = 0.0;
      for( UInt i = 0; i < size_; ++i ) {
        ret += data_[i] * rhs(i);
      }
      return ret;
    }
    
    //! Returns the number of entries
    inline unsigned int size() const { return size_; }
    
    //! Access operator used for expression templates (const)
    inline T &  operator()( unsigned i ) { 
      return data_[i]; 
    }

    //! Access operator used for expression templates 
    inline T operator()(int i) const { 
      return data_[i]; 
    }
    
    //@}

#else

    // =======================================================================
    // MATHEMATICAL OPERATORS
    // =======================================================================

    //! \name Mathematical Operators
    //! \note Due to problems in Doxygen the binary operators +,-,*,/ 
    //!       (which use type promotion) are not shown, although they exist!

    //@{


    //! Assignment operator
    Vector<T> &operator=(const Vector<T> &x);

    //! Unary plus operator (this = +this)
    Vector<T> operator+() const;
    
    //! Create new vector by addition (new = this + y) (type promotion)
    template <typename T2>
    Vector<PROMOTE(T,T2)> operator+( const Vector<T2> & y ) const;
    
    //! Adds a second vector to own one (this += y)
    Vector<T> &operator+= ( const Vector<T> & y );
    
    //! Unary minus operator 
    Vector<T> operator-() const;
    
    //! Create new vector by subtraction (new = this - y) (type promotion)
    template <typename T2>
    Vector<PROMOTE(T,T2)> operator- ( const Vector<T2> &y ) const;

    //! Subtract a second vector from own one (this -= y)
    Vector<T> &operator-= ( const Vector<T> &y );
        
    //! Create new vector by division with scalar (new = this / y) 
    //! (type promotion)
    template <typename T2>
    Vector<PROMOTE(T,T2)> operator/( const T2 &y ) const;

    //! Divide entries of own vector by y
    Vector<T> &operator/= (T x);

    //! Create new vector by multiplication with scalar (new = this * y)
    //! (type promotion)
    template <typename T2>
    Vector<PROMOTE(T,T2)> operator* ( const T2 &y ) const;

    //! Multiply entries of own vector by y
    Vector<T> &operator*= (T x);
    //@}

#endif // USE_EXPRESSION_TEMPLATES


    //! Equality operator - outside USE_EXPRESSION_TEMPLATES
    bool operator==(const Vector<T> &x) const;

    /** comparison is done via memcmp */
    bool operator!=( const Vector<T>& x) const;



    // =======================================================================
    // I/O OPERATIONS
    // =======================================================================

    //! \name Methods for I/O operations

    /** print content for logging or nice output
     * @param format one TS_PLAIN (no brackets), TS_MATLAB, TS_PYTHON, TS_INFO (summary only)
     * @param sep separator string, when empty meaningful default is used
     * @param digits when given enforces format */
    std::string ToString(ToStringFormat format = TS_PYTHON, const std::string& sep="", int digits=-1) const;

    //! Export vector to file

    //! This method can be used to export the vector to an ascii file. The
    //! format is extremely simple. If the vector is of dimension
    //! \f$n\times 1\f$, then the output file will contain \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$a_k\f$.
    virtual void Export(const std::string& fname, BaseMatrix::OutputFormat format ) const;

    /** @see BaseVector::Import()*/
    void Import(const std::string& fname, bool checkSize=true);

    /** writes the content of the vector to a numpy array which needs to have proper size and type.
     * @see Fill() */
    void Export(PyObject* obj);


    // =======================================================================
    // OBTAIN / MANIPULATE MATRIX ENTRIES
    // =======================================================================

    //! \name Obtain/manipulate matrix entries
    //@{

    //! Initialize vector
    void Init( ) {
      Init( 0 );
    }

    //! Initialize vector with a given entry

    //! Initializes the vector with a given entry
    //! \param entry (input) Entry vector gets initialized with
    //! \note this method does not change the size of the vector!
    void Init( const T entry  );

    //! Return a reference to i-th entry
    inline T & operator[]( UInt i) {
#ifdef CHECK_INDEX
    if (i >= size_)
      EXCEPTION("Vector: invalid access to index " << i << " when size is " << size_);
#endif
      
      return data_[i];
    };

    //! Return a reference to i-th entry
    inline const T& operator[]( UInt i ) const {
#ifdef CHECK_INDEX
    if (i >= size_)
      EXCEPTION("Vector: invalid access to index " << i << " when size is " << size_);
#endif
      
      return data_[i];
    };

    //! Set the value of a vector entry

    //! This method sets the entry of the vector at position i to the
    //! specified value.
    void SetEntry( UInt i, const T &val );

    //! Get a reference to a vector entry

    //! This method sets val to the value of vector entry i
    void GetEntry( UInt i, T &val ) const;

    T Last() const;

    //! Get entries at specified indices from in
    const Vector<T> GetEntries( const StdVector<UInt>& in) const;

    //! Add val to the value of a vector entry

    //! This method adds to the entry of the vector at position i the
    //! specified value val.
    void AddToEntry( UInt i, const T &val );

    //! Add val to the value of a vector entry (atomic/thread-safe version)
    void AddToEntryAtomic( UInt i, const T &val );

    //! Add element of the same type at the end of the vector
    void Push_back( const T & y );

    //! Return a special part ( real, imag, amplitude, phase) of a vector
    Vector<Double> GetPart( Global::ComplexPart part ) const;

    //! Set special part ( real, imag, amplitude, phase) of a vector
    void SetPart( Global::ComplexPart part, const Vector<Double> & partVector,
                  bool zeroOtherPart = false );

    //@}

    // =======================================================================
    // DIRECT ACCESS TO INTERNAL DATA STRUCTURE
    // =======================================================================

    //! \name Direct access to internal data structure
    //! The following methods allow a direct manipulation (read or write) of
    //! the internal data structure. Though this does violate the concept of
    //! data encapsulation, these methods are provided for the sake of
    //! efficiency.

    //@{

    //! Return pointer to internal data array data_
    T* GetPointer(){
      return data_;
    }

    //! Return pointer to internal data array data_ (const)
    const T* GetPointer() const{
      return data_;
    }

    /** Check if the vector contains the specific value. A way to assert() that is is completely set */
    bool Contains(const T val) const;

    /** Check if the vector contains NAN. To be used by asserts() */
    bool ContainsNaN() const;

    /** Check if the vector contains +/- INF. To be used by asserts() */
    bool ContainsInf() const;

    //@}

    // the following code makes Vector iterable by ranges for(double v : vec)
    T* begin() { return data_; }
    T* end()   { return data_ + size_; }

    const T* begin() const { return data_; }
    const T* end()   const { return data_ + size_; }


  protected:

    //! 0-based array storing the vector entries.
    T *data_;
    
    //! capacity of the vector
    UInt capacity_;

    //! Flag signaling whether management of data array is done by this object

    //! This attribute is used to keep track on the fact whether the object
    //! is responsible for managing the memory of the data_ array, especially
    //! its deallocation.
    bool memBelongsToMe_;

  };
  
  // ***********************************************************
  // * Additional free functions related with handling vectors *
  // ***********************************************************
  
  //! Overloading << for class vector
  template<typename T>  std::ostream& operator << ( std::ostream & , 
                                                    const Vector<T> &);


  // *******************************************
  //   Vector product, but not the scalar product, which would use the complex conjugate as Inner() does!
  // *******************************************
  template<typename T> template<typename T2>
  PROMOTE(T,T2) Vector<T>::
  operator* (const Vector<T2> &x) const
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.GetSize() == 0))
      EXCEPTION( "Vector: undefined Vector in operator *(vector)" );
#endif  

#ifdef CHECK_INDEX
    if (size_ != x.GetSize())
      EXCEPTION( "Vector: incompatible dimension for operator *(vector)" );
#endif

    PROMOTE(T,T2) ret;
    ret = data_[0] * x[0];
    for (UInt i = 1; i < size_; i++)
      ret += data_[i] * x[i];

    return ret;
  }


  // ************************************************************
  //  INLINE MEMBER DEFINITIONS FOR NON-TEMPLATE EXPRESSION CASE
  // ************************************************************
#ifndef USE_EXPRESSION_TEMPLATES
  
  template<typename T> template<typename T2>
  Vector<PROMOTE(T,T2)> Vector<T>::
  operator+(const Vector<T2> &x) const
  {       

#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.GetSize() == 0))
      EXCEPTION( "Vector: undefined Vector in operator +(vector)" );
#endif  

#ifdef CHECK_INDEX
    if (size_ != x.GetSize())
      EXCEPTION( "Vector: incompatible dimension for operator +(vector)" );
#endif

    Vector<PROMOTE(T,T2)> ret(size_);

    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] + x[i];

    return ret;
  }

  template<typename T> template<typename T2>
  Vector<PROMOTE(T,T2)> Vector<T>::
  operator-(const Vector<T2> &x) const
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.GetSize() == 0))
      EXCEPTION( "Vector: undefined Vector in operator -(vector)" );
#endif  

#ifdef CHECK_INDEX
    if (size_ != x.GetSize())
      EXCEPTION( "Vector: incompatible dimension for operator -(vector)" );
#endif
    Vector<PROMOTE(T,T2)> ret(size_);

    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] - x[i];

    return ret;
  }

  template<typename T> template<typename T2>
  Vector<PROMOTE(T,T2)> Vector<T>::
  operator* (const T2 &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION( "Vector: undefined Vector in operator *(number)" );
#endif

    Vector<PROMOTE(T,T2)> ret(size_);

    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] * x;

    return ret;
  }

  template<typename T> template<typename T2>
  Vector<PROMOTE(T,T2)> Vector<T>::
  operator/(const T2 &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION( "Vector: undefined Vector in operator /(number)" );
#endif

    Vector<PROMOTE(T,T2)> ret(size_);

    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] / x;

    return ret;
  }
#endif // EXPR_TEMPLATE
}


#endif // OLAS_VECTOR_HH
