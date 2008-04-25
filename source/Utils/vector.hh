// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_VECTOR_2004
#define FILE_VECTOR_2004

#include "Utils/boost-serialization.hh"
#include "cfsvector.hh"
#include "DataInOut/WriteInfo.hh"
#include "promote.hh"
#include "General/exception.hh"
#include <def_build_type_options.hh>

#ifdef EXPR_TEMPLATES
#include "exprt/xpr1.hh"
#endif

namespace CoupledField {
  
  
  // !Forward class declarations
  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;
  template<class TYPE> class NodeStoreSol;
  template<class TYPE> class ElemStoreSol;
 
  //! Concrete Template class for a general dense vector
  template<class TYPE>
#ifdef EXPR_TEMPLATES
  class Vector : public CFSVector, public Dim1<TYPE, Vector<TYPE> >
#else
  class Vector : public CFSVector
#endif
  {
  public:
    
    //! Friend declaration for Matrix
    friend class Matrix<TYPE>;
    //! Friend declaration for NodeStoreSol
    friend class NodeStoreSol<TYPE>;
    //! Friend declaration for ElemStoreSol
    friend class ElemStoreSol<TYPE>;
    

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION, INITIALIZATION, RESIZING
    // =======================================================================
        
    //! \name Construction, Destruction, Initialization and Resizing

    //@{ 
    //! Default Constructor

    //! Creates an empty 0-sized vector
    Vector<TYPE>();

    //! Constructor with initial size
    
    //! Creates an vector with size \a size and gets initialized
    //! with \a entry. If no \a entry is provided, the vector
    //! is initialized with zeroes.
    explicit Vector<TYPE>( const UInt size, const TYPE entry = TYPE() );

    //! Copy constructor
    Vector( const Vector<TYPE> & vec );

    //! Converting a 3-dimensional point to a vector
    Vector( const Point & p );

    //! Destructor
    virtual ~Vector();  
    
    //! Set the length of the vector

    //! Set the length of the vector
    //! \param size (input) Length of vector
    //! \note the entries are set to zero afterwards!
    void Resize( const UInt size );
    
    //! Deletes the content of the vector
    void Clear();
    
    //! Initialize vector with a given entry

    //! Initializes the vector with a given entry
    //! \param entry (input) Entry vector gets initialized with
    //! \note this method does not change the size of the vector!
    void Init( const TYPE entry = TYPE() );
    
    //@}

    // =======================================================================
    // GENERAL INFORMATION
    // =======================================================================
    
    //! \name General Vector Information
    
    //@{

    //! Get entry type of vector
    EntryType::ScalarType GetEntryType() const {
      return  EntryTypeMap<TYPE>::S_TYPE;
    }

    //! Get the length of the vector
    inline UInt GetSize() const { return size_; }
    
    //! Return size of space memory of this vector
    UInt Memory() const;

    /** converts the content to a string vector */
    void ToString(StdVector<std::string>& out) const;
    
    /** Dumps the vector to a string - is slow, use only for debugging! */
    std::string ToString() const;

    
    
    //@}
    

    // =======================================================================
    // OBTAIN / MANIPULATE VECTOR ENTRIES
    // =======================================================================
    
    //! \name Obtain / Manipulate Vector Entries

    //@{
    
    //! General access operator (const)

    //! Returns the \a i-th entry of the vector
    //! \param i (input) Vector index
    inline TYPE  operator[] ( const UInt i ) const; 
    
    //! General access operator

    //! Returns the \a i-th entry of the vector
    //! \param i (input) Vector index
    inline TYPE &operator[] ( const UInt i );
      
    //! Return pointer p to plain c-array 
    inline TYPE* GetPointer() const;
    
    /** gives the the pointer to the plain c array */
    virtual void GetPointer(TYPE* &ptr_out) const;
    

    //! Add a new element at position \a pos

    //! Add element of the same type at position pos, 
    //! by default to the beginning 
    void AddElement( const TYPE & y, UInt pos = 0 );

    //! Add element of the same type at the end of the vector
    void Push_back( const TYPE & y );
  
    
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
    //! \param transferMem flag signaling transfer of responsibility for
    //!                    memory management
    void Replace( UInt length, TYPE* entries, bool transferMem );

    //! Set the entry i of the vector to the given value (x[i] = s)
   
    //! Set the entry i of the vector to the given value (x[i] = s)
    //! \param i (input) Index of entry s
    //! \param s (input) Entry to be set on position i
    void SetEntry( const UInt i, const TYPE &s );

    //! Get the entry i of the vector on the given value (ret = x[i])  

    //! Get the entry i of the vector on the given value (ret = x[i])
    //! \param i (input) Index of entry 
    //! \param ret (output) Entry on position i
    void GetEntry( const UInt i, TYPE &ret ) const;


    //! Return a special part ( real, imag, amplitude, phase) of a vector
    Vector<Double> GetPart(  DataType part ) const;

    //! Set special part ( real, imag, amplitude, phase) of a vector
    void SetPart( DataType part, const Vector<Double> & partVector );

    /** Fills the vector by the external content. Adjustes 
     * the size to confirm by the content! */
    void FillVector(const TYPE* data, unsigned int size);

    //! Add s to i-th vector entry (x[i] += s)

    //! Add s to i-th vector entry (x[i] += s)
    //! \param i (input) Index of entry s
    //! \param s (input) Value to be added to x[i]
    void AddEntry( const UInt i, const TYPE &s );
  
    //! Mult the i-th vector entry with s (x[i] *= s)
    
    //! Mult the i-th vector entry with s (x[i] *= s)
    //! \param i (input) Index of entry s
    //! \param s (input) Factor, which i-the entry gets multiplied with
    void MultEntry( const UInt i, const TYPE &s );
  
    //! Multiply the i-th vector entry with a and add s on 
    
    //! Multiply the i-th vector entry with a and add s on 
    //! it (x[i] = a*x[i]+s)
    //! \param i (input) Index of entry s
    //!  \param a (input) Factor the i-the entry gets multiplied with
    //!  \param s (input) Value to be added to a*x[i]
    void MultAddEntry( const UInt i, const TYPE &a, const TYPE &s );

    //@}
    
    // =======================================================================
    // NAMED ARITHMETIC OPERATIONS
    // =======================================================================
    
    //! \name Named Arithmetic Operations
    //@{

    //! Adds another vector to itself (x = x+y)

    //! Adds another vector to itself (x = x+y)
    //! \param y (input) Addend to the vector
    void Add( const CFSVector& y );

    //! Adds the multiple of another vector to itself (x = x +a*y)

    //! Adds the multiple of another vector to itself (x = x +a*y)
    //! \param a (input) Factor for scaling vector y
    //! \param y (input) Addend to the vector
    void Add( const TYPE a, const CFSVector &y );

    //! Replaces the vector by the sum of two scaled vectors (x = a*y+b*z)

    //! Replaces the vector by the sum of two scaled vectors (x = a*y+b*z)
    //! \param a (input) Factor for vector y
    //! \param y (input) Vector scaled by factor a
    //! \param b (input) Factor for vector z
    //! \param z (input) Vector scaled by factor b
    void Add( const TYPE a, const CFSVector& y,
              const TYPE b, const CFSVector& z );
    
    //! Scales the vector itself and adds the multiple of another one y 

    //! Scales the vector itself and adds the multiple of another one y 
    //! (x = a*x + y)
    //! \param a (input) Factor for scaling the vector itself
    //! \param y (input) Addend for the vector (gets not scaled)
    void Axpy( const TYPE a, const CFSVector &y );
 
    //! Calculate inner product of vector

    //! Performes the dot/inner product of the vector itself with y (=x^T*y)
    //! \param y (input) Vector to perform inner product with
    //! \param result (output) Result of x^T * y
    void Inner( const CFSVector &y, TYPE &result ) const;
  
    //! Calculates the Euclidean L2-norm
    Double NormL2() const;

    //! Return vector as separated string
    std::string Serialize( Char separator = ',') const;

    //! Create scalar product of vector and second one
    template <class TYPE2>
    PROMOTE(TYPE,TYPE2)  operator*( const Vector<TYPE2> & x ) const;
    //@}
 

#ifdef EXPR_TEMPLATES
    
    // =======================================================================
    // INTERFACE TO EXPRESSION TEMPLATES
    // =======================================================================
    //! \name Interface To Expression Template Headers
    //@{ 
    
    //! Vector assignment operator using expression templates
    template <class X> Vector<TYPE>&  operator=( const Xpr1<TYPE,X>& rhs ) {
      return assignFrom(rhs);
    }
    //! Abstract vector assignment operator using expression templates
    template <class V> Vector<TYPE>&  operator=( const Dim1<TYPE,V>& rhs ) {
      return assignFrom(rhs);
    }

    //! Vector assignment operator
    Vector<TYPE>&  operator=( const Vector<TYPE>& rhs ) { 
      return assignFrom(rhs); 
    }
    
    //! Returns the number of entries
    inline unsigned int size() const { return size_; }
    
    //! Access operator used for expression templates (const)
    inline TYPE &  operator()( unsigned i ) { 
      return data_[i]; 
    }

    //! Access operator used for expression templates 
    inline TYPE operator()(int i) const { 
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
    Vector<TYPE>  &operator=( const Vector<TYPE> & y );

    //! Unary plus operator (this = +this)
    Vector<TYPE> operator+() const;
    
    //! Create new vector by addition (new = this + y) (type promotion)
    template <class TYPE2>
    Vector<PROMOTE(TYPE,TYPE2)> operator+( const Vector<TYPE2> & y ) const;
    
    //! Adds a second vector to own one (this += y)
    Vector<TYPE> &operator+= ( const Vector<TYPE> & y );
    
    //! Unary minus operator 
    Vector<TYPE> operator-() const;
    
    //! Create new vector by subtraction (new = this - y) (type promotion)
    template <class TYPE2>
    Vector<PROMOTE(TYPE,TYPE2)> operator- ( const Vector<TYPE2> &y ) const;

    //! Subtract a second vector from own one (this -= y)
    Vector<TYPE> &operator-= ( const Vector<TYPE> &y );
        
    //! Create new vector by division with scalar (new = this / y) 
    //! (type promotion)
    template <class TYPE2>
    Vector<PROMOTE(TYPE,TYPE2)> operator/( const TYPE2 &y ) const;

    //! Divide entries of own vector by y
    Vector<TYPE> &operator/= ( const TYPE &y );

    //! Create new vector by multiplication with scalar (new = this * y)
    //! (type promotion)
    template <class TYPE2>
    Vector<PROMOTE(TYPE,TYPE2)> operator* ( const TYPE2 &y ) const;

    //! Multiply entries of own vector by y
    Vector<TYPE> &operator*= (const TYPE &y);
    //@}

#endif // EXPR_TEMPLATES


    // =======================================================================
    // SPECIAL OPERATIONS FOR THE 3-DIMENSIONAL CASE
    // =======================================================================
  
    //! Calc cross product of two vectors v = this x b
    void CrossProduct( const Vector<TYPE>& b, Vector<TYPE>& v );
    
    // =======================================================================
    // BOOLEAN OPERATORS
    // =======================================================================

    //! \name bool operators

    //@{
    
    //! Overloading of operation equal for Vector
    bool operator==( const Vector<TYPE> & ) const;
    
    //! Overloading of operation not equal for Vector
    bool operator!=( const Vector<TYPE> & ) const;
    
    //@}
    
    // =======================================================================
    // MISCELLANEOUS METHODS
    // =======================================================================

    //! \name Miscellaneous methods
    //@{

    //! Return part of Vector from index i to ii
    Vector<TYPE> Part( const UInt i, const UInt ii ) const;

    //! Insert a vector to this vector at position pos
    void InsertVector( const Vector<TYPE> & y, UInt pos = 0 );
    
    //! Constructs the unit vector of length n, which only non-zero
    //! entry is a 1 at the i-th position
    //!\param n (input) length of vector
    //!  \param i (input) component, which is 1
    //!  \f[ \left( \begin{array}{c} 0  \\ \cdots \\ 0 \\ 1 \\ 0 \\ \cdots \\ 0 
    //!  \end{array} \right) \f]
    static Vector<TYPE> Unit( const UInt n, const UInt i );

    //! Delete element from vector on position pos
    void Cut ( const UInt pos );

    //! Delete elements from position pos1 to pos2, on pos1, pos2 too
    void Cut ( const UInt pos1, const UInt pos2 );

    //! Sort vector in ascending order
    void Sort();
 
    //@}

  protected:

    //! Length of the vector
    UInt size_;

    //! Data of the vector
    TYPE*  data_;
  
    //! Capacity of the vector
    UInt capacity_;

    //! Flag signaling whether management of data array is done by this object
  
    //! This attribute is used to keep track on the fact whether the object
    //! is responsible for managing the memory of the data_ array, especially
    //! its deallocation.
    bool memBelongsToMe_;

    // =======================================================================
    // SERIALIZATION FUNCTIONS
    // =======================================================================
    // These functions allow us to write a vector directly
    // into an boost::archive, for saving on a disk or in a 
    // iostream object

    //! allow serialization class to access vector entries
    friend class boost::serialization::access;
    
    //! Saving internal state into a boost::archive
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;
    
    //! Reading internal state from a boost::archive
    template<class Archive>
    void load(Archive & ar, const unsigned int version);
    
    //! The following statement is needed for boost
    BOOST_SERIALIZATION_SPLIT_MEMBER()
      

  
  };
#ifdef DOXYGEN_DETAILED_DOC
  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================
  //! \class Vector
  //! 
  //! \purpose This class implements a general, templatized dense vector with 
  //! the following additional features:
  //! - If the macro USE_EXPR_TEMPLATES is defined, it utilizes an expression
  //! template library, which evaluates mathematical expression at compile 
  //! time.
  //! In this case, the operators(+,-,*,/,-=,+=,...) do not have to be defined
  //! in this class, but are evaluated by the expression template library.
  //! The used library is a modified version of the MET-library 
  //! (<a href="http://met.sourceforge.net">met.sourceforge.net</a>).
  //! - In order to be able to handle mixed Double-Complex valued mathematical
  //! expressions, the concept of Type Promotion / Traits is utilized (see also
  //! <a href="http://osl.iu.edu/~tveldhui/papers/techniques/"> Techniques for
  //! Scientific C++ </a>). This means, expressions like
  //! \verbatim
  //! Matrix<Double> realMat1;
  //! Matrix<Complex> complexMat1;
  //! Vector<Double> realVec1, realVec2;
  //! Vector<Double> complexVec1, complexVec2;
  //! Double realFactor = 1.0;
  //! Complex complexFactor = Complex(1.0, 1.0);
  //! Complex result;
  //!
  //! complexVec1 = complexVec2 * realFactor;
  //! complexVec1 = realVec1 * complexFactor;
  //!
  //! result = complexVec1 * realVec1; 
  //!
  //! complexVec1 = realMat1 * complexVec2;
  //! complexVec2 = complexMat1 * realVec1;
  //! \endverbatim
  //! can be written and the conversion is done automatically.
  //! \note - Multiple Double <-> Complex conversion in one statement are
  //!       not possible!
  //! 
  //! \note -If expression templates are used, statements like
  //! \verbatim
  //! Vector<Double> vec = vec1 * 5.0;
  //! \endverbatim 
  //! have to be replaced by
  //! \verbatim
  //! Vector<Double> vec;
  //! vec = vec * 5.0;
  //! \endverbatim 
  //! 
  //! \collab The Vector class can be used together with the templatized Matrix
  //! class.
  //! 
  //! \implement This class uses the concept of type promotion / traits and
  //! can additionally utilize expression templates.
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

    
  // ******************************************************
  // * Additional functions related with handling vectors *
  // ******************************************************
    
  //! Function for swap  number a and b (Integer, Double)
  template <class S> void Swap(S &, S &);
    
  //! Function for swap Vector<T>;
  template<class TYPE>
  void Swap(Vector<TYPE> & a, Vector<TYPE> & b);
    
  //! Overloading << for class vector
  template<class TYPE>  std::ostream& operator << ( std::ostream & , 
                                                    const Vector<TYPE> &);
    

  // ******************************************************
  //  INLINE MEMBER DEFINITIONS 
  // ******************************************************
  template <class TYPE> 
  inline TYPE & Vector<TYPE>::operator[] ( const UInt i ) {
#ifdef CHECK_INDEX
    if ( i >= size_ ) {
      EXCEPTION( "Vector: Invalid access, index " << i
                 << " not in range [0," << size_ - 1 << "]" );
    }
#endif
    return  data_[i];
  }

  template <class TYPE> 
  inline TYPE Vector<TYPE>::operator[] ( const UInt i ) const {
#ifdef CHECK_INDEX
    if ( i >= size_ ) {
      EXCEPTION( "Vector: Invalid access, index " << i
                 << " not in range [0," << size_ - 1 << "]" );
    }
#endif
    return  data_[i];
  }

  template <class TYPE> 
  inline TYPE * Vector<TYPE>::GetPointer() const 
  {
    if (!data_)
      EXCEPTION( "Vector: undefined Vector" );
    return data_;
  }

  template <class TYPE> 
  void Vector<TYPE>::GetPointer(TYPE* &ptr_out) const 
  {
    if (!data_)
      EXCEPTION( "Vector: undefined Vector" );

    ptr_out = data_;
  }
  

  // ************************************************************
  //  INLINE MEMBER DEFINITIONS FOR NON-TEMPLATE EXPRESSION CASE
  // ************************************************************
  template<class TYPE> template<class TYPE2>
  PROMOTE(TYPE,TYPE2) Vector<TYPE>::
  operator* (const Vector<TYPE2> &x) const
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.GetSize() == 0))
      EXCEPTION( "Vector: undefined Vector in operator *(vector)" );
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.GetSize())
      EXCEPTION( "Vector: incompatible dimension for operator *(vector)" );
#endif

    PROMOTE(TYPE,TYPE2) ret;
 
    ret = data_[0] * x[0];
    for (UInt i = 1; i < size_; i++)
      ret += data_[i] * x[i];
  
    return ret;
  }

#ifndef EXPR_TEMPLATES
  template<class TYPE> template<class TYPE2>
  Vector<PROMOTE(TYPE,TYPE2)> Vector<TYPE>::
  operator+(const Vector<TYPE2> &x) const
  {       

#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.GetSize() == 0))
      EXCEPTION( "Vector: undefined Vector in operator +(vector)" );
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.GetSize())
      EXCEPTION( "Vector: incompatible dimension for operator +(vector)" );
#endif
  
    Vector<PROMOTE(TYPE,TYPE2)> ret(size_);
  
    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] + x[i];
  
    return ret;
  }

  template<class TYPE> template<class TYPE2>
  Vector<PROMOTE(TYPE,TYPE2)> Vector<TYPE>::
  operator-(const Vector<TYPE2> &x) const
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.GetSize() == 0))
      EXCEPTION( "Vector: undefined Vector in operator -(vector)" );
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.GetSize())
      EXCEPTION( "Vector: incompatible dimension for operator -(vector)" );
#endif
    Vector<PROMOTE(TYPE,TYPE2)> ret(size_);

    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] - x[i];
  
    return ret;
  }

  template<class TYPE> template<class TYPE2>
  Vector<PROMOTE(TYPE,TYPE2)> Vector<TYPE>::
  operator* (const TYPE2 &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION( "Vector: undefined Vector in operator *(number)" );
#endif
  
    Vector<PROMOTE(TYPE,TYPE2)> ret(size_);
  
    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] * x;
  
    return ret;
  }

  template<class TYPE> template<class TYPE2>
  Vector<PROMOTE(TYPE,TYPE2)> Vector<TYPE>::
  operator/(const TYPE2 &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION( "Vector: undefined Vector in operator /(number)" );
#endif
  
    Vector<PROMOTE(TYPE,TYPE2)> ret(size_);
  
    for (UInt i = 0; i < size_; i++)
      ret[i] = data_[i] / x;
  
    return ret;
  }
#endif // EXPR_TEMPLATE


  
  template <class TYPE> 
  void Vector<TYPE>::CrossProduct( const Vector<TYPE>& b, 
                                   Vector<TYPE>& v ) {
    if( size_ != 3 || b.size_ != 3 )
      EXCEPTION("CrossProduct can only be calculated for vector of size 3!");

    v.Resize(3);
    v[0] = data_[1] * b.data_[2] - data_[2] * b.data_[1];
    v[1] = data_[2] * b.data_[0] - data_[0] * b.data_[2];
    v[2] = data_[0] * b.data_[1] - data_[1] * b.data_[0];
  }

  template<class TYPE> template< class Archive>
  void Vector<TYPE>::save(Archive & ar, const unsigned int version) const {

    // invoke serialization of the base class 
    ar & boost::serialization::base_object<CFSVector>(*this);

    ar & size_;
    for( UInt i = 0; i < size_; i++ ) 
      ar & data_[i];
  }
  
  template<class TYPE> template <class Archive>
  void Vector<TYPE>::load(Archive & ar, const unsigned int version) {

    // invoke serialization of the base class 
    ar & boost::serialization::base_object<CFSVector>(*this);

    if( data_ != NULL ) {
      delete[] data_;
    }
    
    ar & size_;
    data_ = new TYPE[size_];
    for( UInt i = 0; i < size_; i++ ) {
      ar & data_[i];
    }
  }

} // end of namespace

#endif
