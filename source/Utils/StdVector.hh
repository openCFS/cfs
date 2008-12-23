// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STDVECTOR_2004
#define FILE_STDVECTOR_2004

#include <boost/serialization/split_member.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <vector>
#include "General/exception.hh"
#include <def_build_type_options.hh>

namespace CoupledField {


  // forward class declaration
  template <class TYPE> class StdVectorListInitializer;

  //! General template class for storing elements of same type
  //! in sequential order. In contrast to the CFS-Vector, this class
  //! requires no arithmetic functionality for its elements.
  //! This class replaces the use of std::vector<>-class and therefore
  //! provides a similar interface
  template<class TYPE>
  class StdVector{
  public:

    // =======================================================================
    //  STL-COMPATIBLE ITERATOR DEFINITIONS
    // =======================================================================
    // The implementation of the iterator interface utilizes boost's iterator-
    // facade concept.

    //! Define iterator class
    class iterator 
      :  public boost::iterator_facade
      < iterator, 
        TYPE,
        boost::random_access_traversal_tag
        > {
    public: 
      
      //! default constructor
      iterator() : vec_(NULL), pos_(0) {}
      
    private:

      friend class boost::iterator_core_access;
      friend class StdVector<TYPE>;

      //! iterator with pointer to vector
      iterator( StdVector<TYPE>* p, unsigned int pos = 0 ) 
        : vec_( p ), pos_( pos ) {}
    
      //! dereferencing
      TYPE& dereference() const { 
        return (*vec_)[pos_]; 
      }
      
      //! check for equality
      bool equal( iterator const & other ) const {
        return ( this->vec_ == other.vec_ &&
                 this->pos_ == other.pos_ );
      }
      
      //! increment position
      void increment() { 
        pos_++;
      }
      
      //! decrement position
      void decrement() { pos_--; }
      
      //! advance position by n
      void advance( ptrdiff_t n) { 
        pos_ += (unsigned int) n; }
      
      //! measure distance
      ptrdiff_t distance_to( iterator const & other ) const {
        return other.pos_ - this->pos_;
      }
      
      // references to StdVector
      StdVector<TYPE> * vec_;
      size_t pos_;
    };
    
    //! Define CONST iterator class
    class const_iterator 
    :  public boost::iterator_facade
    < const_iterator, 
    TYPE const,
    boost::random_access_traversal_tag
    > {
    public: 

      //! default constructor
      const_iterator() : vec_(NULL), pos_(0) {}

    private:

      friend class boost::iterator_core_access;
      friend class StdVector<TYPE>;

      //! iterator with pointer to vector
      const_iterator( const StdVector<TYPE>* p, unsigned int pos = 0 ) 
      : vec_( p ), pos_( pos ) {}

      //! dereferencing
      const TYPE& dereference() const { 
        return (*vec_)[pos_]; 
      }

      //! check for equality
      bool equal( const_iterator const & other ) const {
        return ( this->vec_ == other.vec_ &&
            this->pos_ == other.pos_ );
      }

      //! increment position
      void increment() { 
        pos_++;
      }

      //! decrement position
      void decrement() { pos_--; }

      //! advance position by n
      void advance( ptrdiff_t n) { 
        pos_ += (unsigned int) n; }

      //! measure distance
      ptrdiff_t distance_to( const_iterator const & other ) const {
        return other.pos_ - this->pos_;
      }

      // references to StdVector
      const StdVector<TYPE> * vec_;
      size_t pos_;
    };
  
    //! Return iterator pointing to first element
    iterator Begin() {
        return iterator(this, 0);
    }
       
    //! Return iterator pointing to first element
    const_iterator Begin() const {
      return const_iterator(this, 0);
    }
    
    //! Return iterator pointing to last element
    iterator End() {
        return iterator(this, size_);
    }
    
    //! Return iterator pointing to last element
    const_iterator End() const {
        return const_iterator( this, size_ );
    }
    
    // =======================================================================
    //  VECTOR METHODS
    // =======================================================================

    //! Constructor
    StdVector();

    //! Constructor with inital size.
    //! All entries are filled with zeroes
    explicit StdVector(unsigned int size);

    //! Copy constructor
    StdVector(const StdVector<TYPE> & vec);

    //! Copy constructor with std::vector
    StdVector(const std::vector<TYPE> & vec);

    //! Destructor
    ~StdVector();  

    //! Deletes all entries in the vector
    //! and frees its memory
    void Clear();

    //! Initalizes the vector with a given entry
    /*!
      \param entry (input) Entry vector gets inialized with
    */
    //! \note this method does not change the size of the vector!
    void Init(const TYPE entry = TYPE());

    //! True, if vector is empty
    inline bool IsEmpty() const {return (size_? false : true);}
  
    //! Returns capacity of the vector
    inline unsigned int Capacity() const {return capacity_;}

    //! 
    void Reserve(unsigned int size);

    //! Get the length of the vector
    inline unsigned int GetSize() const {return size_;}

    /** Set the length of the vector.
     * @note The vector is NOT initialized, use the Resize with init parameter for this purpose */
    void Resize(const unsigned int size);

    /** Set the length of the vector and initilize
     * @note init Init() is called with this value */
    void Resize(const unsigned int size, TYPE entry);
    
    //! Overloading of operation =
    StdVector     &operator=      (const StdVector &);

    //! Overloading for operator=

    //! This method is needed to intialize a StdVector like this
    //! StdVector<Integer> A;<br>
    //! A = 1,2,5,10
    inline StdVectorListInitializer<TYPE> operator=(const TYPE x);

    //! Build vector from std::vector
    StdVector & operator= (const std::vector<TYPE> & vec);
  
    //! Returns las entry of vector
    inline TYPE & Last();

    //! Returns last entry of vector (read-only)
    inline TYPE Last() const;
 

    //! General access operator
    inline TYPE &operator[] (const unsigned int i);

    //! Read-Only access operator
    inline const TYPE & operator[] (const unsigned int i) const;

    //! Return pointer p to array 
    inline TYPE*  GetPointer() const;

    /** Imports data from extern, adjusts internal size and capacity.
     * Any existing data is overwritten. */
    void Import(const TYPE* source, unsigned int size);
    
    //! Add element of the same type at position pos
    void Insert(const TYPE & y, unsigned int pos);

    //! Add element of the same type at position pos
    //! numCopies times
    void Insert(const TYPE & y, unsigned int pos, unsigned int numCopies);

    //! Add element of the same type at the end of the vector
    void Push_back(const TYPE & y = TYPE() );

    //! Delete element from vector on position pos
    void Erase (const unsigned int pos);

    //! Delete elements from position pos1 to pos2, on pos1, pos2 too
    void Erase (const unsigned int pos1, const unsigned int pos2);


    //! Return the position number of element x in the vector

    //! Finds the element x in the vector and returns the 
    //! position. If no element was found, it returns -1.
    int Find(const TYPE &x) const;

    /** Checks if an element exists */ 
    bool Contains(const TYPE &x) const { 
      return Find(x) != -1; 
    } 
    
    //! Overloading of operation equal for Vector
    bool operator== (const StdVector &) const;

    //! Overloading of operation equal for data array
    bool operator== (const TYPE* other) const;
    
    
    //! Overloading of operation not equal for Vector
    bool operator!= (const StdVector &) const;
  

    // ******************************************************
    // MISCELANEOUS FUNCTIONS
    // ******************************************************
  
    //! Sort of vector: v - vec.p, n - vec.size
    template <class S> void Sort(S* v, unsigned int n);
  
    //! Swap 2 elements in vector Ex Swap(v[i],v[j])
    template<class T2> void Swap(T2& a, T2 & b);

    //! Return vector as separated string
    std::string Serialize( char separator = ',') const;

   
     /** Lists the content comma seperated.
      * @param level 0=all data, 1 is summary, the higher, the less 
      * @param stride on level=0 every element(1), every second (2), ...*/ 
     std::string ToString(int level=0, int stride=1) const;
     
     /** List the content or summary of an external source 
      * @see ToString(int)*/
     static std::string ToString(int size, const TYPE* data, int level=0, int stride=1);

     /** converts the content to a string vector */
     void ToString(StdVector<std::string>& out) const;

     /** Reas the content from a string list */
     void Parse(const StdVector<std::string>& in);
     
  protected:

    //! Length of the vector
    unsigned int size_;
    
    //! Data of the vector
    TYPE* data_;
  
    //! Capacity of the vector
    unsigned int capacity_;

    // ******************************************************
    // SERIALIZATION FUNCTIONS
    // ******************************************************
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

  // ******************************************************
  // HELPER CLASS FOR INITALIZING StdVector
  // (ref. 'Techniques for Scientific C++' 
  // by Todd Veldhuizen, page. 43ff
  // ******************************************************
  template <class TYPE>
  class StdVectorListInitializer
  {
  public:

    //! Constructor
    StdVectorListInitializer(StdVector<TYPE> * vec)
      :vec_(vec) {};
    
    //! Overloading of comma operator
    StdVectorListInitializer<TYPE> operator, (const TYPE x)
    {
      vec_->Push_back(x);
      return StdVectorListInitializer<TYPE>(vec_);
    }

  private:
    //! pointer to vector
    StdVector<TYPE> * vec_;
  };

  // ******************************************************
  // INLINE FUNCTION DECLARATION
  // ******************************************************
  
  //! Element can be referred to as v[i]
  template<class TYPE>
  TYPE & StdVector<TYPE>::operator[] (const unsigned int i)
  {     
    #ifdef CHECK_INDEX
    if (i >= size_){
      EXCEPTION( "Invalid access to element "
                 << i+1
                 << " \n Length of vector: " << size_ );
    }
    #endif
    return  data_[i];
  }

  //! Element can be referred to as v[i]
  template<class TYPE>
  const TYPE & StdVector<TYPE>::operator[] (const unsigned int i) const
  {     
#ifdef CHECK_INDEX
     if (i >= size_){
       EXCEPTION( "Invalid access to element "
                  << i+1
                  << " \n Length of vector: " << size_ );
    }
#endif
    return  data_[i];
  }
  

  //! 
  template<class TYPE>
  StdVectorListInitializer<TYPE> StdVector<TYPE>::operator= (const TYPE x)
  {
    Clear();
    Push_back(x);
    return StdVectorListInitializer<TYPE>(this);
  }

  //! Return pointer p to array 
  template<class TYPE>
  TYPE* StdVector<TYPE>::GetPointer() const
  {
    if (!data_)
      EXCEPTION("Vector: undefined Vector" );
    return data_;
  }

  template<class TYPE>
  TYPE &StdVector<TYPE>::Last() {
#ifdef CHECK_INITIALIZED
    if ( size_ == 0)
      EXCEPTION("undefined Vector");
#endif
      return data_[size_-1] ;
  }
  
  template<class TYPE>
  TYPE StdVector<TYPE>::Last() const {
#ifdef CHECK_INITIALIZED
    if ( size_ == 0)
      EXCEPTION("Vector: undefined Vector");
#endif
    return data_[size_-1] ;
  }

  // ******************************************************
  // * Additional functions related with handling vectors *
  // ******************************************************
  
  //! Overloading << for class vector
  template<class TYPE>  std::ostream& operator << ( std::ostream & , const StdVector<TYPE> &);

} // end of namespace

// Inclusion of Implementation
#include "StdVector.cc"

#endif
