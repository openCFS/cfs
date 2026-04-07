#ifndef FILE_STDVECTOR_2004
#define FILE_STDVECTOR_2004

#include <array>
#include <boost/iterator/iterator_facade.hpp>
#include <vector>
#include "General/Exception.hh"
#include "General/EnvironmentTypes.hh"
#include <def_build_type_options.hh>
#include <limits>

namespace CoupledField {


  // forward class declaration
  template <class TYPE> class StdVectorListInitializer;

  //! General template class for storing elements of same type
  //! in sequential order. In contrast to the CFS-Vector, this class
  //! requires no arithmetic functionality for its elements.
  //! This class replaces the use of std::vector<>-class and therefore
  //! provides a similar interface. 
  //! One advantage is that we have range checks in debug mode
  //! Has some constructors and Replace() to wrap external data.
  //! Resizing with external data is allowed as long as we stay within the capacity
  template<class TYPE>
  class StdVector{
  public:

    //! public typedef for STL compatibility
    typedef TYPE value_type;
    typedef TYPE& reference;
    typedef const TYPE& const_reference;
    typedef ptrdiff_t difference_type;
    typedef unsigned int size_type;
    
    // =======================================================================
    //  STL-COMPATIBLE ITERATOR DEFINITIONS
    // =======================================================================
    // The implementation of the iterator interface utilizes boost's iterator-
    // facade concept.

    //! Define iterator class
    class iterator : public boost::iterator_facade<iterator, TYPE, boost::random_access_traversal_tag> 
    {
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
        if( n > 0 )
          pos_ += std::abs(n);
        else
          pos_ -= std::abs(n);
      }
      
      //! measure distance
      ptrdiff_t distance_to( iterator const & other ) const {
        return other.pos_ - this->pos_;
      }
      
      // references to StdVector
      StdVector<TYPE> * vec_;
      size_t pos_;
    };
    
    //! Define CONST iterator class
    class const_iterator : public boost::iterator_facade<const_iterator, TYPE const, boost::random_access_traversal_tag>
    {
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
        if( n > 0 )
          pos_ += std::abs(n);
        else
          pos_ -= std::abs(n);
      }
      
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
    
    //! Return iterator pointing to first element
    //! otherwise we can't use C++11 auto
    iterator begin() {
        return iterator(this, 0);
    }

    //! Return iterator pointing to first element
    //! otherwise we can't use C++11 auto
    const_iterator begin() const {
      return const_iterator(this, 0);
    }

    //! Return iterator pointing to (undefined) element after last element
    iterator End() {
        return iterator(this, size_);
    }
    
    //! Return iterator pointing to (undefined) element after last element
    const_iterator End() const {
        return const_iterator( this, size_ );
    }
    

    //! Return iterator pointing to (undefined) element after last element
    //! otherwise we can't use C++11 auto
    iterator end() {
        return iterator(this, size_);
    }

    //! Return iterator pointing to (undefined) element after last element
    //! otherwise we can't use C++11 auto
    const_iterator end() const {
        return const_iterator( this, size_ );
    }

    // =======================================================================
    //  VECTOR METHODS
    // =======================================================================

    //! Constructor
    StdVector() { }; // default is fine

    //! Constructor with initial size.
    //! All entries are filled with zeroes
    explicit StdVector(size_type size) {
      size_ = size;
      capacity_ = size;
      data_ = new TYPE[size];
      Init();
    }

    //! Copy constructor - we copy to own memory. 
    StdVector(const StdVector<TYPE>& vec) {
       size_ = vec.size_;
       capacity_ = vec.size_;
       data_ = new TYPE [vec.size_];
       std::copy(vec.data_, vec.data_+size_, data_);
    }

    /** will copy data to own memory */
    StdVector(const std::vector<TYPE>& vec) : StdVector(vec.data(), vec.size()) {}
    
    /** UseStdVector as wrapper to existing data, e.g. std::array<double, 4>.
     * Avoids dynamic memory allocation on head. Saves > 10 times for small vectors. 
     * @see heap_vs_stack() in testbed.cc
     * We do not copy nor own the data. Resizing() is allowed (assert in debug)
     * This is a wrapper, it acts directly with the external data
     * @see Replace() */
    StdVector(TYPE* data, size_type size, ExternalDataMode mode = WRAP) {
      Replace(data, size, mode);
    }

    /** const data is copied by default. It must not be wrapped.
     * If you want the other way, const_cast yourself */
    StdVector(const TYPE* data, size_type size) {
      Replace(const_cast<TYPE*>(data), size, COPY); 
    }

    /** Use this if you want to provide hast stack data to be wrapped. */
    template<std::size_t N>
    explicit StdVector(std::array<TYPE, N>& arr) {
      size_ = N;
      capacity_ = N;
      data_ = arr.data(); 
      wrapped_ = true;
    }

    //! STL-compatible version
    template <class InputIterator>
    StdVector(InputIterator first, InputIterator last);

    //! Destructor
    ~StdVector() 
    {
      if(!wrapped_)
        delete[] data_;
    }

    /** helper for easy usage with a stack array as buffer to be wrapped by a StdVector. Use like
     * auto [myvec, buffer] = StdVector<double>::CreateWithStackData<4>(size);
     * myvec will be a vector, buffer is just there to be the underlaying data. */ 
    template<std::size_t N>
    static std::pair<StdVector<TYPE>, std::array<TYPE, N>> CreateWithStackData() {
      std::array<TYPE, N> arr; // compiler magic will allocate the objects on the callers stack
      return std::make_pair(StdVector<TYPE>(arr), arr);
    }
    
    //! Clear the vector

    //! This method clears the vector, i.e. it sets its size
    //! to zero. If keepCapacity is set to false
    //! also the memory is freed. Otherwise, the old data
    //! pointer is kept, so successive Resize-operations are
    //! faster.
    //! By default the capacity is kept for performance reasons. When
    //! the space is not needed any more often the vector itself will be destructed?!
    inline void Clear(bool keepCapacity = true);

    //! Initalizes the vector with a given entry
    /*!
      \param entry (input) Entry vector gets initialized with
    */
    //! \note this method does not change the size of the vector!
    void Init(const TYPE entry = TYPE()) {
      std::fill(data_, data_+size_, entry);
    }


    //! True, if vector is empty
    inline bool IsEmpty() const {return (size_? false : true);}
  
    //! True, if vector is empty
    inline bool empty() const {return (size_? false : true);}
  
    /** Returns capacity of the vector
      * @see Reserve() */
    inline size_type GetCapacity() const {return capacity_;}

    /** Is there space for another Push_back() ? */
    bool HasSpace() const { return size_ < capacity_; }

    /** Increases the capacity.
     * If the new capacity is smaller than the current one nothing happens.
     * Otherwise the old data is copied up to size and the rest is uninitialized!
     * @see GetCapacity() */
    void Reserve(size_type capacity);

    //! Delete non-used reserved memory
    
    //! This methods ensures, that the size and capacity of the vector have the
    //! size, i.e. that no memory is reserved any more for future Push_back()
    //! operations. This method is especially useful, if a vector gets 
    //! successively filled by Push_back-operations and at the end the non-used
    //! capacity is not needed anymore.
    void Trim();
    
    //! Get the length of the vector
    inline size_type GetSize() const {return size_;}

    //! Get the length of the vector
    inline size_type size() const {return size_;}

    /** Get maximal length of vector.
     * On Windows max() could be a macro, so call #undef max in some header or #define NOMINMAX #include <windows.h> */
    inline size_type max_size() const {return std::numeric_limits<size_type>::max();}

    /** Set the length of the vector but keep the capacity.
     * When we shrink, capacity will remain. If we grow beyond capacity we copy 
     * When we are wrapped it is allowed as long as we are within the capacity
     * @param size if smaller capacity only the internal size parameter is adjusted.
     *        If larger than the current capacity the old data is copied!
     * @note Additional data is NOT initialized, use the other constructor Resize with init parameter sets ALL data */
    inline void Resize(size_type size)
    {
      if(size > capacity_)
      {
        assert(!wrapped_);
        TYPE* help = new TYPE[size];
        for (unsigned int i=0; i<size_; i++)
          help[i] = data_[i];
        delete[] data_;
        data_ = help;
        capacity_ = size;
      }
      size_ = size;
    }

    /** Resize by keeping the capacity. In case we grow beyond capacity, the data is not copied. */ 
    void ResizeNoCopy(size_type size)
    {
      if(size > capacity_)
      {
        assert(!wrapped_);
        delete[] data_;
        data_ = new TYPE[size];
        capacity_ = size;
      }
      size_ = size;
    }

    /** Set the length of the vector and initialize
     * @note Init() is called with this value */
    inline void Resize(size_type size, TYPE entry) {
      ResizeNoCopy(size);
      Init(entry);
    }

    /** extract part of the vector by a pattern
     * @param start this is index of the first item to be returned
     * @param stride the returned values are for start, start+stride, start+2*stride, ... */
    StdVector Slice(size_type start, size_type stride);

    //! Overloading of operation =
    StdVector& operator= (const StdVector &);

    //! Overloading for operator=

    //! This method is needed to initialize a StdVector like this
    //! StdVector<Integer> A;<br>
    //! A = 1,2,5,10
    inline StdVectorListInitializer<TYPE> operator=(const TYPE x);

    //! Build vector from std::vector
    StdVector& operator= (const std::vector<TYPE> & vec);
  
    //! Returns the last entry of vector
    inline TYPE& Last();
    inline TYPE Last() const;
 
    //! Returns the first entry of vector
    inline TYPE& First();
    inline TYPE First() const;

    /** convenience when having a pointer. In debug the operator checks bounds */
    inline TYPE& At(size_type i) { return data_[i]; }
    inline const TYPE& At(size_type i) const { return data_[i]; }

    //! General access operator
    inline TYPE& operator[] (const unsigned int i);

    //! Read-Only access operator
    inline const TYPE& operator[] (const unsigned int i) const;

    //! Return pointer p to array 
    inline TYPE*  GetPointer() const;

    /** Imports data from external, adjusts internal size and capacity. The data is copied.
     * Any existing data is overwritten.
     * @see Replace(COPY) what is actually done */
    void Import(const TYPE* source, size_type size) { 
      Replace(const_cast<TYPE*>(source), size, COPY); 
    }

    /** Import data from an external stl container. Resizes this vector accordingly.
     * Not all containers have a size, e.g. boost::tokenizer. If possible provide it, otherwise Push_back might have to
     * copy the data multiple time.
       @param begin assign your begin iterator
     * @param end your end iterator
     * @param size if known for you container. Sets capacity to the value */
    // https://stackoverflow.com/questions/7728478/c-template-class-function-with-arbitrary-container-type-how-to-define-it
    template <typename Iter>
    void Import(Iter it, Iter end, int size = -1) {
      Resize(0);
      if(size > 0)
        Reserve(size);

      for(; it!=end; ++it)
         Push_back(*it);
    }

    /** Make a wrapper or copy external data to internal data (will be resized).
     * Own memory will be deleted first when we wrap. Otherwise we leave the object deletion to the assignment operator.
     * Usefull e.g. to assign a stack buffer std::array<T,N> 
     * @see Import() */
    inline void Replace(TYPE* data, size_type length, ExternalDataMode mode = WRAP)
    {
      assert(mode == WRAP || mode == COPY);
      if(mode == WRAP) {
        if(!wrapped_)
          delete[] data_;
        data_ = data;
        size_ = length;
        capacity_ = length;
        wrapped_ = true;
      } else {
        assert(data != nullptr && length > 0);
        if(wrapped_) { // when we are wrapped but shall copy, we allocate own memory and are not wrapped any more
          data_ = nullptr;
          size_ = 0;
          capacity_ = 0;
          wrapped_ = false; // reset from wrapped before calling ResizeNoCopy()
        }
        ResizeNoCopy(length); 
        std::copy_n(data, length, data_);
        assert(wrapped_ == false);
      }
    } 

    /** Add element of the same type at the end of the vector.
     * If there is not enough capacity (GetCapacity()) the data is copied to a new data of doubled size.
     * Note that this invalidates all pointers to data entries! If you have such an issue, make sure you Reserve() sufficient
     * space before Push_back() the first element.
     * @param no_expand throws an exception if the capacity is too small and data expanding would be necessary
     * @see HasSpace() to check if Push_back() does not need expansion */
    inline TYPE& Push_back(const TYPE & y = TYPE(), bool no_expand = false) 
    {
      if ( size_ < capacity_ ) {
        data_[size_++] = y;
      } else {
        if(no_expand)
          EXCEPTION("Capacity " << capacity_ << " for Push_back() too small but data expansion not allowed");
        Push_back_expand(y);
      }
      return data_[size_-1];
    }

    /** Lower-case version for STL compatibility */
    inline void push_back(const TYPE &y = TYPE()) {
      Push_back(y);
    }

    /** convenience function adding two elements */
    inline void Push_back(const TYPE& a, const TYPE& b) {
      Push_back(a);
      Push_back(b);
    }

    //! Delete element from vector on position pos
    void Erase(size_type pos);

    //! Delete elements from position pos1 to pos2, on pos1, pos2 too
    void Erase(size_type pos1, size_type pos2);

    /** Insert the element at the given position, anything form that position on is moved one back.
     * @param pos the last valid pos is size which corresponds to an append */
    void Insert(size_type pos, const TYPE& dat);
    
    //! Append 2nd vector at the end
    void Append(const StdVector<TYPE>& vec );

    /** Return the position number of element x in the vector.
     * Finds the element x in the vector and returns the
     * position. If no element was found, it returns -1.
     * @param quiet if not, throw exception in failed case */
    int Find(const TYPE &x, bool quiet = true) const;

    /** Return the position number of element x in the vector. Faster with a guess
     * Takes a guess on where to start linear searching.
     * Consider if you need a thread local storage for the guess!
     * @param guess when found, guess will set to index+1 or 0 if at end */
    int Find(const TYPE &x, unsigned int& guess, bool quiet = true) const;

    //! Finds all elements x in the vector and returns the
    //! position. If no element was found, the length of
    //! the returned vector is zero.
    StdVector<unsigned int> FindAll(const TYPE &x) const;


    /** Checks if an element exists */ 
    bool Contains(const TYPE &x) const { 
      return Find(x) != -1; 
    } 
    
    /** Checks if there are non-unique elements */
    bool IsUnique() const;

    //! Overloading of operation equal for Vector
    bool operator== (const StdVector &) const;

    //! Overloading of operation equal for data array
    bool operator== (const TYPE* other) const;
    
    
    //! Overloading of operation not equal for Vector
    bool operator!= (const StdVector &) const;
  

    // ******************************************************
    // MISCELANEOUS FUNCTIONS
    // ******************************************************

    /** Swap to entries */
    void Swap(size_type idx1, size_type idx2);

    //! Return vector as separated string
    std::string Serialize( char separator = ',') const;

     /** print content for logging or nice output.
      * See unittests.cc Matrix_ToString for test cases
      * TS_PLAIN is space separated
      * TS_MATLAB is with brackets, no comma between elements
      * TS_PYTHON is with brackets, commas between elements
      * TS_INFO gives summary
      * TS_NONZEROS invalid option
      * @param sep separator string, when empty meaningful default is used
      * @param in_window only for window. what requires the window to be set */
     std::string ToString(ToStringFormat format = TS_PLAIN, const std::string& sep="", bool in_window = false) const;

     /** comma separated output for window only  */
     std::string ToString(bool in_window) const {
       return ToString(TS_PLAIN, ", ", in_window);
     }

     /** List the content or summary of an external source 
      * @see ToString() */
     static std::string ToString(int size, const TYPE* data, ToStringFormat format = TS_PLAIN, const std::string& sep="");

     /** converts the content to a string vector */
     void ToString(StdVector<std::string>& out) const;

     /** Reads the content from a string list */
     void Parse(const StdVector<std::string>& in);

     /** This is a little helper to define a range within the data.
      * It is just a rucksack of information with no further implication.
      * One could add verification of the ranges but this is also ensured by
      * the vector itself. */
     class Window
     {
     public:

       /** default construcor */
       Window();

       /** Sets automatically acive */
       void Set(unsigned int start, unsigned int size);

       /** Scale to the whole vector dimension */
       void Set(StdVector& vec);

       /** has the window valid properties */
      bool Initialized() const { return active_; }

       /** index of the first element */
       unsigned int GetStart() const;

       /** size of the window */
       unsigned int GetSize() const;

     private:

       unsigned int start_;
       unsigned int size_;
       bool active_;
     };

     /** Out rucksack data */
     Window window;
     
     /** Verifies if we are in the window. Error also if there is no window initialized (debug mode only) */
     bool InWindow(unsigned int index);

  protected:

    //! Internal push back command for resizing
    void Push_back_expand(const TYPE & y);

    //! Length of the vector
    size_type size_ = 0;
    
    //! Capacity of the vector
    size_type capacity_ = 0;

    //! Data of the vector
    TYPE* data_ = nullptr;

    /** there is the case of being a wrapper for external data, as in Vector.hh,
     * however we may Resize() within the capacity - so we need a flag */
    bool wrapped_ = false;
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
    StdVector<TYPE>* vec_;
  };

  // ******************************************************
  // INLINE FUNCTION DECLARATION
  // ******************************************************
  
  //! Element can be referred to as v[i]
  template<class TYPE>
  TYPE& StdVector<TYPE>::operator[] (size_type i)
  {     
    #ifdef CHECK_INDEX
    if (i >= size_){
      EXCEPTION( "Invalid access to element " << i+1 << " Length of vector: " << size_ );
    }
    #endif
    return  data_[i];
  }

  //! Element can be referred to as v[i]
  template<class TYPE>
  const TYPE& StdVector<TYPE>::operator[] (size_type i) const
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
    Clear(false); // don't keep old capacity
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
  TYPE& StdVector<TYPE>::Last() {
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

  template<class TYPE>
  TYPE& StdVector<TYPE>::First() {
#ifdef CHECK_INITIALIZED
    if ( size_ == 0)
      EXCEPTION("undefined Vector");
#endif
      return data_[0] ;
  }

  template<class TYPE>
  TYPE StdVector<TYPE>::First() const {
#ifdef CHECK_INITIALIZED
    if ( size_ == 0)
      EXCEPTION("Vector: undefined Vector");
#endif
    return data_[0] ;
  }


  // ******************************************************
  // * Additional functions related with handling vectors *
  // ******************************************************
  
  //! Overloading << for class vector
  template<class TYPE>  std::ostream& operator << ( std::ostream & , const StdVector<TYPE> &);

} // end of namespace

// Inclusion of Implementation
// We do it for StdVector since it may be a container for
// many different types unlike Vector which only holds doubles,
// UInts and a few more.
#include "StdVector.cc"

#endif
