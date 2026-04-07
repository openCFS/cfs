#ifndef FILE_STDVECTOR_IMPLEMENTATION_2004
#define FILE_STDVECTOR_IMPLEMENTATION_2004

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <limits>

namespace CoupledField {

  template<class TYPE>
  template<class InputIterator>
  StdVector<TYPE>::StdVector(InputIterator first, InputIterator last) {
    size_ = 0;
    capacity_ = size_;
    data_ = NULL;
    for (; first != last; ++first)
      Push_back(*first);
  }

  template<class TYPE>
  void StdVector<TYPE>::Clear(bool keepCapacity )
  {
    size_ = 0;
    if( keepCapacity == false ) {
      capacity_ = 0;
      if(!wrapped_)
        delete[] data_;
      data_ = nullptr;
      wrapped_ = false; // default case
    }
  }

  template<class TYPE>
  void StdVector<TYPE>::Reserve(size_type capacity)
  {
    unsigned int old_size = size_;

    ResizeNoCopy(capacity);

    size_ = old_size;
  }
  
  
  template<class TYPE>
  void StdVector<TYPE>::Trim()
  {
    // only perform trimming, if size_ and capacity_ differ
    if( capacity_ == size_ )
      return;

    assert(!wrapped_);  
    TYPE* help = new TYPE[size_];
    for (unsigned int i=0; i<size_; i++)
      help[i] = data_[i];
    delete[] data_;
    data_ = help;
    capacity_ = size_;
  }

  template<class TYPE>
  StdVector<TYPE> StdVector<TYPE>::Slice(size_type start, size_type stride)
  {
    assert(stride > 0);
    size_type size = start >= size ? 0 : (size_type) ((size_ - start) / stride) + 1;
    StdVector<TYPE> res(size);
    for(size_type i = 0; i < size; i++)
      res[i] = data_[start + i * stride];
    return res;
  }
  
  // *************
  //   operator=
  // *************
  template<class TYPE>
  StdVector<TYPE>& StdVector<TYPE>::operator= ( const StdVector &vec ) 
  {
    if ( this == &vec ) 
      return *this;
    
    // If there is not enough space to copy the entries
    // perform a re-allocation
    if(capacity_ < vec.size_)
    {
      // Delete old buffer
      assert(!wrapped_);
      delete[] data_;

      // Allocate new buffer
      capacity_ = vec.size_;
      data_     = new TYPE[capacity_];
    }

    // Copy entries
    size_ = vec.size_;
    if(data_)
      std::copy(vec.data_, vec.data_+size_, data_);

    return *this;
  }


  // *************
  //   operator=
  // *************
  template<class TYPE>
  StdVector<TYPE>& StdVector<TYPE>::operator=(const std::vector<TYPE> &vec)
  {
    if(capacity_ < vec.size())
    {
      assert(!wrapped_);
      delete[] data_;
      capacity_ = vec.size();
      data_ = new TYPE [capacity_];
    }

    size_ = vec.size();
    
    if(data_)
      std::copy(vec.begin(), vec.end(), data_);

    return *this;
  }


/*
  template<class TYPE>
  void StdVector<TYPE>::Assign(TYPE* source, size_type size, bool delete_previous_data)
  {
    assert(!(source == NULL && size != 0));
    assert(!(source != NULL && size == 0));

    // we may overwrite our possibly not owned data source.
    if(delete_previous_data && data_ != NULL) {
      assert(memBelongsToMe_);
      delete[] data_;
    }
    memBelongsToMe_ = true;
    data_ = source;
    size_ = size;
    capacity_ = size;
  }
*/
  // *************
  //   Push_back
  // *************
  template<class TYPE>
  void StdVector<TYPE>::Push_back_expand( const TYPE &y )
  {
    assert(!wrapped_);
    // perform memory allocation
    capacity_ = (size_ == 0) ? 1 : 2 * size_;
    TYPE* help = new TYPE[ capacity_ ];

    // copy old entries into new buffer if we have data
    std::copy(data_, data_+size_, help);

    // Perform push-back and increase size
    // note, that y might point to data, so copy before delete
    help[size_] = y;
    size_++;

    // delete old buffer and re-set pointer
    delete[] data_;
    data_ = help;
  }


  // *********
  //   Erase
  // *********
  template<class TYPE> void StdVector<TYPE>::Erase( size_type pos ) {
#ifdef CHECK_INITIALIZED
    if (size_ == 0) EXCEPTION( "Vector: Undefined Vector in function Erase" );
#endif

#ifdef CHECK_INDEX
    if (pos >=size_) EXCEPTION( "Invalid index for Erase" );
#endif

    if(pos == size_ - 1) // erasing last position
    {
      --size_;
      return;
    }

    assert(!wrapped_);
    TYPE* help=new TYPE[size_-1];
    
    for(unsigned int i=0; i < pos; i++)
      help[i] = data_[i];
    for(unsigned int i=pos+1; i < size_; i++)
      help[i-1] = data_[i];
 
    delete [] data_; 
    data_ = help;
    size_--; 
    capacity_ = size_;
  }

  template<class TYPE>
  void StdVector<TYPE>::Insert(size_type pos, const TYPE& dat)
  {
#ifdef CHECK_INDEX
    if (pos > size_) EXCEPTION( "Invalid index for Insert" );
#endif

    Resize(size_ + 1);

    for(unsigned int i = size_-1; i > pos; i--)
      data_[i] = data_[i-1];

    data_[pos] = dat;
  }

  template<class TYPE>
  void StdVector<TYPE>::Append(const StdVector<TYPE>& vec ) {

    unsigned int vecSize = vec.GetSize();
    unsigned int oldSize = size_;
    Resize(size_ + vecSize);
    std::copy(vec.data_, vec.data_+vecSize, data_+oldSize);
  }

  template<class TYPE>
  void StdVector<TYPE>::Erase(const unsigned int pos1, const unsigned int pos2)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0) {
      EXCEPTION( "Vector: Undefined Vector in function Erase" );
    }
#endif

#ifdef CHECK_INDEX
    if (pos1 < 0 || pos1 >= size_ || pos2 < 0 || pos2 >= size_) 
      EXCEPTION( "Invalid index for cut" );
    if (pos1 > pos2)
      EXCEPTION( "First index is bigger than second one in function Erase()" );
#endif
    assert(!wrapped_); 
    unsigned int l=pos2-pos1+1;
    TYPE* help=new TYPE[size_-l];
    for (unsigned int i=0; i < pos1; i++) 
      help[i]=data_[i];
    for (unsigned int i=size_-1; i > pos2 ; i--) 
      help[i-pos2+pos1-1]=data_[i];

    delete [] data_; 
    data_=help;
       
    size_ = size_ - l; 
    capacity_ = size_;
  }


  template<class TYPE>
  int StdVector<TYPE>::Find(const TYPE &x, bool quiet) const
  {
    for(unsigned int i = 0; i < size_; ++i)
      if(data_[i] == x)
        return i;

    // not found
    if(!quiet)
      EXCEPTION("cannot find " << x << " in vector of size " << size_);
    return -1;
  }

  template<class TYPE>
  int StdVector<TYPE>::Find(const TYPE &x, unsigned int& guess, bool quiet) const
  {
    int idx = -1;
    // optimistically the guess is right, search upwards
    for(unsigned int i = guess; idx < 0 && i < size_; i++)
      if(data_[i] == x)
        idx = i;

    // in case we did not find, search from the beginning
    for(unsigned int i = 0; idx < 0 && i < guess; i++)
      if(data_[i] == x)
        idx = i;

    if(idx > -1) {
      guess = idx+1 == (int) size_ ? 0 : idx+1;
      return idx;
    }

    if(!quiet)
      EXCEPTION("cannot find " << x << " in vector of size " << size_);
    return -1;

  }


  template<class TYPE>
  StdVector<unsigned int> StdVector<TYPE>::FindAll(const TYPE &x) const
  {
    StdVector<unsigned int> t (0);
    for(unsigned int i = 0; i < size_; ++i){
      if(data_[i] == x) t.Push_back(i);
    }
    return t;
  }


  template<class TYPE>
  bool StdVector<TYPE>::IsUnique() const
  {
    // possibly not the fastest algorithm as we check any pair
    for(unsigned int s = 0; s < size_; s++) // slow variable
      for(unsigned int f = s+1; f < size_; f++) // fast variable
        if(data_[s] == data_[f])
          return false;

    return true;
  }


  template<class TYPE>
  bool StdVector<TYPE>::operator== (const StdVector<TYPE> & vec) const
  {
    if (size_ == 0 && vec.size_ == 0)
      return true;
  
    for (unsigned int i = 0; i < size_; i++)
      if (data_[i] != vec.data_ [i])
        return false;
  
    return true;
  }
  
  template<class TYPE>
  bool StdVector<TYPE>::operator== (const TYPE* other) const
  {
    if (size_ == 0 && other == NULL)
      return true;
  
    if(other == NULL)
      EXCEPTION("Cannot compare non-empty StdVector with NULL");
    
    if(size_ == 0)
      return false;
    
    for (unsigned int i = 0; i < size_; i++)
      if (data_[i] != other[i])
        return false;
  
    return true;
  }

  
  
  template<class TYPE>
  bool StdVector<TYPE>::operator!= (const StdVector<TYPE> & vec) const
  {
    return ! ( *this==vec );
  }

  template<class TYPE>
  std::string StdVector<TYPE>::Serialize( char separator ) const {
    std::stringstream out;

    if( size_ > 0 ) {
      for( unsigned int i = 0; i < size_-1; i++ ) {
        out << data_[i] << separator << " ";
      }
      out << data_[size_-1];
    }
    return out.str();
  }


  
  template<class TYPE>
  std::string StdVector<TYPE>::ToString(int size, const TYPE* data, ToStringFormat format, const std::string& sep_in)
  {
    std::ostringstream os;

    std::string sep = sep_in != "" ? sep_in : (format == TS_PYTHON ? ", " : " ");

    switch(format)
    {
    case TS_INFO:
      os << " size=" << size;
      break;

    case TS_MATLAB:
    case TS_PYTHON:
      os << "[";
      // intentionally no break;
    case TS_PLAIN:
      for(int i = 0; i < size; i++)
      {
        os << data[i];

        if(i < size-1)
          os << sep;
      }
      if(format != TS_PLAIN)
        os << "]";
      break;
      
    case TS_NONZEROS:
      throw Exception("TS_NONZEROS is no valid format for StdVector::ToString()");
      break;
    }
    return os.str();
  }


  template<class TYPE>
  bool StdVector<TYPE>::InWindow(unsigned int index)
  {
    return index >= window.GetStart() && index < window.GetStart() + window.GetSize();
  }

  template<class TYPE>
  StdVector<TYPE>::Window::Window()
  {
    this->active_ = false;
#ifdef _WIN32
    this->start_  = UINT_MAX;
    this->size_   = UINT_MAX;
#else
    this->start_  = std::numeric_limits<unsigned int>::max();
    this->size_   = std::numeric_limits<unsigned int>::max();
#endif
  }

  template<class TYPE>
  void StdVector<TYPE>::Window::Set(unsigned int start, unsigned int size)
  {
    this->active_ = true;
    this->start_  = start;
    this->size_   = size;
  }

  template<class TYPE>
  void StdVector<TYPE>::Window::Set(StdVector<TYPE>& vec)
  {
    this->active_ = true;
    this->start_  = 0;
    this->size_   = vec.GetSize();
  }


  template<class TYPE>
  unsigned int StdVector<TYPE>::Window::GetStart() const
  {
    #ifdef CHECK_INITIALIZED
      if(!active_)
        EXCEPTION("Vector: window not initialized." );
    #endif
    return start_;
  }

  template<class TYPE>
  unsigned int StdVector<TYPE>::Window::GetSize() const
  {
    #ifdef CHECK_INITIALIZED
      if(!active_)
        EXCEPTION("Vector: window not initialized." );
    #endif
    return size_;
  }

  template<class TYPE>
  std::string StdVector<TYPE>::ToString(ToStringFormat format, const std::string& sep, bool in_window) const
  {
    if(!in_window)
      return StdVector<TYPE>::ToString(size_, data_, format, sep);
    else
    {
      #ifdef CHECK_INDEX
        if(window.GetStart() + window.GetSize() > size_)
          EXCEPTION("Vector: window bounds violated." );
      #endif
      return StdVector<TYPE>::ToString(window.GetSize(), data_ + window.GetStart(), format, sep);
    }
  }
  
  template<class TYPE>  
  void StdVector<TYPE>::ToString(StdVector<std::string>& out) const
  {
    out.ResizeNoCopy(size_);
    for(unsigned int i = 0; i < size_; i++)
      out[i] = boost::lexical_cast<std::string>(data_[i]);
  }

  template<class TYPE>
  void StdVector<TYPE>::Parse(const StdVector<std::string>& in)
  {
    ResizeNoCopy(in.GetSize());

    for(unsigned int i = 0; i < in.GetSize(); i++)
      data_[i] = boost::lexical_cast<TYPE>(in[i]);
  }
  
  template<class TYPE>
  void StdVector<TYPE>::Swap(size_type idx1, size_type idx2)
  {
    TYPE tmp = data_[idx1];
    data_[idx1] = data_[idx2];
    data_[idx2] = tmp;
  }


  template<class S>
  std::ostream & operator << ( std::ostream & out, const StdVector<S> & vc)
  {
    for (unsigned int i=0; i < vc.GetSize(); i++)
      out << vc[i] << " " << std::endl;
    return out;
  }  
  
} // end of namespace

#endif
