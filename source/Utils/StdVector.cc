// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STDVECTOR_IMPLEMENTATION_2004
#define FILE_STDVECTOR_IMPLEMENTATION_2004

#include <iostream>
#include <boost/lexical_cast.hpp>

namespace CoupledField {

  template<class TYPE>
  StdVector<TYPE>::StdVector() :
    size_(0),
    capacity_(0),
    data_(NULL)
  { }

  template<class TYPE>
  StdVector<TYPE>::StdVector(unsigned int size) :
    size_(size),
    capacity_(size),
    data_(new TYPE [size])
  {
    Init();
  }

  template<class TYPE>
  StdVector<TYPE>::StdVector(const StdVector<TYPE> & vec) :
    size_(vec.size_),
    capacity_(vec.size_),
    data_(new TYPE [vec.size_])
  {
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] =  vec.data_[i];
  }

  template<class TYPE>
  StdVector<TYPE>::StdVector(const std::vector<TYPE> & vec) :
    size_(vec.size()),
    capacity_(vec.size()),
    data_(new TYPE [vec.size()])
  {
    //Warning(" This function should not be used anymore!");
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] =  vec[i];
  }

  template<class TYPE>
  StdVector<TYPE>::~StdVector()
  {
    delete[] data_;
  }

  template<class TYPE>
  void StdVector<TYPE>::Clear()
  {
    size_ = 0;
    capacity_ = 0;
    delete[] data_;
    data_ = NULL;
  }

  template<class TYPE>
  void StdVector<TYPE>::Init(const TYPE entry)
  {
    for(unsigned int i = 0; i < size_; ++i) 
      data_[i] = entry;
  }

  template<class TYPE>
  void StdVector<TYPE>::Reserve(unsigned int capacity)
  {
    unsigned int old_size = size_;

    Resize(capacity);

    size_ = old_size;
  }

  template<class TYPE>
  void StdVector<TYPE>::Resize(const unsigned int size)
  {
#ifdef CHECK_INITIALIZED
    if(capacity_ < size_)
      EXCEPTION("capacity " << capacity_ << " smaller size " << size_);
#endif

    // the cheap case, e.g. Resize(0)
    if(size <= capacity_)
    {
      size_ = size;
    }
    else
    {
      TYPE* help = new TYPE[size];

      for (unsigned int i=0; i<size_; i++)
        help[i] = data_[i];
      
      // interestingly std::copy(data_, data_ +  sizeof(TYPE) * size_, help)
      // does not work. Due to copy constructor stuff for complex types??

      delete[] data_;
      data_ = help;
      size_ = size;
      capacity_ = size;
    }
  }

  template<class TYPE>
  void StdVector<TYPE>::Resize(const unsigned int size, TYPE entry)
  {
    Resize(size);
    Init(entry);
  }

  
  // *************
  //   operator=
  // *************
  template<class TYPE>
  StdVector<TYPE>& StdVector<TYPE>::operator= ( const StdVector &vec ) {

    // Avoid self-assignements
    if ( this == &vec ) {
      return *this;
    }

    // Vectors are different

    // If there is not enough space to copy the entries
    // perform a re-allocation
    if(capacity_ < vec.size_)
    {
      // Delete old buffer
      delete [] data_;

      // Allocate new buffer
      capacity_ = vec.size_;
      data_     = new TYPE[capacity_];
    }

    // Copy entries
    size_ = vec.size_;
    for(unsigned int i = 0; i < size_; ++i)
      data_[i] = vec.data_[i];

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
      delete [] data_;
      capacity_ = vec.size();
      data_ = new TYPE [capacity_];
    }

    size_ = vec.size();
    
    for(unsigned int i = 0; i < size_; ++i)
      data_ [i] = vec[i];

    return *this;
  }

  template<class TYPE>
  void StdVector<TYPE>::Import(const TYPE* source, unsigned int size)
  {
    if(source == NULL) EXCEPTION("cannot import NULL");

    // Reserve and Resize do stuff we don't need
    if(data_ != NULL) delete[] data_;
    data_ = new TYPE[size];
    size_ = size;
    capacity_ = size;
    
    for(unsigned int i = 0; i < size; i++)
      data_[i] = source[i];
  }
  

  // *************
  //   Push_back
  // *************
  template<class TYPE>
  void StdVector<TYPE>::Push_back( const TYPE &y )
  {
    // Check whether capacity is sufficiently large to perform
    // a push-back operation. If not allocate memory according
    // to the following simply scheme: Each time capacity is
    // exceeded allocate twice as much memory as there was before.
    if ( size_ >= capacity_ ) {

      TYPE *help;

      // perform memory allocation
      capacity_ = size_ == 0 ? 1 : 2 * size_;
      help = new TYPE[ capacity_ ];

      // copy old entries into new buffer
      for ( unsigned int i = 0; i < size_; i++ ) {
        help[i] = data_[i];
      }

      // Perform push-back and increase size
      // note, that y might point to data, so copy before delete
      help[size_] = y;
      size_++;

      // delete old buffer and re-set pointer
      delete[] data_;
      data_ = help;

    } else {
      // Perform push-back and increase size
      data_[size_] = y;
      size_++;
    }
  }


  // *********
  //   Erase
  // *********
  template<class TYPE> void StdVector<TYPE>::Erase( const unsigned int pos ) {
#ifdef CHECK_INITIALIZED
    if (size_ == 0) {
      EXCEPTION( "Vector: Undefined Vector in function Erase" );
    }
#endif

#ifdef CHECK_INDEX
    if (pos >=size_) 
      EXCEPTION( "Invalid index for cut" );
#endif

    if(pos == size_ - 1) // erasing last position
    {
      --size_;
      return;
    }
    
    unsigned int i;
 
    TYPE * help=new TYPE[size_-1];
    
    for (i=0; i < pos; i++) 
      help[i] = data_[i];
    for (i=pos+1; i < size_; i++) 
      help[i-1] = data_[i];
 
    delete [] data_; 
    data_ = help;
    size_--; 
    capacity_ = size_;
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
    unsigned int i;
 
    unsigned int l=pos2-pos1+1;
    TYPE * help=new TYPE[size_-l];
    for (i=0; i < pos1; i++) 
      help[i]=data_[i];
    for (i=size_-1; i > pos2 ; i--) 
      help[i-pos2+pos1-1]=data_[i];

    delete [] data_; 
    data_=help;
   
    size_ = size_ - l; 
    capacity_ = size_;
  }

  template<class TYPE>
  int StdVector<TYPE>::Find(const TYPE &x) const
  {
    for(unsigned int i = 0; i < size_; ++i)
      if(data_[i] == x) return i;
    
    // not found
    return -1;
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
  std::string StdVector<TYPE>::ToString(int size, const TYPE* data, int level, int stride)
  {
    std::ostringstream os;

    switch(level)
    {
    case 0:

      for(int i = 0; i < size; i += stride)
      {
        os << data[i];
        if(i < size-1) os << ", ";
      }
      break;

    default:
      
      os << "size=" << size;
      if(size > 0)
      {
        // todo: We have no min/max in the complex case!

        /*
        TYPE min = data[0];
        TYPE max = data[0];

        for(int i = 0; i < size; i++)
        {
          min = std::min(min, data[i]);
          max = std::max(max, data[i]);
        }
        os << " min=" << min << " max=" << max;
        */        
      }
      break;
    }

    return os.str();
  }

  template<class TYPE>
  std::string StdVector<TYPE>::ToString(int level, int stride) const
  {
    return StdVector<TYPE>::ToString(size_, data_, level, stride);
  }
  
  template<class TYPE>  
  void StdVector<TYPE>::ToString(StdVector<std::string>& out) const
  {
    out.Resize(size_);
    for(unsigned int i = 0; i < size_; i++)
      out[i] = boost::lexical_cast<std::string>(data_[i]);
  }

  template<class TYPE>
  void StdVector<TYPE>::Parse(const StdVector<std::string>& in)
  {
    Resize(in.GetSize());

    for(unsigned int i = 0; i < in.GetSize(); i++)
      data_[i] = boost::lexical_cast<TYPE>(in[i]);
  }
  
  template<class TYPE>
  void StdVector<TYPE>::Swap(unsigned int idx1, unsigned int idx2)
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
