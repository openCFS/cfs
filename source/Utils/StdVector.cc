#ifndef FILE_STDVECTOR_IMPLEMENTATION_2004
#define FILE_STDVECTOR_IMPLEMENTATION_2004

#include "DataInOut/WriteInfo.hh"

namespace CoupledField {

  template<class TYPE>
  StdVector<TYPE>::StdVector()
  {
    ENTER_IFCN( "StdVector::StdVector()" );

    size_ = 0;
    capacity_ = 0;
    data_ = NULL;
  }

  template<class TYPE>
  StdVector<TYPE>::StdVector(UInt size)
  {
    ENTER_IFCN( "StdVector::StdVector(UInt)" );

    size_ = size;
    capacity_ = size;
    data_ = new TYPE [size];
  
    for (UInt i = 0; i < size; i++)
      data_ [i] = TYPE();
  }

  template<class TYPE>
  StdVector<TYPE>::StdVector(const StdVector<TYPE> & vec)
  {
    ENTER_IFCN( "StdVector::StdVector(StdVector)" );
  
    size_ = vec.size_;
    capacity_ = vec.size_;
    data_ = new TYPE [size_];

    for (UInt i = 0; i < size_; i++)
      data_[i] =  vec.data_[i];
  }

  template<class TYPE>
  StdVector<TYPE>::StdVector(const std::vector<TYPE> & vec)
  {
    ENTER_IFCN( "StdVector::StdVector(std::vector)" );

    //Warning(" This function should not be used anymore!");
  
    size_ = vec.size();
    capacity_ = vec.size();
    data_ = new TYPE[size_];

    for (UInt i = 0; i < size_; i++)
      data_[i] =  vec[i];
  }

  template<class TYPE>
  StdVector<TYPE>::~StdVector()
  {
    ENTER_IFCN( "StdVector::~StdVector" );

    if (data_)
      delete[] data_;
  }

  template<class TYPE>
  void StdVector<TYPE>::Clear()
  {
    ENTER_IFCN( "StdVector::Clear" );

    size_ = 0;
    capacity_ = 0;
    if (data_)
      delete[] data_;
    data_ = NULL;
  }

  template<class TYPE>
  void StdVector<TYPE>::Init(const TYPE entry)
  {
    ENTER_IFCN( "StdVector::Init" );

    // #ifdef CHECK_INITIALIZED
    //   if (size_ == 0) 
    //     Warning("Don't use Init() to undefined vector", 
    //      __FILE__, __LINE__);
    // #endif
  
    for (UInt i=0; i<size_; i++) 
      data_[i]=entry;
  }

  template<class TYPE>
  void StdVector<TYPE>::Reserve(UInt capacity)
  {
    ENTER_IFCN( "StdVector::Reserve" );

#ifdef CHECK_INDEX
    if (capacity < 0)
      EXCEPTION( "Invalid value for StdVector::Reserve" );
#endif
  
    if (capacity > capacity_)
      {
        capacity_ = capacity;
      
        TYPE * help = new TYPE[capacity];
      
        for (UInt i=0; i<size_; i++)
          help[i] = data_[i];
      
        delete[] data_;
        data_ = help;
      }
  }

  template<class TYPE>
  void StdVector<TYPE>::Resize(const UInt size)
  {
    ENTER_IFCN( "StdVector::Resize" );

#ifdef CHECK_INDEX
    if (size < 0) 
      EXCEPTION( "Invalid dimension for Resize" );
#endif  
  
    if (size != size_ ||  capacity_ < size)
      {

        TYPE* help = new TYPE[size];
        
        if(size > size_) {
          for( UInt i = 0; i < size_; i++ ) 
            help[i] = data_[i];
        } else {
          for( UInt i = 0; i < size; i++ ) 
            help[i] = data_[i];
        }
        
        delete[] data_;
        data_ = help;
        size_ = size;
        capacity_ = size;
      }
  

  }


  // *************
  //   operator=
  // *************
  template<class TYPE>
  StdVector<TYPE>& StdVector<TYPE>::operator= ( const StdVector &vec ) {

    ENTER_IFCN( "StdVector::operator=(StdVector)" );

    // Avoid self-assignements
    if ( this == &vec ) {
      return *this;
    }

    // Vectors are different
    else {

      // If there is not enough space to copy the entries
      // perform a re-allocation
      if ( capacity_ < vec.size_ ) {

        // Delete old buffer
        if ( data_ != NULL ) {
          delete [] data_;
        }

        // Allocate new buffer
        capacity_ = vec.size_;
        data_     = new TYPE[capacity_];
      }

      // Copy entries
      size_ = vec.size_;
      for ( UInt i = 0; i < size_; i++ ) {
        data_[i] = vec.data_[i];
      }
    }

    return *this;

  }


  // *************
  //   operator=
  // *************
  template<class TYPE>
  StdVector<TYPE>& StdVector<TYPE>::operator= ( const std::vector<TYPE> &vec){
    ENTER_IFCN( "StdVector::operator=(const std::vector)" );

    if (size_ != vec.size())
      { 
        if (data_)
          delete [] data_;
      
        size_ = vec.size();
        capacity_ = vec.size();
        data_ = new TYPE [size_];
      }
  
    for (UInt i = 0; i < size_; i++)
      data_ [i] = vec[i];
  
    return *this;
  
  }
  
  template<class TYPE>
  void StdVector<TYPE>::Insert(const TYPE & y, UInt pos)
  {
    ENTER_IFCN( "StdVector::Insert" );
  
#ifdef CHECK_INDEX
    if (pos >= size_)
      EXCEPTION( "StdVector::Insert: Index out of bounds" );
#endif

    TYPE * help = new TYPE[size_+1];
    size_ = size_+1;
    capacity_ = size_;

    for (UInt i=0; i<pos; i++)
      help[i] = data_[i];

    help[pos] = y;

    for (UInt i=pos+1; i<size_; i++)
      help[i] = data_[i-1];

    delete[] data_;
    data_ = help;
  }

  template<class TYPE>
  void StdVector<TYPE>::Insert(const TYPE & y, UInt pos, UInt numCopies)
  {
    ENTER_IFCN( "StdVector::Insert" );
  
#ifdef CHECK_INDEX
    if (pos >= size_)
      EXCEPTION( "StdVector::Insert: Index out of bounds" );

    if (numCopies < 1)
      EXCEPTION( "StdVector::Insert: NumCopies must be greater 1" );
#endif

    TYPE * help = new TYPE[size_+numCopies];
    size_ = size_+numCopies;
    capacity_ = size_;

    for (UInt i=0; i<pos; i++)
      help[i] = data_[i];

    for (UInt i=pos; i<pos+numCopies; i++)   
      help[i] = y;

    for (UInt i=pos+numCopies; i<size_; i++)
      help[i] = data_[i-1];

    delete[] data_;
    data_ = help;
  }


  // *************
  //   Push_back
  // *************
  template<class TYPE>
  void StdVector<TYPE>::Push_back( const TYPE &y ) {

    ENTER_IFCN( "StdVector::Push_back" );
  
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
      for ( UInt i = 0; i < size_; i++ ) {
        help[i] = data_[i];
      }

      // Perform push-back and increase size
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
  template<class TYPE> void StdVector<TYPE>::Erase( const UInt pos ) {
    ENTER_IFCN( "StdVector::Erase" );

#ifdef CHECK_INITIALIZED
    if (size_ == 0) {
      EXCEPTION( "Vector: Undefined Vector in function Erase" );
    }
#endif

#ifdef CHECK_INDEX
    if (pos<0 || pos >=size_) 
      EXCEPTION( "Invalid index for cut" );
#endif

    UInt i;
 
    TYPE * help=new TYPE[size_-1];
    for (i=0; i < pos; i++) 
      help[i]=data_[i];
    for (i=pos+1; i < size_; i++) 
      help[i-1]=data_[i];
 
    delete [] data_; 
    data_=help;
    size_--; 
    capacity_--;
  }


  template<class TYPE>
  void StdVector<TYPE>::Erase(const UInt pos1, const UInt pos2)
  {
    ENTER_IFCN( "StdVector::Erase" );
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
    UInt i;
 
    UInt l=pos2-pos1+1;
    TYPE * help=new TYPE[size_-l];
    for (i=0; i < pos1; i++) 
      help[i]=data_[i];
    for (i=size_-1; i > pos2 ; i--) 
      help[i-pos2+pos1-1]=data_[i];

    delete [] data_; 
    data_=help;
   
    size_ = size_ - l; 
  }

  template<class TYPE>
  Integer StdVector<TYPE>::Find(const TYPE &x) const
  {
    ENTER_IFCN( "StdVector::Find" );

    Integer pos = -1;

    for (UInt i=0; i<size_; i++)
      if (data_[i] == x)
      {
        pos = i;
        break;
      }
    
    return pos;
  }

  template<class TYPE>
  bool StdVector<TYPE>::operator== (const StdVector<TYPE> & vec) const
  {
    ENTER_IFCN( "StdVector::operator==" );
    // #ifdef CHECK_INITIALIZED 
    //   if ((size_ == 0) || (vec.size_ == 0))
    //     Warning("Vector: undefined Vector in operator==",__FILE__,__LINE__);
    // #endif

    if (size_ == 0 && vec.size_ == 0)
      return true;
  
    for (UInt i = 0; i < size_; i++)
      if (data_[i] != vec.data_ [i])
        return false;
  
    return true;
  }
  

  template<class TYPE>
  bool StdVector<TYPE>::operator!= (const StdVector<TYPE> & vec) const
  {
    ENTER_IFCN( "StdVector::operator!=" );
    return ! ( *this==vec );
  }

  template<class TYPE>
  std::string StdVector<TYPE>::Serialize( Char separator ) const {
    std::stringstream out;

    if( size_ > 0 ) {
      for( UInt i = 0; i < size_-1; i++ ) {
        out << data_[i] << separator << " ";
      }
      out << data_[size_-1];
    }
    return out.str();
  }
        
  
  template <class S> 
  void Sort(S* v, UInt n)
  {
    ENTER_IFCN( "StdVector::Sort" );
  }

  template<class S>
  std::ostream & operator << ( std::ostream & out, const StdVector<S> & vc)
  {
    ENTER_IFCN( "operator <<(StdVector)" );

    for (UInt i=0; i < vc.GetSize(); i++)
      out << vc[i] << " " << std::endl;
    return out;
  }  
  
  template<class TYPE> template< class Archive>
  void StdVector<TYPE>::save(Archive & ar, const unsigned int version) const {
    ar & size_;
    ar & capacity_;
    for( UInt i = 0; i < size_; i++ ) 
      ar & data_[i];
  }
  
  template<class TYPE> template <class Archive>
  void StdVector<TYPE>::load(Archive & ar, const unsigned int version) {
    if( data_ != NULL ) {
      delete[] data_;
    }
    
    ar & size_;
    ar & capacity_;
    data_ = new TYPE[capacity_];
    for( UInt i = 0; i < size_; i++ ) {
      ar & data_[i];
    }
  }
  
} // end of namespace

#endif
