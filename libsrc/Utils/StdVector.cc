#ifndef FILE_STDVECTOR_IMPLEMENTATION_2004
#define FILE_STDVECTOR_IMPLEMENTATION_2004

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
StdVector<TYPE>::StdVector(Integer size)
{
  ENTER_IFCN( "StdVector::StdVector(Integer)" );

  size_ = size;
  capacity_ = size;
  data_ = new TYPE [size];
  
  for (Integer i = 0; i < size; i++)
    data_ [i] = TYPE();
}

template<class TYPE>
StdVector<TYPE>::StdVector(const StdVector<TYPE> & vec)
{
  ENTER_IFCN( "StdVector::StdVector(StdVector)" );
  
  size_ = vec.size_;
  capacity_ = vec.size_;
  data_ = new TYPE [size_];

  for (Integer i = 0; i < size_; i++)
    data_[i] =  vec.data_[i];
}

template<class TYPE>
StdVector<TYPE>::StdVector(const std::vector<TYPE> & vec)
{
  ENTER_IFCN( "StdVector::StdVector(std::vector)" );

  Warning(" This function should not be used anymore!");
  
  size_ = vec.size();
  capacity_ = vec.size();
  data = new TYPE[size_];

  for (Integer i = 0; i < size_; i++)
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
// 	    __FILE__, __LINE__);
// #endif
  
  for (Integer i=0; i<size_; i++) 
    data_[i]=entry;
}

template<class TYPE>
void StdVector<TYPE>::Reserve(Integer capacity)
{
  ENTER_IFCN( "StdVector::Reserve" );

#ifdef CHECK_INDEX
  if (capacity < 0)
    Info->Error("Invalid value for StdVector::Reserve",
		__FILE__, __LINE__);
#endif
  
  if (capacity > capacity_)
    {
      capacity_ = capacity;
      
      TYPE * help = new TYPE[capacity];
      
      for (Integer i=0; i<size_; i++)
	help[i] = data_[i];
      
      delete[] data_;
      data_ = help;
    }
}

template<class TYPE>
void StdVector<TYPE>::Resize(const Integer size)
{
  ENTER_IFCN( "StdVector::Resize" );

#ifdef CHECK_INDEX
  if (size < 0) 
    Error("Invalid dimension for Resize", __FILE__, __LINE__);
#endif  
  
  if (size != size_ ||  capacity_ < size)
    {
      if (data_) delete[] data_;
      
      size_ = size;
      capacity_ = size;
      data_ = new TYPE[size_];
    }
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] = TYPE();
}
 
template<class TYPE>
StdVector<TYPE>& StdVector<TYPE>::operator= (const StdVector & vec)
{
  ENTER_IFCN(" StdVector::operator=(StdVector)" );

  if (this == &vec)
    return *this;


  if (size_ != vec.size_)
    {	
      if (data_)
	delete [] data_;
      
      size_ = vec.size_;
      capacity_ = vec.size_;
      data_ = new TYPE [size_];
    }
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] = vec.data_[i];
  
  return *this;
}

template<class TYPE>
StdVector<TYPE>& StdVector<TYPE>::operator= (const std::vector<TYPE> & vec)
{
  ENTER_IFCN( "StdVector::operator=(const std::vector)" );

  if (size_ != vec.size())
    {	
      if (data_)
	delete [] data_;
      
      size_ = vec.size();
      capacity_ = vec.size();
      data_ = new TYPE [size_];
    }
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] = vec[i];
  
  return *this;
  
}
  
template<class TYPE>
void StdVector<TYPE>::Insert(const TYPE & y, Integer pos)
{
  ENTER_IFCN( "StdVector::Insert" );
  
#ifdef CHECK_INDEX
  if (pos >= size_)
    Error( "StdVector::Insert: Index out of bounds", __FILE__, __LINE__);
#endif

  TYPE * help = new TYPE[size_+1];
  size_ = size_+1;
  capacity_ = size_;

  for (Integer i=0; i<pos; i++)
    help[i] = data_[i];

  help[pos] = y;

  for (Integer i=pos+1; i<size_; i++)
    help[i] = data_[i-1];

  delete[] data_;
  data_ = help;
}

template<class TYPE>
void StdVector<TYPE>::Insert(const TYPE & y, Integer pos, Integer numCopies)
{
  ENTER_IFCN( "StdVector::Insert" );
  
#ifdef CHECK_INDEX
  if (pos >= size_)
    Error( "StdVector::Insert: Index out of bounds", __FILE__, __LINE__);

  if (numCopies < 1)
    Error( "StdVector::Insert: NumCopies must be greater 1",
	   __FILE__, __LINE__);
#endif

  TYPE * help = new TYPE[size_+numCopies];
  size_ = size_+numCopies;
  capacity_ = size_;

  for (Integer i=0; i<pos; i++)
    help[i] = data_[i];

  for (Integer i=pos; i<pos+numCopies; i++)   
    help[i] = y;

  for (Integer i=pos+numCopies; i<size_; i++)
    help[i] = data_[i-1];

  delete[] data_;
  data_ = help;
}

template<class TYPE>
void StdVector<TYPE>::Push_back(const TYPE & y)
{
  ENTER_IFCN( "StdVector::Push_back" );
  
  if (size_ >= capacity_)
    {
      TYPE * help;
      if (size_ == 0)
	{
	  help = new TYPE[1];
	  capacity_ = 1;
	}
      else 
	{
	  help = new TYPE[size_*2];
	  capacity_ = size_*2;
	}
      
      for (Integer i=0; i<size_; i++)
	help[i] = data_[i];
      
      delete[] data_;
      data_ = help;
    } 
  
  data_[size_] = y;
  size_++;
  
}

template<class TYPE>
void StdVector<TYPE>::Erase(const Integer pos)
{
  ENTER_IFCN( "StdVector::Erase" );
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Warning("Vector: undefined Vector in function Erase", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
   if (pos<0 || pos >=size_) 
     Error("Invalid index for cut");
#endif

  Integer i;
 
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
void StdVector<TYPE>::Erase(const Integer pos1, const Integer pos2)
{
  ENTER_IFCN( "StdVector::Erase" );
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Warning("Vector: undefined Vector in function Erase()", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
   if (pos1 < 0 || pos1 >= size_ || pos2 < 0 || pos2 >= size_) 
     Error("Invalid index for cut");
   if (pos1 > pos2)
     Error("First index is bigger than second one in function Erase()",
	   __FILE__, __LINE__);
#endif
   Integer i;
 
   Integer l=pos2-pos1+1;
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

  for (Integer i=0; i<size_; i++)
    if (data_[i] == x)
      pos = i;
  return pos;
}

template<class TYPE>
Boolean StdVector<TYPE>::operator== (const StdVector<TYPE> & vec) const
{
  ENTER_IFCN( "StdVector::operator==" );
// #ifdef CHECK_INITIALIZED 
//   if ((size_ == 0) || (vec.size_ == 0))
//     Error("Vector: undefined Vector in operator==",__FILE__,__LINE__);
// #endif

  if (size_ == 0 && vec.size_ == 0)
    return TRUE;
  
  for (Integer i = 0; i < size_; i++)
    if (data_[i] != vec.data_ [i])
      return FALSE;
  
  return TRUE;
}
  

template<class TYPE>
Boolean StdVector<TYPE>::operator!= (const StdVector<TYPE> & vec) const
{
  ENTER_IFCN( "StdVector::operator!=" );
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator !=", __FILE__, __LINE__);
#endif
  
  for (Integer i = 0; i < size_; i++)
    if (data_[i] != vec.data_[i])
      return FALSE;
  
  return TRUE;
}
  
template <class S> 
void Sort(S* v, Integer n)
{
  ENTER_IFCN( "StdVector::Sort" );
}

template<class S>
std::ostream & operator << ( std::ostream & out, const StdVector<S> & vc)
{
  ENTER_IFCN( "operator <<(StdVector)" );

for (Integer i=0; i < vc.GetSize(); i++)
  out << vc[i] << " " << std::endl;
 return out;
}  
  
} // end of namespace

#endif
