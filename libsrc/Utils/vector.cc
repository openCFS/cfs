#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

#include "vector.hh"
#include <Matrix/matrix.hh>


namespace CoupledField
{ 

template<class TYPE>
Vector<TYPE>::Vector()
  : CFSVector()
{
  ENTER_IFCN("Vector::Vector");
  data_ = NULL;
  size_ = 0;
}

template<class TYPE> 
Vector<TYPE>::Vector(Integer size)
  : CFSVector(size)
{
  ENTER_IFCN("Vector::Vector");
  size_ = size;
  data_ = new TYPE [size];

#ifdef CHECK_MEMORY
  if (!data_)
    Error("Vector: out of memory",__FILE__,__LINE__);
#endif

  for (Integer i = 0; i < size; i++)
    data_ [i] = 0;
}

template<class TYPE> 
Vector<TYPE>::Vector(const Vector<TYPE> &vec)
{
  ENTER_IFCN("Vector::Vector");    
  size_ = vec.size_;
  data_ = new TYPE [vec.size_];

#ifdef CHECK_MEMORY
  if (!data_)
    Error("Vector: out of memory",__FILE__,__LINE__);
#endif
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] =  vec.data_[i];
}
template<class TYPE> 
Vector<TYPE>::Vector(const std::vector<TYPE>  &vec)
{
  ENTER_IFCN("Vector::Vector");    

  // ** TO IMPLEMENT **
}

template<class TYPE>
Vector<TYPE>::~Vector()
{
  ENTER_IFCN("Vector::~Vector");
  if (data_)
    delete[] data_;
}

// Two different methods for getting a double pointer to the data
// in the vector:
// 1.) Double-Vectors: Here simply the pointer data_ is passed
// 2.) Complex-Vectors: A new array is created, which holds the 
//                      complex values as a serial sequence of
//                      Double values.
Double* Vector<Double>::GetDoublePointer()
{
  ENTER_IFCN("Vector::GetDoublePointer");

  return data_;
}

Double* Vector<Integer>::GetDoublePointer()
{
  ENTER_IFCN("Vector::GetDoublePointer");

  Error("Vector<Integer>::GetDoublePointer: Not implemented!",__FILE__,__LINE__);
  return NULL;
}

Double* Vector<Complex>::GetDoublePointer()
{
  ENTER_IFCN("Vector::GetDoublePointer");

  Double * help = new Double[size_ * 2];

#ifdef CHECK_MEMORY  
  if (help == NULL) Error("out of memory",__FILE__,__LINE__); 
#endif
  
  for (Integer i=0; i<size_; i++)
    {
     help[2*i]   = data_[i].real();
     help[2*i+1] = data_[i].imag(); 
    }

  return help;
}

template<class TYPE>
void Vector<TYPE>::Init(const TYPE entry)
{
  ENTER_IFCN("Vector::Init");
#ifdef CHECK_INITIALIZED
  if (size_ == 0) Error("Don't use Init() to undefined vector", __FILE__, __LINE__);
#endif

 for (Integer i=0; i<size_; i++) 
   data_[i]=entry;
}


template<class TYPE> 
void Vector<TYPE>::Resize(const Integer size)
{
  ENTER_IFCN("Vector::Resize");
#ifdef CHECK_INDEX
  if (size <= 0) Error("invalid dimension for Resize", __FILE__, __LINE__);
#endif  
  
  if (size != size_)
    {
      if (data_) delete[] data_;
      
      size_=size;
      data_ = new TYPE[size_];
#ifdef CHECK_MEMORY     
      if (!data_) Error("Out of memory in vector.Resize()");
#endif
    }
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] = 0;
}


template<class TYPE> 
void Vector<TYPE>::SetEntry(const Integer i, const TYPE &s)
{
  ENTER_IFCN("Vector::SetEntry");
#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector: invalid access to element " << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__,__LINE__);
#endif
    data_[i] = s;
}
  
template<class TYPE> 
void Vector<TYPE>::GetEntry(const Integer i, TYPE &ret) const
{
  ENTER_IFCN("Vector::GetEntry");  
#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector: invalid access to element " << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__,__LINE__);
#endif
    ret=  data_[i];
}

template<class TYPE> 
void Vector<TYPE>::AddEntry(const Integer i, const TYPE &s)
{
  ENTER_IFCN("Vector::AddEntry");
#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector: invalid access to element " << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__,__LINE__);
#endif
    data_[i]+=s;
}

template<class TYPE> 
void Vector<TYPE>::MultEntry(const Integer i, const TYPE &s)
{
  ENTER_IFCN("Vector::MultEntry");
#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector::MultEntry: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__,__LINE__);
#endif
    data_[i]+=s;
}

template<class TYPE> 
void Vector<TYPE>::MultAddEntry(const Integer i, const TYPE &a, const TYPE &s)
{
  ENTER_IFCN("Vector::MuldAddEntry");
#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector::MultEntry: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__,__LINE__);
#endif
    data_[i]=a*data_[i] + s;
}

template<class TYPE> 
void Vector<TYPE>::Add(const CFSVector& y)
{
  ENTER_IFCN("Vector::Add");
  TRY_CAST

  CONSTREFCAST(y,Vector<TYPE>,vec);
#ifdef CHECK_INDEX
  if (size_ != vec.size_)
    Error("Vector: incompatible dimension for operator Add(Basevector)",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<size_; i++)
    data_[i]+=vec.data_[i];

  CATCH_CAST
}

template<class TYPE> 
void Vector<TYPE>::Add(const TYPE a, const CFSVector &y)
{
  ENTER_IFCN("Vector::Add"); 

  TRY_CAST
  CONSTREFCAST(y,Vector<TYPE>,vec);
#ifdef CHECK_INDEX
  if (size_ != vec.size_)
    Error("Vector: incompatible dimension for operator Add(TYPE,Basevector)",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<size_; i++)
    data_[i]+=a*vec.data_[i];

  CATCH_CAST
}

template<class TYPE> 
void Vector<TYPE>::Add( const TYPE a, const CFSVector& y,
		     const TYPE b, const CFSVector& z )
{
  ENTER_IFCN("Vector::Add");  

  TRY_CAST
  CONSTREFCAST(y,Vector<TYPE>,vec);
#ifdef CHECK_INDEX
  if (size_ != vec.size_)
    Error("Vector: incompatible dimension for operator Add(T,Basevector,T,Basevector)",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<size_; i++)
    data_[i] = a*data_[i]+ b*vec.data_[i];

  CATCH_CAST
}

template<class TYPE> 
void Vector<TYPE>::Axpy(const TYPE a, const CFSVector &y)
{
  ENTER_IFCN("Vector::Axpy"); 
  
  TRY_CAST
  CONSTREFCAST(y,Vector<TYPE>,vec);
#ifdef CHECK_INDEX
  if (size_ != vec.size_)
    Error("Vector: incompatible dimension for operator Add(TYPE,Basevector)",__FILE__,__LINE__);
#endif
  for (Integer i=0; i<size_; i++)
    data_[i] = a*data_[i]+ vec.data_[i];

  CATCH_CAST
}
 
template<class TYPE> 
void Vector<TYPE>::Inner(const CFSVector &y, TYPE &result) const
{
  ENTER_IFCN("Vector::Inner");

  TRY_CAST
  CONSTREFCAST(y,Vector<TYPE>,vec);
#ifdef CHECK_INDEX
  if (size_ != vec.size_)
    Error("Vector: incompatible dimension for operator Add(T,Basevector)",__FILE__,__LINE__);
#endif
  result = 0;
  for (Integer i=0; i<size_; i++)
    result += data_[i] * vec.data_[i];

  CATCH_CAST
}




//*************************************************
//* old interface which is ocmpatible to previous *
//* version of Vector<TYPE>                          *
//*************************************************


template<class TYPE>
void Vector<TYPE>::ToStdVector(std::vector<TYPE> &vec) const
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0) 
    Error("Don't use toStdVector() to undefined vector", __FILE__, __LINE__);
#endif

  vec.resize(size_);
  for (Integer i=0; i<size_; i++)
    vec[i] = data_[i];

}


template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator=(const Vector<TYPE> &x)
{
  
  if (this == &x)
    return *this;
  
  if (size_ != x.size_)
    {	
      if (data_)
	delete [] data_;
      
      size_ = x.size_;
      data_ = new TYPE [size_];
#ifdef CHECK_MEMORY
      if (!data_)
	std::cerr << "Vector: out of memory";
#endif
    }
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] = x.data_[i];
  
  return *this;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator=(const std::vector<TYPE> &x)
{
  
  if (size_ != x.size())
    {	
      if (data_)
	delete [] data_;
      
      size_ = x.size();
      data_ = new TYPE [size_];
#ifdef CHECK_MEMORY
      if (!data_)
	Error("Vector: out of memory",__FILE__,__LINE__);
#endif
    }
  
  for (Integer i = 0; i < size_; i++)
    data_ [i] = x[i];
  
  return *this;
}




template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator+(const Vector<TYPE> &x) const
{	
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator +(vector)",__FILE__,__LINE__);
#endif  
  
#ifdef CHECK_INDEX
  if (size_ != x.size_)
    Error("Vector: incompatible dimension for operator +(vector)",__FILE__,__LINE__);
#endif
  
  Vector ret(size_);
  
  for (Integer i = 0; i < size_; i++)
    ret.data_[i] = data_[i] + x.data_ [i];
  
  return ret;
}


template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator+=(const Vector<TYPE> &x)
{
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator +=(vector)",__FILE__,__LINE__);
#endif  
  
#ifdef CHECK_INDEX
  if (size_ != x.size_)
    Error("Vector: incompatible dimension for operator +=(vector)",__FILE__,__LINE__);
#endif

  for (Integer i = 0; i < size_; i++)
    data_[i] += x.data_[i];
  
  return *this; 
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator- () const
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in oprator -()",__FILE__,__LINE__); 
#endif

  Vector ret(size_);
  
  for (Integer i = 0; i < size_; i++)
    ret.data_ [i] = -data_ [i];
  
  return ret;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator-(const Vector<TYPE> &x) const
{
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator -(vector)",__FILE__,__LINE__);
#endif  
  
#ifdef CHECK_INDEX
  if (size_ != x.size_)
    Error("Vector: incompatible dimension for operator -(vector)",__FILE__,__LINE__);
#endif
  Vector ret(size_);

  for (Integer i = 0; i < size_; i++)
    ret.data_ [i] = data_[i] - x.data_ [i];
  
  return ret;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator-=(const Vector<TYPE> &x)
{
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator -=(vector)",__FILE__,__LINE__);
#endif  
  
#ifdef CHECK_INDEX
  if (size_ != x.size_)
    Error("Vector: incompatible dimension for operator -=(vector)",__FILE__,__LINE__);
#endif

  for (Integer i = 0; i < size_; i++)
    data_[i] -= x.data_[i];
  
  return *this;
}


template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator* (const TYPE &x) const
{	
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in operator *(number)",__FILE__,__LINE__); 
#endif
  
  Vector ret(size_);
  
  for (Integer i = 0; i < size_; i++)
    ret.data_[i] = data_[i] * x;
  
  return ret;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator/ (const TYPE &x) const
{       
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in operator /(number)",__FILE__,__LINE__); 
#endif
  
  Vector ret(size_);
  
  for (Integer i = 0; i < size_; i++)
    ret.data_[i] = data_[i] / x;
  
  return ret;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator/= (const TYPE &x)
{	
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in operator /=(number)",__FILE__,__LINE__); 
#endif

  TYPE y = x;
  
  for (Integer i = 0; i < size_; i++)
    data_[i] /= y;
  
  return *this;
}

template<class TYPE>
TYPE Vector<TYPE>::operator* (const Vector<TYPE> &x) const
{	
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator *(vector)",__FILE__,__LINE__);
#endif  
  
#ifdef CHECK_INDEX
  if (size_ != x.size_)
    Error("Vector: incompatible dimension for operator *(vector)",__FILE__,__LINE__);
#endif

  TYPE ret;
 
  ret = data_[0] * x.data_[0];
  for (Integer i = 1; i < size_; i++)
    ret += data_[i] * x.data_[i];
  
  return ret;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator*(const Matrix<TYPE> &x) const
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in operator *(Matrix)",__FILE__,__LINE__);
  if (!x.data_)
    Error("Vector: undefined Matrix in operator *(Matrix)",__FILE__,__LINE__);
#endif  
  
#ifdef CHECK_INDEX
  if (size_ != x.size_row_)
    Error("Vector: incompatible dimension for operator *(Matrix)",__FILE__,__LINE__);
#endif

  TYPE  a;
  Vector ret(x.size_col_);
  
  for (Integer i = 0; i < x.size_col_; i++)
    {	
      a = data_ [0] * x.data_ [0] [i];
      for (Integer j = 1; j < size_; j++)
	a += data_ [j] * x.data_ [j] [i];
      ret.data_ [i] = a;
    }
  
  return ret;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator=(const Matrix<TYPE> &x) const
{	
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error( "undefined Vector in operator =", __FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (!x.data_)
    Error( "undefined Matrix in operator = ", __FILE__,__LINE__);

  if (x.size_col_ != 1)
    Error( "matrix has more tha one row. No assignment to vector possible ", __FILE__,__LINE__);
#endif
  
  Vector ret(x.size_row_);
  
  for (Integer i = 0; i < x.size_col_; i++)
    ret.data_[i] = x[i][0];

  return ret;
}


template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator*= (const TYPE &x)
{	
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in operator*=",__FILE__,__LINE__); 
#endif
  
  TYPE y = x;
  
  for (Integer i = 0; i < size_; i++)
    data_[i] *= y;
  
  return *this;
}

template<class TYPE>

Vector<TYPE> &Vector<TYPE>::operator*=(const Matrix<TYPE> &x)
{	
  *this = *this * x;
  
  return *this;
}

template<class TYPE>
Integer Vector<TYPE>::operator== (const Vector<TYPE> &x) const
{	
#ifdef CHECK_INITIALIZED 
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator==",__FILE__,__LINE__);
#endif
  
  for (Integer i = 0; i < size_; i++)
    if (data_[i] != x.data_ [i])
      return FALSE;
  
  return TRUE;
}

template<class TYPE>
Integer Vector<TYPE>::operator!= (const Vector<TYPE> &x) const
{	
#ifdef CHECK_INITIALIZED
  if ((size_ == 0) || (x.size_ == 0))
    Error("Vector: undefined Vector in operator !=", __FILE__, __LINE__);
#endif
  
  for (Integer i = 0; i < size_; i++)
    if (data_[i] != x.data_[i])
      return FALSE;
  
  return TRUE;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::Part(const Integer i1,const Integer i2) const
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    std::cerr << "Vector: undefined Vector";
#endif
  
#ifdef CHECK_INDEX
  if (i1 < 0 || i1 > i2 || i2 >= size_)
    std::cerr << "Vector: invalid index";
#endif
  
  Vector ret (i2 - i1 + 1);
  
  for (Integer i = 0; i < ret.size_; i++)
    ret.data_[i] = data_[i1 + i];
  
  return ret;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::Unit (const Integer n,const Integer i)
{
#ifdef CHECK_INDEX	
  if (n <= 0)
    std::cerr << "Vector::unit() invalid dimension";
  
  if (i < 0 || i >= n)
    std::cerr << "Vector::unit() invalid index";
#endif
  
  Vector<TYPE> ret(n);
  
  ret.data_[i] = 1;
  
  return ret;
}


template<class TYPE>
Double Vector<TYPE>::NormL2() const
{	
#ifdef CHECK_INITIALIZED
if (size_ == 0)
  Error("Vector: undefined Vector in function norm_2()", __FILE__, __LINE__);
#endif
	
 TYPE ret = (*this)*(*this);
 return sqrt (ret);
}

Double Vector<Complex>::NormL2() const
{
  Complex temp = (*this)*(*this);
  Double ret = temp.real()*temp.real();
  ret += temp.imag()*temp.imag();
  
  return sqrt(ret);
  
 
}


template<class TYPE>
void Vector<TYPE>::Push_back(const TYPE & y)
{
  // add y at the end of the vector (size_=length of vec)
  AddElement(y,size_);
}



template<class TYPE>
void  Vector<TYPE>:: AddElement (const TYPE & y, Integer pos)
{
#ifdef CHECK_INDEX
  if (pos < 0 || pos > size_)
    Error("Vector::AddElemen(): Index out of bounds", __FILE__,__LINE__);
#endif
  
  Integer i;
  
  TYPE * help=new TYPE[size_+1];
  
  for (i=0; i < pos; i++) 
    help[i]=data_[i];
  
  help[pos]=y;
  for (i=pos+1; i<size_+1; i++) 
    help[i]=data_[i-1];
  
  delete [] data_;
  data_=help;
  size_++;
}


template<class TYPE>
void  Vector<TYPE>::InsertVector (const Vector<TYPE> & y, Integer pos)
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in function InsertVector()", __FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (pos < 0)
    Error("Vector: index is smaller than zero in function InsertVector()", __FILE__,__LINE__);
#endif

   Integer i;
   
   Integer l=y.size_;
   TYPE * help=new TYPE[size_+l];

   for (i=0; i < pos; i++) 
     help[i]=data_[i];
   
   for (i=pos; i < pos+l; i++) 
     help[i]=y.data_[i-pos];

   for (i=pos+l; i < size_+l; i++) 
     help[i]=data_[i-l];
  
   delete [] data_; 
   data_=help;
   size_+=l;
}

template<class TYPE>
void  Vector<TYPE>:: Cut (const Integer pos)
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in function Cut()", __FILE__,__LINE__);
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
}

template<class TYPE>
void  Vector<TYPE>:: Cut (const Integer pos1, const Integer pos2)
{
#ifdef CHECK_INITIALIZED
  if (size_ == 0)
    Error("Vector: undefined Vector in function Cut()", __FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
   if (pos1 < 0 || pos1 >= size_ || pos2 < 0 || pos2 >= size_) 
     Error("Invalid index for cut");
   if (pos1 > pos2)
     Error("First index is bigger than second one in function Cut()",__FILE__,__LINE__);
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
Integer  Vector<TYPE>:: Memory() const
{
 return size_*sizeof(TYPE);
}

template<class TYPE>
void  Vector<TYPE>::TransformInVector(const Integer nsize, TYPE * ptdata) 
{
  size_ = nsize;
#ifdef CHECK_MEMORY
  if (ptdata == 0)
    Error("Pointer to data is undefined in function TransformInVector()",__FILE__,__LINE__);
#endif

  data_ = ptdata;
}

template<class TYPE>
void Swap(Vector<TYPE> & a, Vector<TYPE> & b) 
{
 Integer tmp_size;

 tmp_size = a.size_;
 a.size_ = b.size_;
 b.size_ = tmp_size;

 TYPE * tmp_data_;

 tmp_data = a.data_;
 a.data_ = b.data_;
 b.data_ = temp_data;
}

template<class TYPE> 
void Swap(TYPE & a, TYPE & b)
{
   TYPE tmp=a;
   a=b;
   b=tmp;
}

template void Swap<Integer>(Integer & ,Integer &);
template void Swap<Double>(Double & ,Double &);

template<class S> void Sort(S* v, Integer n)
{
  for (Integer gap=n/2; gap >0; gap/=2)
    for (Integer i=gap; i<n; i++)
       for (Integer j=i-gap; j>=0; j-=gap)
          if (v[j+gap]<v[j]) Swap<S>(v[j],v[j+gap]);

}
template void Sort<Integer>(Integer * ,Integer);
template void Sort<Double>(Double * ,Integer);

// Overloading << for Vector

template<class S>
std::ostream & operator << ( std::ostream & out, const Vector<S> & vc)
{
#ifdef CHECK_MEMORY
  if (vc.GetSize()==0) out << "vector is undefined" ;
#endif

for (Integer i=0; i < vc.GetSize(); i++)
  out << vc[i] << " " << std::endl;
 return out;
}

template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);
template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);
}
