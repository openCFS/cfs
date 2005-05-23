#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

#include "vector.hh"
#include "Matrix/matrix.hh"

namespace CoupledField
{ 

  template<class TYPE>
  Vector<TYPE>::Vector()
    : CFSVector()
  {
    ENTER_IFCN("Vector::Vector");
    data_ = NULL;
    size_ = 0;
    capacity_ = 0;
    memBelongsToMe_ = TRUE;
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(Integer size, const TYPE entry)
    : CFSVector(size)
  {
    ENTER_IFCN("Vector::Vector");
    size_ = size;
    data_ = new TYPE [size];
    memBelongsToMe_ = TRUE;

    for (Integer i = 0; i < size; i++)
      data_ [i] = entry;
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(const Vector<TYPE> &vec)
  {
    ENTER_IFCN("Vector::Vector");    

    size_ = vec.size_;
    data_ = new TYPE [vec.size_];
    memBelongsToMe_ = TRUE;
  
    for (Integer i = 0; i < size_; i++)
      data_ [i] =  vec.data_[i];
  }

  template<class TYPE>
  Vector<TYPE>::~Vector()
  {
    ENTER_IFCN("Vector::~Vector");
    if (data_ && memBelongsToMe_ == TRUE )
      delete[] data_;
  }


  template<class TYPE>
  void Vector<TYPE>::Clear()
  {
    ENTER_IFCN( "Vector::Clear()" );
  
    if (memBelongsToMe_ == FALSE ) {
      (*error) << "Refusing to clear vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
    
    size_ = 0;
    if (data_)
      delete[] data_;
    data_ = NULL;
  }

  template<class TYPE>
  Boolean Vector<TYPE>::IsComplex()
  {
    return FALSE;
  }

  template<>
  Boolean Vector<Complex>::IsComplex()
  {
    return TRUE;
  }

  template<class TYPE>
  void Vector<TYPE>::Init(const TYPE entry)
  {
    ENTER_IFCN("Vector::Init");

    //#ifdef CHECK_INITIALIZED
    //  if (size_ == 0) 
    //    Warning("Don't use Init() to undefined vector", __FILE__, __LINE__);
    //#endif

    for (Integer i=0; i<size_; i++) 
      data_[i]=entry;
  }

  template<class TYPE>
  void Vector<TYPE>::Replace( Integer length, TYPE* entries, Boolean transferMem ) {
  
    ENTER_FCN( "Vector::Replace" );

    // De-allocate old array, if required
    if ( memBelongsToMe_ == TRUE ) {
      delete[] data_ ;
    }

    // Re-set internal attributes
    data_           = entries;
    size_           = length;
    capacity_       = length;
    memBelongsToMe_ = transferMem;  
  }



  template<class TYPE> 
  void Vector<TYPE>::Resize(const Integer size)
  {
    ENTER_IFCN("Vector::Resize");

#ifdef CHECK_INDEX
    if (size <= 0) 
      Error("invalid dimension for Resize", __FILE__, __LINE__);
#endif  
  
    if (memBelongsToMe_ == FALSE ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }

    if (size != size_)
      {
        if (data_) delete[] data_;
      
        size_=size;
        data_ = new TYPE[size_];
      }
  
    for (Integer i = 0; i < size_; i++)
      data_ [i] = TYPE();
  }


  template<class TYPE> 
  void Vector<TYPE>::SetEntry(const Integer i, const TYPE &s)
  {
    ENTER_IFCN("Vector::SetEntry");

#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__, __LINE__);
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
    errorMsgStream << "Vector: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__, __LINE__);
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
    errorMsgStream << "Vector: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Error(errorMsg.c_str(),__FILE__, __LINE__);
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
      Error(errorMsg.c_str(),__FILE__, __LINE__);
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
      Error(errorMsg.c_str(),__FILE__, __LINE__);
#endif

    data_[i]=a*data_[i] + s;
  }

  template<class TYPE> 
  void Vector<TYPE>::Add(const CFSVector& y)
  {
    ENTER_IFCN("Vector::Add");

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>& >(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Error("Vector: incompatible dimension for operator Add(Basevector)",
            __FILE__, __LINE__);
#endif

    for (Integer i=0; i<size_; i++)
      data_[i]+=vec.data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::Add(const TYPE a, const CFSVector &y)
  {
    ENTER_IFCN("Vector::Add"); 

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Error("Vector: incompatible dimension for operator Add(TYPE,Basevector)",
            __FILE__, __LINE__);
#endif

    for (Integer i=0; i<size_; i++)
      data_[i]+=a*vec.data_[i];

  }

  template<class TYPE> 
  void Vector<TYPE>::Add( const TYPE a, const CFSVector& y,
                          const TYPE b, const CFSVector& z )
  {
    ENTER_IFCN("Vector::Add");  

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Error("Vector: incompatible dimension for operator \
Add(T,Basevector,T,Basevector)",__FILE__, __LINE__);
#endif

    for (Integer i=0; i<size_; i++)
      data_[i] = a*data_[i]+ b*vec.data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::Axpy(const TYPE a, const CFSVector &y)
  {
    ENTER_IFCN("Vector::Axpy"); 
    
    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Error("Vector: incompatible dimension for operator Add(TYPE,Basevector)",
            __FILE__, __LINE__);
#endif

    for (Integer i=0; i<size_; i++)
      data_[i] = a*data_[i]+ vec.data_[i];
  }
 
  template<class TYPE> 
  void Vector<TYPE>::Inner(const CFSVector &y, TYPE &result) const
  {
    ENTER_IFCN("Vector::Inner");

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>& >(y);
#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Error("Vector: incompatible dimension for operator Add(T,Basevector)",
            __FILE__, __LINE__);
#endif

    result = 0;
    for (Integer i=0; i<size_; i++)
      result += data_[i] * vec.data_[i];
  }




  //*************************************************
  //* old interface which is ocmpatible to previous *
  //* version of Vector<TYPE>                          *
  //*************************************************


  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator=(const Vector<TYPE> &x)
  {
    ENTER_IFCN( "Vector::operator=(const Vector)");
  
    if (this == &x)
      return *this;
  
    if (size_ != x.size_)
      {   
        if (memBelongsToMe_ == FALSE ) {
          (*error) << "Refusing to resize vector, since memory does not " 
                   << "belong to me!";
          Error( __FILE__, __LINE__ );
        }
        if (data_)
          delete [] data_;
      
        size_ = x.size_;
        data_ = new TYPE [size_];
      }
  
    for (Integer i = 0; i < size_; i++)
      data_ [i] = x.data_[i];
  
    return *this;
  }

  template<class TYPE>
  CFSVector & Vector<TYPE>::operator= (const CFSVector & vec)
  {
    ENTER_IFCN( "Vector::operator=(const CFSVector)");

    Vector<TYPE> const & temp = dynamic_cast<const Vector<TYPE>&>(vec);
  
    if (this == &temp)
      return *this;
  
    if (size_ != temp.size_)
      {   
        if (memBelongsToMe_ == FALSE ) {
          (*error) << "Refusing to resize vector, since memory does not " 
                   << "belong to me!";
          Error( __FILE__, __LINE__ );
        }
      
        if (data_)
          delete [] data_;
      
        size_ = temp.size_;
        data_ = new TYPE [size_];
      }
  
    for (Integer i = 0; i < size_; i++)
      data_ [i] = temp.data_[i];
  
    return dynamic_cast<CFSVector &>(*this);

  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator+(const Vector<TYPE> &x) const
  {       
    ENTER_IFCN( "Vector::operator+" );

#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator +(vector)",
              __FILE__, __LINE__);
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      Error("Vector: incompatible dimension for operator +(vector)",
            __FILE__, __LINE__);
#endif
  
    Vector ret(size_);
  
    for (Integer i = 0; i < size_; i++)
      ret.data_[i] = data_[i] + x.data_ [i];
  
    return ret;
  }


  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator+=(const Vector<TYPE> &x)
  {
    ENTER_IFCN( "Vector::operator+" );
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator +=(vector)",
              __FILE__, __LINE__);
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      Error("Vector: incompatible dimension for operator +=(vector)",
            __FILE__, __LINE__);
#endif

    for (Integer i = 0; i < size_; i++)
      data_[i] += x.data_[i];
  
    return *this; 
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator- () const
  {
    ENTER_IFCN( "Vector::operator-" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in oprator -()",__FILE__, __LINE__); 
#endif

    Vector ret(size_);
  
    for (Integer i = 0; i < size_; i++)
      ret.data_ [i] = -data_ [i];
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator-(const Vector<TYPE> &x) const
  {
    ENTER_IFCN( "Vector::operator-" );
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator -(vector)",
              __FILE__, __LINE__);
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      Error("Vector: incompatible dimension for operator -(vector)",
            __FILE__, __LINE__);
#endif
    Vector ret(size_);

    for (Integer i = 0; i < size_; i++)
      ret.data_ [i] = data_[i] - x.data_ [i];
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator-=(const Vector<TYPE> &x)
  {
    ENTER_IFCN( "Vector::operator-=" );
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator -=(vector)",
              __FILE__, __LINE__);
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      Error("Vector: incompatible dimension for operator -=(vector)",
            __FILE__, __LINE__);
#endif

    for (Integer i = 0; i < size_; i++)
      data_[i] -= x.data_[i];
  
    return *this;
  }


  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator* (const TYPE &x) const
  {
    ENTER_IFCN( "Vector::operator*" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in operator *(number)",
              __FILE__, __LINE__); 
#endif
  
    Vector ret(size_);
  
    for (Integer i = 0; i < size_; i++)
      ret.data_[i] = data_[i] * x;
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator/ (const TYPE &x) const
  {
    ENTER_IFCN( "Vector::operator/" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in operator /(number)",
              __FILE__, __LINE__); 
#endif
  
    Vector ret(size_);
  
    for (Integer i = 0; i < size_; i++)
      ret.data_[i] = data_[i] / x;
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator/= (const TYPE &x)
  {
    ENTER_IFCN( "Vector::operator/=" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in operator /=(number)",
              __FILE__, __LINE__); 
#endif

    TYPE y = x;
  
    for (Integer i = 0; i < size_; i++)
      data_[i] /= y;
  
    return *this;
  }

  template<class TYPE>
  TYPE Vector<TYPE>::operator* (const Vector<TYPE> &x) const
  {
    ENTER_IFCN( "Vector::operator*" );
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator *(vector)",
              __FILE__, __LINE__);
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      Error("Vector: incompatible dimension for operator *(vector)", 
            __FILE__, __LINE__);
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
    ENTER_IFCN( "Vector::operator*" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in operator *(Matrix)",
              __FILE__, __LINE__);
    if (!x.data_)
      Warning("Vector: undefined Matrix in operator *(Matrix)",
              __FILE__, __LINE__);
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_row_)
      Error("Vector: incompatible dimension for operator *(Matrix)",
            __FILE__, __LINE__);
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
    ENTER_IFCN( "Vector::operator=(const Matrix)" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning( "undefined Vector in operator =", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
    if (!x.data_)
      Error( "undefined Matrix in operator = ", __FILE__, __LINE__);

    if (x.size_col_ != 1)
      Error( "matrix has more tha one row. No assignment to vector possible ", 
             __FILE__, __LINE__);
#endif
  
    Vector ret(x.size_row_);
  
    for (Integer i = 0; i < x.size_col_; i++)
      ret.data_[i] = x[i][0];

    return ret;
  }


  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator*= (const TYPE &x)
  {
    ENTER_IFCN( "Vector::operator*=" );   
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in operator*=",__FILE__, __LINE__); 
#endif
  
    TYPE y = x;
  
    for (Integer i = 0; i < size_; i++)
      data_[i] *= y;
  
    return *this;
  }

  template<class TYPE>

  Vector<TYPE> &Vector<TYPE>::operator*=(const Matrix<TYPE> &x)
  {       
    ENTER_IFCN( "Vector::operator*=" );
    *this = *this * x;
  
    return *this;
  }

  template<class TYPE>
  Integer Vector<TYPE>::operator== (const Vector<TYPE> &x) const
  {       
#ifdef CHECK_INITIALIZED 
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator==",
              __FILE__, __LINE__);
#endif
  
    for (Integer i = 0; i < size_; i++)
      if (data_[i] != x.data_ [i])
        return FALSE;
  
    return TRUE;
  }

  template<class TYPE>
  Integer Vector<TYPE>::operator!= (const Vector<TYPE> &x) const
  {
    ENTER_IFCN( "Vector::operator!=" );
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator !=", 
              __FILE__, __LINE__);
#endif
  
    for (Integer i = 0; i < size_; i++)
      if (data_[i] != x.data_[i])
        return FALSE;
  
    return TRUE;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::Part(const Integer i1,const Integer i2) const
  {
    ENTER_IFCN( "Vector::Part" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning( "Vector: undefined Vector", __FILE__, __LINE__);
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
    ENTER_IFCN( "Vector::Unit" );
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
    ENTER_IFCN( "Vector::NormL2" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function norm_2()",  
              __FILE__, __LINE__);
#endif
        
    TYPE ret = (*this)*(*this);
    return sqrt ((Double)ret);
  }

  template<>
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
    ENTER_IFCN( "Vector::Push_back" );
    // add y at the end of the vector (size_=length of vec)
    AddElement(y,size_);
  }



  template<class TYPE>
  void  Vector<TYPE>:: AddElement (const TYPE & y, Integer pos)
  {
    ENTER_IFCN( "Vector::AddElement" );
#ifdef CHECK_INDEX
    if (pos < 0 || pos > size_)
      Error("Vector::AddElemen(): Index out of bounds", __FILE__, __LINE__);
#endif
  
    if (memBelongsToMe_ == FALSE ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
  
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
    ENTER_IFCN( "Vector::InsertVector" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function InsertVector()", 
              __FILE__, __LINE__);
#endif
  
#ifdef CHECK_INDEX
    if (pos < 0)
      Error("Vector: index is smaller than zero in function InsertVector()", 
            __FILE__, __LINE__);
#endif
  
    if (memBelongsToMe_ == FALSE ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
  
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
    ENTER_IFCN( "Vector::Cut" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function Cut()", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
    if (pos<0 || pos >=size_) 
      Error("Invalid index for cut");
#endif
   
    if (memBelongsToMe_ == FALSE ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
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
    ENTER_IFCN( "Vector::Cut" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function Cut()", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
    if (pos1 < 0 || pos1 >= size_ || pos2 < 0 || pos2 >= size_) 
      Error("Invalid index for cut");
    if (pos1 > pos2)
      Error("First index is bigger than second one in function Cut()",
            __FILE__, __LINE__);
#endif
 
    if (memBelongsToMe_ == FALSE ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
   
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
  Integer  Vector<TYPE>::Memory() const
  {
    ENTER_IFCN( "Vector::Memory" );
    return size_*sizeof(TYPE);
  }

  template<class TYPE>
  void Swap(Vector<TYPE> & a, Vector<TYPE> & b) 
  {
    ENTER_IFCN( "Swap" );
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
  void Vector<TYPE>::Sort()
  {
    ENTER_IFCN( "Vector::Sort" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector::Sort: Use of uninitialized vector");
#endif
    TYPE temp; 
 
    for (Integer gap=size_/2; gap >0; gap/=2)
      for (Integer i=gap; i<size_; i++)
        for (Integer j=i-gap; j>=0; j-=gap)
          if (data_[j+gap] < data_[j]) {
            temp = data_[j];
            data_[j] = data_[j+gap];
            data_[j+gap] = temp;
          }
  }

  template<> 
  void Vector<Complex>::Sort(){
    Error("Vector<Complex>::Sort: Not defined for complex numbers",
          __FILE__, __LINE__);
  }

  // Overloading << for Vector
  template<class S>
  std::ostream & operator << ( std::ostream & out, const Vector<S> & vc)
  {
    ENTER_IFCN( "operator <<(Vector)" );

    for (Integer i=0; i < vc.GetSize(); i++)
      out << vc[i] << " " << std::endl;
    return out;
  }

#ifdef __GNUC__
  template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);
  template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);
  template std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> &);
#endif


  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Vector<Integer>
#pragma instantiate Vector<Double>
#pragma instantiate Vector<Complex>
#pragma instantiate std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &)
#pragma instantiate std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &)
#pragma instantiate std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> &)
#endif
}
