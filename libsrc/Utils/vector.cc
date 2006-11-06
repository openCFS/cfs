#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

#include "vector.hh"
#include "Matrix/matrix.hh"
#include "Utils/tools.hh"

#include "boost-serialization.hh"

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
    memBelongsToMe_ = true;
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(UInt size, const TYPE entry)
    : CFSVector(size)
  {
    ENTER_IFCN("Vector::Vector");
    size_ = size;
    capacity_ = size;
    data_ = new TYPE [size];
    memBelongsToMe_ = true;

    for (UInt i = 0; i < size; i++)
      data_ [i] = entry;
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(const Vector<TYPE> &vec)
  {
    ENTER_IFCN("Vector::Vector");    

    size_ = vec.size_;
    capacity_ = size_;
    
    if ( size_ > 0 ) {
      data_ = new TYPE [vec.size_];
    } else {
      data_ = NULL;
    }
    
    memBelongsToMe_ = true;
  
    for (UInt i = 0; i < size_; i++)
      data_ [i] =  vec.data_[i];
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(const Point<2> & p)
  {
    ENTER_IFCN("Vector::Vector");
    size_ = 2;
    capacity_ = 2;
    data_ = new TYPE[2];
    memBelongsToMe_ = true;
    data_[0] = (TYPE) p[0];
    data_[1] = (TYPE) p[1];
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(const Point<3> & p)
  {
    ENTER_IFCN("Vector::Vector");
    size_ = 3;
    capacity_ = 3;
    data_ = new TYPE[3];
    memBelongsToMe_ = true;
    data_[0] = (TYPE) p[0];
    data_[1] = (TYPE) p[1];
    data_[2] = (TYPE) p[2];
  }
  

  template<class TYPE>
  Vector<TYPE>::~Vector()
  {
    ENTER_IFCN("Vector::~Vector");
    if (data_ && memBelongsToMe_ == true )
      delete[] data_;
  }


  template<class TYPE>
  void Vector<TYPE>::Clear()
  {
    ENTER_IFCN( "Vector::Clear()" );
  
    if (memBelongsToMe_ == false ) {
      (*error) << "Refusing to clear vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
    
    size_ = 0;
    capacity_ = 0;
    if (data_)
      delete[] data_;
    data_ = NULL;
  }

  template<class TYPE>
  void Vector<TYPE>::Init(const TYPE entry)
  {
    ENTER_IFCN("Vector::Init");

    //#ifdef CHECK_INITIALIZED
    //  if (size_ == 0) 
    //    Warning("Don't use Init() to undefined vector", __FILE__, __LINE__);
    //#endif

    for (UInt i=0; i<size_; i++) 
      data_[i]=entry;
  }

  template<class TYPE>
  void Vector<TYPE>::Replace( UInt length, TYPE* entries, bool transferMem ) {
  
    ENTER_FCN( "Vector::Replace" );

    // De-allocate old array, if required
    if ( memBelongsToMe_ == true ) {
      delete[] data_ ;
    }

    // Re-set internal attributes
    data_           = entries;
    size_           = length;
    capacity_       = length;
    memBelongsToMe_ = transferMem;  
  }



  template<class TYPE> 
  void Vector<TYPE>::Resize(const UInt size)
  {
    ENTER_IFCN("Vector::Resize");

#ifdef CHECK_INDEX
    if (size < 0) 
      Error("invalid dimension for Resize", __FILE__, __LINE__);
#endif  
  

    if (size != size_)
      {
        if (memBelongsToMe_ == false ) {
          (*error) << "Refusing to resize vector, since memory does not " 
                   << "belong to me!";
          Error( __FILE__, __LINE__ );
        }  
        
        if (data_) delete[] data_;
        
        size_=size;
        data_ = new TYPE[size_];
      }
    
 //    for (UInt i = 0; i < size_; i++)
//       data_ [i] = TYPE();
  }


  template<class TYPE> 
  void Vector<TYPE>::SetEntry(const UInt i, const TYPE &s)
  {
    ENTER_IFCN("Vector::SetEntry");

#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream warnMsgStream(errorMsg);
    warnMsgStream << "Vector: invalid access to element ";
    warnMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Warning(errorMsg.c_str(),__FILE__, __LINE__);
#endif
    data_[i] = s;
  }
  
  template<class TYPE> 
  void Vector<TYPE>::GetEntry(const UInt i, TYPE &ret) const
  {
    ENTER_IFCN("Vector::GetEntry");  

#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream warnMsgStream(errorMsg);
    warnMsgStream << "Vector: invalid access to element ";
    warnMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Warning(errorMsg.c_str(),__FILE__, __LINE__);
#endif
    ret=  data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::AddEntry(const UInt i, const TYPE &s)
  {
    ENTER_IFCN("Vector::AddEntry");

#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream warnMsgStream(errorMsg);
    warnMsgStream << "Vector: invalid access to element ";
    warnMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Warning(errorMsg.c_str(),__FILE__, __LINE__);
#endif

    data_[i]+=s;
  }

  template<class TYPE> 
  void Vector<TYPE>::MultEntry(const UInt i, const TYPE &s)
  {
    ENTER_IFCN("Vector::MultEntry");

#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector::MultEntry: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Warning(errorMsg.c_str(),__FILE__, __LINE__);
#endif

    data_[i]+=s;
  }

  template<class TYPE> 
  void Vector<TYPE>::MultAddEntry(const UInt i, const TYPE &a, const TYPE &s)
  {
    ENTER_IFCN("Vector::MuldAddEntry");

#ifdef CHECK_INDEX
    std::string errorMsg;
    std::stringstream errorMsgStream(errorMsg);
    errorMsgStream << "Vector::MultEntry: invalid access to element ";
    errorMsgStream << i << "\n Length of vector: " << size_;
    if (i >= size_)
      Warning(errorMsg.c_str(),__FILE__, __LINE__);
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
      Warning("Vector: incompatible dimension for operator Add(Basevector)",
            __FILE__, __LINE__);
#endif

    for (UInt i=0; i<size_; i++)
      data_[i]+=vec.data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::Add(const TYPE a, const CFSVector &y)
  {
    ENTER_IFCN("Vector::Add"); 

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Warning("Vector: incompatible dimension for operator Add(TYPE,Basevector)",
            __FILE__, __LINE__);
#endif

    for (UInt i=0; i<size_; i++)
      data_[i]+=a*vec.data_[i];

  }

  template<class TYPE> 
  void Vector<TYPE>::Add( const TYPE a, const CFSVector& y,
                          const TYPE b, const CFSVector& z )
  {
    ENTER_IFCN("Vector::Add");  

    const Vector<TYPE> & yvec = dynamic_cast<const Vector<TYPE>&>(y);
    const Vector<TYPE> & zvec = dynamic_cast<const Vector<TYPE>&>(z);

#ifdef CHECK_INDEX
    if (size_ != yvec.size_ || size_ != zvec.size_)
      Warning("Vector: incompatible dimension for operator \
Add(T,Basevector,T,Basevector)",__FILE__, __LINE__);
#endif

    for (UInt i=0; i<size_; i++)
      data_[i] = a*yvec.data_[i]+ b*zvec.data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::Axpy(const TYPE a, const CFSVector &y)
  {
    ENTER_IFCN("Vector::Axpy"); 
    
    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Warning("Vector: incompatible dimension for operator Add(TYPE,Basevector)",
            __FILE__, __LINE__);
#endif

    for (UInt i=0; i<size_; i++)
      data_[i] = a*data_[i]+ vec.data_[i];
  }
 
  template<class TYPE> 
  void Vector<TYPE>::Inner(const CFSVector &y, TYPE &result) const
  {
    ENTER_IFCN("Vector::Inner");

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>& >(y);
#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      Warning("Vector: incompatible dimension for operator Add(T,Basevector)",
            __FILE__, __LINE__);
#endif

    result = 0;
    for (UInt i=0; i<size_; i++)
      result += data_[i] * vec.data_[i];
  }

#ifndef EXPR_TEMPLATES


  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator=(const Vector<TYPE> &x)
  {
    ENTER_IFCN( "Vector::operator=(const Vector)");
  
    if (this == &x)
      return *this;
  
    if (size_ != x.size_)
      {   
        if (memBelongsToMe_ == false ) {
          (*error) << "Refusing to resize vector, since memory does not " 
                   << "belong to me!";
          Error( __FILE__, __LINE__ );
        }
        if (data_)
          delete [] data_;
      
        size_ = x.size_;
        data_ = new TYPE [size_];
      }
  
    for (UInt i = 0; i < size_; i++)
      data_ [i] = x.data_[i];
  
    return *this;
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
      Warning("Vector: incompatible dimension for operator +=(vector)",
            __FILE__, __LINE__);
#endif

    for (UInt i = 0; i < size_; i++)
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
  
    for (UInt i = 0; i < size_; i++)
      ret.data_ [i] = -data_ [i];
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator+ () const
  {
    ENTER_IFCN( "Vector::operator+" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in oprator +()",__FILE__, __LINE__); 
#endif
    
    return *this;
    
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
      Warning("Vector: incompatible dimension for operator -=(vector)",
            __FILE__, __LINE__);
#endif

    for (UInt i = 0; i < size_; i++)
      data_[i] -= x.data_[i];
  
    return *this;
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
  
    for (UInt i = 0; i < size_; i++)
      data_[i] /= y;
  
    return *this;
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
  
    for (UInt i = 0; i < size_; i++)
      data_[i] *= y;
  
    return *this;
  }


#endif

  template<class TYPE>
  bool Vector<TYPE>::operator== (const Vector<TYPE> &x) const
  {       
#ifdef CHECK_INITIALIZED 
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator==",
              __FILE__, __LINE__);
#endif
  
    for (UInt i = 0; i < size_; i++)
      if (data_[i] != x.data_ [i])
        return false;
  
    return true;
  }

  template<class TYPE>
  bool Vector<TYPE>::operator!= (const Vector<TYPE> &x) const
  {
    ENTER_IFCN( "Vector::operator!=" );
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      Warning("Vector: undefined Vector in operator !=", 
              __FILE__, __LINE__);
#endif
  
    for (UInt i = 0; i < size_; i++)
      if (data_[i] != x.data_[i])
        return false;
  
    return true;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::Part(const UInt i1,const UInt i2) const
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
  
    for (UInt i = 0; i < ret.size_; i++)
      ret.data_[i] = data_[i1 + i];
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::Unit (const UInt n,const UInt i)
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
    (*error) << "Vector<TPYE>::NormL2 only defined for TYPE=Complex/Double";
    Error( __FILE__, __LINE__ );
    return TYPE();
  }

 template<>
  Double Vector<Double>::NormL2() const
  {
    ENTER_IFCN( "Vector::NormL2" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function norm_2()",  
              __FILE__, __LINE__);
#endif
        
    Double ret = (*this)*(*this);
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
  void  Vector<TYPE>:: AddElement (const TYPE & y, UInt pos)
  {
    ENTER_IFCN( "Vector::AddElement" );
#ifdef CHECK_INDEX
    if (pos < 0 || pos > size_)
      Warning("Vector::AddElemen(): Index out of bounds", __FILE__, __LINE__);
#endif
  
    if (memBelongsToMe_ == false ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
  
    UInt i;
  
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
  void  Vector<TYPE>::InsertVector (const Vector<TYPE> & y, UInt pos)
  {
    ENTER_IFCN( "Vector::InsertVector" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function InsertVector()", 
              __FILE__, __LINE__);
#endif
  
#ifdef CHECK_INDEX
    if (pos < 0)
      Warning("Vector: index is smaller than zero in function InsertVector()", 
            __FILE__, __LINE__);
#endif
  
    if (memBelongsToMe_ == false ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
  
    UInt i;
  
    UInt l=y.size_;
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
  void  Vector<TYPE>:: Cut (const UInt pos)
  {
    ENTER_IFCN( "Vector::Cut" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function Cut()", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
    if (pos<0 || pos >=size_) 
      Warning("Invalid index for cut");
#endif
   
    if (memBelongsToMe_ == false ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
    UInt i;
   
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
  void  Vector<TYPE>:: Cut (const UInt pos1, const UInt pos2)
  {
    ENTER_IFCN( "Vector::Cut" );
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      Warning("Vector: undefined Vector in function Cut()", __FILE__, __LINE__);
#endif

#ifdef CHECK_INDEX
    if (pos1 < 0 || pos1 >= size_ || pos2 < 0 || pos2 >= size_) 
      Warning("Invalid index for cut");
    if (pos1 > pos2)
      Warning("First index is bigger than second one in function Cut()",
            __FILE__, __LINE__);
#endif
 
    if (memBelongsToMe_ == false ) {
      (*error) << "Refusing to resize vector, since memory does not " 
               << "belong to me!";
      Error( __FILE__, __LINE__ );
    }
   
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
  UInt  Vector<TYPE>::Memory() const
  {
    ENTER_IFCN( "Vector::Memory" );
    return size_*sizeof(TYPE);
  }

  template<class TYPE>
  void Swap(Vector<TYPE> & a, Vector<TYPE> & b) 
  {
    ENTER_IFCN( "Swap" );
    UInt tmp_size;

    tmp_size = a.size_;
    a.size_ = b.size_;
    b.size_ = tmp_size;

    TYPE * tmp_data;

    tmp_data = a.data_;
    a.data_ = b.data_;
    b.data_ = tmp_data;
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
 
    for (UInt gap=size_/2; gap >0; gap/=2)
      for (UInt i=gap; i<size_; i++)
        for (UInt j=i-gap; j>=0; j-=gap)
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

    for (UInt i=0; i < vc.GetSize(); i++)
      out << vc[i] << " " << std::endl;
    return out;
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


// explicit template instantiation for GCC compiler
#if defined(__GNUC__) 
  template class Vector<Integer>;
  template class Vector<Double>;
  template class Vector<UInt>;
  template class Vector<Complex>;
  template class Vector<bool>;
  template std::ostream & operator<<<UInt> (std::ostream & , const Vector<UInt> &);
  template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);
  template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);
  template std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> &);
#endif

  // explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Vector<UInt>
#pragma instantiate Vector<Integer>
#pragma instantiate Vector<Double>
#pragma instantiate Vector<Complex>
#pragme instantiate Vector<bool>;
#pragma instantiate std::ostream & operator<<<UInt> (std::ostream & , const Vector<UInt> &)
#pragma instantiate std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &)
#pragma instantiate std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &)
#pragma instantiate std::ostream & operator<<<Complex> (std::ostream & , const Vector<Complex> &)
#endif
}

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::CFSVector, "CoupledField_CFSVector")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Double>, "CoupledField_Vector_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Complex>, "CoupledField_Vector_Complex")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Integer>, "CoupledField_Vector_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::UInt>, "CoupledField_Vector_UInt")


