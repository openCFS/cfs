// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

#include "vector.hh"
#include "Matrix/matrix.hh"
#include "Utils/tools.hh"

namespace CoupledField
{ 

  template<class TYPE>
  Vector<TYPE>::Vector()
    : CFSVector()
  {
    data_ = NULL;
    size_ = 0;
    capacity_ = 0;
    memBelongsToMe_ = true;
  }

  template<class TYPE> 
  Vector<TYPE>::Vector(UInt size, const TYPE entry)
    : CFSVector(size)
  {
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
  Vector<TYPE>::Vector(const Point & p)
  {
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
    if (data_ && memBelongsToMe_ == true )
      delete[] data_;
  }


  template<class TYPE>
  void Vector<TYPE>::Clear()
  {
  
    if (memBelongsToMe_ == false ) {
      EXCEPTION( "Refusing to clear vector, since memory does not " 
                 << "belong to me!" );
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

    //#ifdef CHECK_INITIALIZED
    //  if (size_ == 0) 
    //    EXCEPTION("Don't use Init() to undefined vector" );
    //#endif

    for (UInt i=0; i<size_; i++) 
      data_[i]=entry;
  }

  template<class TYPE>
  void Vector<TYPE>::Replace( UInt length, TYPE* entries, bool transferMem ) {
  

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

#ifdef CHECK_INDEX
    if (size < 0) 
      EXCEPTION( "invalid dimension for Resize" );
#endif  
  

    if (size != size_)
      {
        if (memBelongsToMe_ == false ) {
          EXCEPTION( "Refusing to resize vector, since memory does not " 
                     << "belong to me!" );
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

#ifdef CHECK_INDEX
    if (i >= size_) {
      EXCEPTION( "Vector: invalid access to element " << i 
                 << "\n Length of vector: " << size_ );
    }
#endif
    data_[i] = s;
  }
  
  template<class TYPE> 
  void Vector<TYPE>::GetEntry(const UInt i, TYPE &ret) const
  {

#ifdef CHECK_INDEX
    if (i >= size_) {
      EXCEPTION( "Vector: invalid access to element "
                 << i << "\n Length of vector: " << size_ );
    }
#endif
    ret=  data_[i];
  }


  template<class TYPE>
  Matrix<Double> Matrix<TYPE>::GetPart(  DataType part ) const {
    EXCEPTION( "Matrix::GetPart: Only Implemented for Real and " 
               << " Complex matrices!" );
    Matrix<Double> temp;
    return temp;
  }

  
  template<>
  Vector<Double> Vector<Double>::GetPart(  DataType part ) const {
    
    if ( part != REAL ) {
      EXCEPTION("Vector<Double>::GetPart: Only possible for REAL part." );
    }
      return *this;
  }

  template<>
  Vector<Double> Vector<Complex>::GetPart(  DataType part ) const {
    
    Vector<Double> ret;
    if ( part == REAL ) {
      ret.Resize( size_ );
      for ( UInt i = 0; i < size_; i++ ) {
	ret[i]  = data_[i].real();
      }
    } 
    else if ( part == IMAG ) {
      ret.Resize( size_ );
      for ( UInt i = 0; i < size_; i++ ) {
	ret[i]  = data_[i].imag();
      }
    } 
    else {
      EXCEPTION("Vector<Complex>::GetPart: Only possible for REAL or IMAG part!" );
    }
    
    return ret;
  }


  template<class TYPE>
  void Vector<TYPE>::SetPart( DataType part, const Vector<Double> & partVector ) {
    EXCEPTION( "Vector::SetPart:" << 
               "Only Implemented for Real and Complex vectors!" );
  }
  
  template<>
    void Vector<Double>::SetPart( DataType part, const Vector<Double> & partVector ) {
    
    if ( size_ != partVector.GetSize() ) {
      EXCEPTION( "Vector<Double>::SetPart: Dimension of vectors do not match!" );
    }
 
    if ( part != REAL ) {
      EXCEPTION( "Vector<Double>::SetPart: Only possible for REAL part." );
    }
    *this = partVector;
  }

  template<>
  void Vector<Complex>::SetPart( DataType part, const Vector<Double> & partVector ) {

    if ( size_ != partVector.GetSize() ) {
      EXCEPTION( "Vector<Complex>::SetPart: Dimension of vectors do not match!" );
    }
        
    if ( part == REAL ) {
      for ( UInt i = 0; i < size_; i++ ) {
	data_[i]  = Complex( partVector[i], data_[i].imag() );
      }
    } 
    else if ( part == IMAG ) {
      for ( UInt i = 0; i < size_; i++ ) {
	data_[i]  = Complex( data_[i].real(), partVector[i] );
      }
    } 
    else {
      EXCEPTION( "Vector<Complex>::SetPart: Only possible for REAL or IMAG part!" );
    }
  }
 

  template<class TYPE> 
  void Vector<TYPE>::AddEntry(const UInt i, const TYPE &s)
  {

#ifdef CHECK_INDEX
    if (i >= size_) {
      EXCEPTION( "Vector: invalid access to element "
                 << i << "\n Length of vector: " << size_ );
    }
#endif
      
      data_[i]+=s;
    }

  template<class TYPE> 
  void Vector<TYPE>::MultEntry(const UInt i, const TYPE &s)
  {

#ifdef CHECK_INDEX
    if (i >= size_) {
      EXCEPTION( "Vector::MultEntry: invalid access to element "
                 << i << "\n Length of vector: " << size_ );
    }
#endif

    data_[i]+=s;
  }

  template<class TYPE> 
  void Vector<TYPE>::MultAddEntry(const UInt i, const TYPE &a, const TYPE &s)
  {

#ifdef CHECK_INDEX
    if (i >= size_) {
      EXCEPTION( "Vector::MultEntry: invalid access to element "
                 << i << "\n Length of vector: " << size_ );
    }
#endif

    data_[i]=a*data_[i] + s;
  }

  template<class TYPE> 
  void Vector<TYPE>::Add(const CFSVector& y)
  {

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>& >(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      EXCEPTION( "Vector: incompatible dimension for operator Add(Basevector)" );
#endif

    for (UInt i=0; i<size_; i++)
      data_[i]+=vec.data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::Add(const TYPE a, const CFSVector &y)
  {

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      EXCEPTION("Vector: incompatible dimension for operator Add(TYPE,Basevector)" );
#endif

    for (UInt i=0; i<size_; i++)
      data_[i]+=a*vec.data_[i];

  }

  template<class TYPE> 
  void Vector<TYPE>::Add( const TYPE a, const CFSVector& y,
                          const TYPE b, const CFSVector& z )
  {

    const Vector<TYPE> & yvec = dynamic_cast<const Vector<TYPE>&>(y);
    const Vector<TYPE> & zvec = dynamic_cast<const Vector<TYPE>&>(z);

#ifdef CHECK_INDEX
    if (size_ != yvec.size_ || size_ != zvec.size_)
      EXCEPTION("Vector: incompatible dimension for operator"
                << "Add(T,Basevector,T,Basevector)" );
#endif

    for (UInt i=0; i<size_; i++)
      data_[i] = a*yvec.data_[i]+ b*zvec.data_[i];
  }

  template<class TYPE> 
  void Vector<TYPE>::Axpy(const TYPE a, const CFSVector &y)
  {
    
    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>&>(y);

#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      EXCEPTION( "Vector: incompatible dimension for operator Add(TYPE,Basevector)" );
#endif

    for (UInt i=0; i<size_; i++)
      data_[i] = a*data_[i]+ vec.data_[i];
  }
 
  template<class TYPE> 
  void Vector<TYPE>::Inner(const CFSVector &y, TYPE &result) const
  {

    const Vector<TYPE> & vec = dynamic_cast<const Vector<TYPE>& >(y);
#ifdef CHECK_INDEX
    if (size_ != vec.size_)
      EXCEPTION("Vector: incompatible dimension for operator Add(T,Basevector)" );
#endif

    result = 0;
    for (UInt i=0; i<size_; i++)
      result += data_[i] * vec.data_[i];
  }

#ifndef EXPR_TEMPLATES


  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator=(const Vector<TYPE> &x)
  {
  
    if (this == &x)
      return *this;
  
    if (size_ != x.size_)
      {   
        if (memBelongsToMe_ == false ) {
          EXCEPTION( "Refusing to resize vector, since memory does not " 
                     << "belong to me!" );
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
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      EXCEPTION("Vector: undefined Vector in operator +=(vector)" );
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      EXCEPTION("Vector: incompatible dimension for operator +=(vector)" );
#endif

    for (UInt i = 0; i < size_; i++)
      data_[i] += x.data_[i];
  
    return *this; 
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator- () const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in oprator -()" );
#endif

    Vector ret(size_);
  
    for (UInt i = 0; i < size_; i++)
      ret.data_ [i] = -data_ [i];
  
    return ret;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::operator+ () const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in oprator +()" );
#endif
    
    return *this;
    
  }

 

  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator-=(const Vector<TYPE> &x)
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      EXCEPTION("Vector: undefined Vector in operator -=(vector)" );
#endif  
  
#ifdef CHECK_INDEX
    if (size_ != x.size_)
      EXCEPTION("Vector: incompatible dimension for operator -=(vector)" );
#endif

    for (UInt i = 0; i < size_; i++)
      data_[i] -= x.data_[i];
  
    return *this;
  }


 

  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator/= (const TYPE &x)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator /=(number)" );
#endif

    TYPE y = x;
  
    for (UInt i = 0; i < size_; i++)
      data_[i] /= y;
  
    return *this;
  }


  template<class TYPE>
  Vector<TYPE> &Vector<TYPE>::operator*= (const TYPE &x)
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in operator*=" );
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
      EXCEPTION("Vector: undefined Vector in operator==" );
#endif
  
    for (UInt i = 0; i < size_; i++)
      if (data_[i] != x.data_ [i])
        return false;
  
    return true;
  }

  template<class TYPE>
  bool Vector<TYPE>::operator!= (const Vector<TYPE> &x) const
  {
#ifdef CHECK_INITIALIZED
    if ((size_ == 0) || (x.size_ == 0))
      EXCEPTION("Vector: undefined Vector in operator !=" );
#endif
  
    for (UInt i = 0; i < size_; i++)
      if (data_[i] != x.data_[i])
        return false;
  
    return true;
  }

  template<class TYPE>
  Vector<TYPE> Vector<TYPE>::Part(const UInt i1,const UInt i2) const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION( "Vector: undefined Vector" );
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
    EXCEPTION( "Vector<TPYE>::NormL2 only defined for TYPE=Complex/Double" );
    return TYPE();
  }

 template<>
  Double Vector<Double>::NormL2() const
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in function norm_2()" );
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
  std::string Vector<TYPE>::Serialize( Char separator ) const {
    std::stringstream out;
    
    if( size_ > 0 ) {
      for( UInt i = 0; i < size_-1; i++ ) {
        out << data_[i] << separator << " ";
      }
      out << data_[size_-1];
    }
    return out.str();
  }


  template<class TYPE>
  void Vector<TYPE>::Push_back(const TYPE & y)
  {
    // add y at the end of the vector (size_=length of vec)
    AddElement(y,size_);
  }



  template<class TYPE>
  void  Vector<TYPE>:: AddElement (const TYPE & y, UInt pos)
  {
#ifdef CHECK_INDEX
    if (pos < 0 || pos > size_)
      EXCEPTION("Vector::AddElemen(): Index out of bounds" );
#endif
  
    if (memBelongsToMe_ == false ) {
      EXCEPTION( "Refusing to resize vector, since memory does not " 
                 << "belong to me!" );
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
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in function InsertVector()" );
#endif
  
#ifdef CHECK_INDEX
    if (pos < 0)
      EXCEPTION("Vector: index is smaller than zero in function InsertVector()" );
#endif
  
    if (memBelongsToMe_ == false ) {
      EXCEPTION( "Refusing to resize vector, since memory does not " 
                 << "belong to me!" );
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
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in function Cut()" );
#endif

#ifdef CHECK_INDEX
    if (pos<0 || pos >=size_) 
      EXCEPTION("Invalid index for cut");
#endif
   
    if (memBelongsToMe_ == false ) {
      EXCEPTION( "Refusing to resize vector, since memory does not " 
                 << "belong to me!" );
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
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector: undefined Vector in function Cut()" );
#endif

#ifdef CHECK_INDEX
    if (pos1 < 0 || pos1 >= size_ || pos2 < 0 || pos2 >= size_) 
      EXCEPTION("Invalid index for cut");
    if (pos1 > pos2)
      EXCEPTION("First index is bigger than second one in function Cut()" );
#endif
 
    if (memBelongsToMe_ == false ) {
      EXCEPTION( "Refusing to resize vector, since memory does not " 
                 << "belong to me!" );
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
    return size_*sizeof(TYPE);
  }

  template<class TYPE>
  void Swap(Vector<TYPE> & a, Vector<TYPE> & b) 
  {
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
  std::string Vector<TYPE>::ToString()
  {
     std::ostringstream os;
     os << "[";
     for(unsigned int i=0; i<size_; i++) {
        os << data_[i];
        if(i < size_-1) os << ", ";
     }
     os << "]";
     return os.str();
  }
     

  template<class TYPE> 
  void Vector<TYPE>::Sort()
  {
#ifdef CHECK_INITIALIZED
    if (size_ == 0)
      EXCEPTION("Vector::Sort: Use of uninitialized vector");
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
    EXCEPTION( "Vector<Complex>::Sort: Not defined for complex numbers" );
  }

  // Overloading << for Vector
  template<class S>
  std::ostream & operator << ( std::ostream & out, const Vector<S> & vc)
  {

    for (UInt i=0; i < vc.GetSize(); i++)
      out << vc[i] << " " << std::endl;
    return out;
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

#include "Utils/boost-serialization.hh"
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Double>, "CoupledField_Vector_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Complex>, "CoupledField_Vector_Complex")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Integer>, "CoupledField_Vector_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::UInt>, "CoupledField_Vector_UInt")


