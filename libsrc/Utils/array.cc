#include "array.hh"
#include <Matrix/matrix.hh>

namespace CoupledField{

template<class TYPE>
Array<TYPE>::Array()
{
  
  dim_ = 0;
  size_ = 0;
  sol_ = 0;

}


template<class TYPE>
Array<TYPE>::Array(ShortInt dim, Integer size)
{
  if (dim <= 0 || size <= 0)
    Error("Array: wrong dimension/size", __FILE__,__LINE__);

  dim_ = dim;
  size_ = size;

  sol_ = new Vector<TYPE>[dim_];
  for (Integer i=0; i<dim; i++)
    sol_[i].Resize(size);

}

template<class TYPE>
Array<TYPE>::Array(const Array& x)
{
  dim_ = x.dim_;
  size_ = x.size_;

  sol_ = new Vector<TYPE>[dim_];
  
  for (Integer i=0; i<dim_; i++)
    sol_[i] = x.sol_[i];

}

template<class TYPE>
Array<TYPE>::Array(std::vector<TYPE> & v)
{
  dim_ = 1;
  size_ = v.size();

  
  sol_ = new Vector<TYPE>[dim_];
  sol_[0].Resize(v.size());

  for (Integer i=0; i<size_; i++)
    sol_[0][i] = v[i];

}

template<class TYPE>
Array<TYPE>::Array(Vector<TYPE> & v)
{
  dim_ = 1;
  size_ = v.size();
  
  sol_ = new Vector<TYPE>[dim_];
  sol_[0] = v;

}

template<class TYPE>
Array<TYPE>::~Array()
{
  if (sol_)
    delete[] sol_;
}

template<class TYPE>
void Array<TYPE>::clear()
{
  dim_ = 0;
  size_ = 0;
  if( sol_)
    delete[] sol_;
  
}

template<class TYPE>
void Array<TYPE>::init()
{
  for (Integer i=0; i<dim_; i++)
    sol_[i].Init();   
}


template<class TYPE>
void Array<TYPE>::redim(Integer dim)
{

  if (dim <= 0)
    Error ("Array: invalid dimension",__FILE__,__LINE__);

  if (dim != dim_)
    {
      dim_ = dim;

      if (sol_)
	delete[] sol_;
      
      sol_ = new Vector<TYPE>[dim_];

      if (size_ > 0)
	for (Integer i=0; i<dim_; i++)
	  sol_[i].Resize(size_);
    }
}

template<class TYPE>
void Array<TYPE>::resize(Integer size)
{
  if (size <= 0)
    Error ("Array: invalid size",__FILE__,__LINE__);

  
   if (size != size_)
    {
      size_ = size;
      
      if (dim_ >0)
	{
	  if (!sol_)
	    sol_ = new Vector<TYPE>[dim_];
	  
	  for (Integer i=0; i<dim_; i++)
	    sol_[i].Resize(size);
	}
    }

}

template<class TYPE>
void Array<TYPE>::reshape(ShortInt dim, Integer size)
{
  if (size <= 0 || size <= 0)
    Error ("Array: invalid size/dimension",__FILE__,__LINE__);


  if ( (size != size_) || (dim != dim_))
    {
      dim_ = dim;
      size_ = size;

      if (sol_)
	delete[] sol_;
      
      sol_ = new Vector<TYPE>[dim_];
      
      for (Integer i=0; i<dim_; i++)
	sol_[i].Resize(size);
    }

}


template<class TYPE>
void Array<TYPE>::setValuesRow(Vector<TYPE> & v, Integer pos)
{

   if (v.size() != dim_)
    Error("Array<TYPE>setValuesRow: vector has wrong dimension",__FILE__,__LINE__);

   if (pos >= size_)
      Error("Array<TYPE>setValuesRow: index out of bounds",__FILE__,__LINE__);


   for (Integer i=0; i<dim_; i++)
     sol_[i].p[pos] = v.p[i];
}

template<class TYPE>
void Array<TYPE>::setValuesColumn(Vector<TYPE> & v, Integer pos)
{

  if (v.size() != size_)
    Error("Array<TYPE>SetValuesColumn: vector has wrong dimension",__FILE__,__LINE__);
  
  if (pos >= dim_)
    Error("Array<TYPE>SetValuesColumn: wrong dimension",__FILE__,__LINE__);
  
  sol_[pos] = v;
  
}

template<class TYPE>
void Array<TYPE>::push_back(Array & arr)
{

  if (arr.dim_ != dim_)
   Error("Array<TYPE>push_back: arrays of different dimensions",__FILE__,__LINE__);

  for (Integer i=0; i<dim_; i++)
    sol_[i].add(arr.sol_[i],size_);
 
  size_ = sol_[0].size();
}

template<class TYPE>
void Array<TYPE>::push_back(Vector<TYPE> & v)
{

  if (v.size() != dim_)
   Error("Array<TYPE>push_back: arrays of different dimensions",__FILE__,__LINE__);

  for (Integer i=0; i<dim_; i++)
    sol_[i].add(v[i],size_);
 
  size_ = sol_[0].size();
}


template<class TYPE>
Array<TYPE>& Array<TYPE>::operator= (const Array & x)
{
  
  if (!sol_)
    sol_ = new Vector<TYPE>[x.dim_];

  if (sol_ && (dim_ != x.dim_) )
    {
      delete[] sol_;
      sol_ = new Vector<TYPE>[x.dim_];
    } 
  
  dim_ = x.dim_;
  size_ = x.size_;
  
  for (Integer i=0; i<dim_; i++)
    sol_[i] = x.sol_[i];
}

template<class TYPE>
Array<TYPE>& Array<TYPE>::operator= (const std::vector<TYPE> x)
{
  
  if (!sol_)
    sol_ = new Vector<TYPE>[1];

  if (sol_ && (dim_ != 1 ))
    {
      delete[] sol_;
      sol_ = new Vector<TYPE>[1];
    } 

  dim_ = 1;
  size_ = x.size(); 
  sol_[0].Resize(size_);
  
  for (Integer i=0; i<size_; i++)
    sol_[0][i] = x[i];
  
}

template<class TYPE>
Array<TYPE>& Array<TYPE>::operator= (const Vector<TYPE> x)
{

  if (!sol_)
    sol_ = new Vector<TYPE>[1];

  if (sol_ && (dim_ != 1) )
    {
      delete[] sol_;
      sol_ = new Vector<TYPE>[1];
    }
  
  dim_ = 1;
  size_ = x.size();

  sol_[0] = x;
  

}

template<class TYPE>
Array<TYPE> Array<TYPE>::operator+ () const
{
  return *this;
}

template<class TYPE>
Array<TYPE> Array<TYPE>::operator+ (const Array & x) const
{
  
  if (x.dim_ != dim_ || x.size_ != size_ || !sol_)
    Error("Array<TYPE>operator+: arrays have different size/dimension",__FILE__,__LINE__);
  
  if (dim_<= 0 || size_ <=0)
    Error("Array: not inialized",__FILE__,__LINE__);

  Array ret(dim_, size_);

  for (Integer i=0; i<dim_; i++)
    ret.sol_[i] = sol_[i] +x.sol_[i];
  
  return ret;
}

template<class TYPE>
Array<TYPE>& Array<TYPE>::operator+= (const Array & x)
{
  
   if (x.dim_ != dim_ || x.size_ != size_ || ! sol_) 
    Error("Array<TYPE>operator+=: arrays have different size/dimension",__FILE__,__LINE__);
 
   if (dim_<= 0 || size_ <=0)
     Error("Array: not inialized",__FILE__,__LINE__);

  for (Integer j=0; j<dim_; j++)
    sol_[j] += x.sol_[j];

  return *this;
}

template<class TYPE>
Array<TYPE> Array<TYPE>::operator- () const
{

  if (dim_<= 0 || size_ <=0 || ! sol_)
     Error("Array: not inialized",__FILE__,__LINE__);
  
  Array ret(dim_, size_);
  
  for (Integer i=0; i<dim_; i++)
    ret.sol_[i] = - sol_[i];

  return ret;
  
}

template<class TYPE>
Array<TYPE> Array<TYPE>::operator- (const Array & x) const
{

  if (x.dim_ != dim_ || x.size_ != size_ || ! sol_)
    Error("Array<TYPE>operator-: arrays have different size/dimension",__FILE__,__LINE__);

  if (dim_<= 0 || size_ <=0)
     Error("Array: not inialized",__FILE__,__LINE__); 
  
  Array ret(dim_, size_);
  
  for (Integer i=0; i<dim_; i++)
    ret.sol_[i] = sol_[i] - x.sol_[i];

  return ret;
  
}

template<class TYPE>
Array<TYPE>& Array<TYPE>::operator-= (const Array & x)
{
  if (x.dim_ != dim_ || x.size_ != size_ || !sol_)
    Error("Array<TYPE>operator-=: arrays have different size/dimension",__FILE__,__LINE__);
 
  if (dim_<= 0 || size_ <=0)
     Error("Array: not inialized",__FILE__,__LINE__);
  
  for (Integer i=0; i<dim_; i++)
    sol_[i] -= x.sol_[i];

  return *this; 
}

template<class TYPE>
TYPE Array<TYPE>::operator* (const Array & x) const
{
   if (x.dim_ != dim_ || x.size_ != size_)
    Error("Array<TYPE>operator*: arrays have different size/dimension",__FILE__,__LINE__);
  
  if ( x.dim_ != 1)
    Error("Array<TYPE>operator*: not defined for arrays with dimension > 1",__FILE__,__LINE__);
  
  if (dim_<= 0 || size_ <=0 || ! sol_)
     Error("Array: not inialized",__FILE__,__LINE__);
  
 TYPE ret = 0;
  
  for (Integer j=0; j<size_; j++) 
    ret += sol_[0][j] * x.sol_[0][j];
  
  return ret;   
}

template<class TYPE>
Array<TYPE> Array<TYPE>::operator* (TYPE x) const
{
  
  if (dim_<= 0 || size_ <=0 || !sol_)
     Error("Array: not inialized",__FILE__,__LINE__);
  
  Array ret(dim_, size_);
  
  for (Integer i=0; i<dim_; i++)
    ret.sol_[i] = sol_[i] * x;
  
  return ret;      
}

template<class TYPE>  
Array<TYPE> Array<TYPE>::operator/ (TYPE x) const
{
  if (dim_<= 0 || size_ <=0 || !sol_)
     Error("Array: not inialized",__FILE__,__LINE__);
  
  Array ret(dim_, size_);
  
  for (Integer i=0; i<dim_; i++)
    ret.sol_[i] = sol_[i] / x;
  
  return ret;    
}

template<class TYPE>
Array<TYPE>& Array<TYPE>::operator/= (TYPE x)
{
  if (dim_<= 0 || size_ <=0 || !sol_)
    Error("Array: not inialized",__FILE__,__LINE__);
  
 for (Integer i=0; i<dim_; i++)
     sol_[i] /= x;

  return *this;    
}

template<class TYPE>  
Double Array<TYPE>::normL2 (Integer pos)
{
  if (!sol_)
    Error("Array not initalized",__FILE__,__LINE__);

  if (pos >= size_)
     Error("Array<TYPE>normL2: index out of bounds",__FILE__,__LINE__);

  Double ret = 0;

  for (Integer i=0; i<dim_; i++)
    ret += sol_[i][pos] * sol_[i][pos];

  return sqrt(ret);

}

template<class TYPE>
Double Array<TYPE>::normL2()
{
  if (!sol_)
    Error("Array not initalized",__FILE__,__LINE__);
  
  Double ret = 0;
  Double temp;

  for (Integer j=0; j<size_; j++)
    {
      temp = 0.0;
      for (Integer i=0; i<dim_; i++)
	temp += sol_[i][j] * sol_[i][j];

      ret += sqrt(temp);
    }

  return ret;
}

template<class TYPE> 
void Array<TYPE>::toMatrix(Matrix<TYPE> & m)
{

  Error("Not implemented yet",__FILE__,__LINE__);
    

}

template<class TYPE>
void Array<TYPE>::toStdVector(std::vector<TYPE> & v, Integer dim)
{

  if (dim >= dim_ || dim <= 0)
    Error("Array: wrong dimension",__FILE__,__LINE__);
 
 v.resize(size_);

  for (Integer i=0; i<size_; i++)
    v[i] = sol_[dim][i];  
}

template<class TYPE> 
void Array<TYPE>::toVector(Vector<TYPE> & v, Integer dim)
{

  if (dim >= dim_ || dim <= 0)
    Error("Array: wrong dimension",__FILE__,__LINE__);

   for (Integer i=0; i<size_; i++)
    v[i] = sol_[dim][i];

}

template<class TYPE>
std::ostream& operator<< (std::ostream &out, const Array<TYPE> &a)
{
    for (Integer i=0; i<a.size(); i++)
    {
      out << i << ": ";
      
      for (Integer j=0; j<a.dim(); j++)
	out << a[j][i] << " ";
      
      out << std::endl;
    }
  
  return out;
}

template std::ostream& operator<<<Integer> (std::ostream &out, const Array<Integer> &);
template std::ostream& operator<<<Double> (std::ostream &out, const Array<Double> &);

} // end of namespace  



