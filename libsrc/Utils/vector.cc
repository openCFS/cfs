#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

#include "general_head.hh"
#include "tools.hh"
#include "vector.hh"
#include "matrix_head.hh"

namespace CoupledField
{ 

template<class TYPE>
Vector<TYPE>::Vector ()
{
	n = 0;
	p = NULL;
}

template<class TYPE> Vector<TYPE>::Vector(Integer i)
{
   if (i <= 0)
Error("Vector: invalid dimension in constructor",__FILE__,__LINE__);
	n = i;
	p = new TYPE [n];
	if (!p)
		Error("Vector: out of memory",__FILE__,__LINE__);

	for (Integer j = 0; j < n; j++)
		p [j] = 0;
}

template<class TYPE> void Vector<TYPE>::Allocate(const Integer i)
{
   if (p)
   Error("use .Resize to defined Vector for changing size",__FILE__,__LINE__);

   n=i;
   p=new TYPE[n];

   if (!p) Error("Out of memory", __FILE__,__LINE__);

}

template<class TYPE> void Vector<TYPE>::Resize(const Integer i)
{
        if (i <= 0) Error("invalid dimension for Resize", __FILE__, __LINE__);

        if (p) delete [] p;

        n=i;
        p = new TYPE [n];
        if (!p) Error("Out of memory in vector.Resize()");
 
        for (Integer j = 0; j < n; j++)
                p [j] = 0;
}

template<class TYPE>
Vector<TYPE>::Vector (const Integer i, const TYPE *const x)
{	if (i <= 0)
  Error("Vector: invalid dimension in constructor",__FILE__,__LINE__);
	n = i;
	p = new TYPE [n];
	if (!p)
		cerr << "Vector: out of memory";

	for (Integer j = 0; j < n; j++)
		p [j] = x [j];
}

template
<class TYPE>
Vector<TYPE>::Vector (const Vector<TYPE> &x)
{  
    if (!x.n)
        Error("Vector: undefined Vector in copy-constructor",__FILE__,__LINE__);

	n = x.n;
	p = new TYPE [n];
	if (!p)
              Error("Vector: out of memory",__FILE__,__LINE__);

	for (Integer i = 0; i < n; i++)
		p [i] = x.p [i];
}

template<class TYPE>
void Vector<TYPE>::Init(const Integer l)
{
 if (!n) Error("Don't use Init() to undefined vector", __FILE__, __LINE__);
 
 Integer i;
 for (i=0; i<n; i++) p[i]=l;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator=(const Vector<TYPE> &x)
{
//	if (x.n==0)
//		{ if (p) delete [] p;
//                  p=NULL;
//                  n=0; 
//                  return * this;}
//if (!x.n) 
//    Error("Undefined vector for copy-constructor",__FILE__,__LINE__);

	if (this == &x)
		return *this;

	if (n != x.n)
	{	if (p)
			delete [] p;
		n = x.n;
		p = new TYPE [n];
		if (!p)
		cerr << "Vector: out of memory";
	}

	for (Integer i = 0; i < n; i++)
		p [i] = x.p [i];

	return *this;
}

template<class TYPE>
TYPE &Vector<TYPE>::operator[] (const Integer i) const
{	if (!n)
Error("Vector: undefined Vector in operator[]",__FILE__,__LINE__);

	if (i < 0 || i >= n)
Error("Vector: invalid index in operator[]",__FILE__,__LINE__);

	return  p[i];
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator+ () const
{	if (!n)
Error("Vector: undefined Vector in operator +()",__FILE__,__LINE__);
	return *this;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator+(const Vector<TYPE> &x) const
{	if (!n || !x.n)
Error("Vector: undefined Vector in operator +(vector)",__FILE__,__LINE__);

	if (n != x.n)
Error("Vector: incompatible dimension for operator +(vector)",__FILE__,__LINE__);

	Vector	z (n);

	for (Integer i = 0; i < n; i++)
		z.p[i] = p[i] + x.p [i];

	return z;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator+=(const Vector<TYPE> &x)
{	if (!n || !x.n)
Error("Vector: undefined Vector in operator+=(vector)",__FILE__,__LINE__);

	if (n != x.n)
Error("Vector: incompatible dimension for operation +=",__FILE__,__LINE__); 

	for (Integer i = 0; i < n; i++)
		p [i] += x.p [i];

	return *this; 
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator- () const
{	if (!n)
Error("Vector: undefined Vector in oprator -()",__FILE__,__LINE__); 

	Vector	z (n);

	for (Integer i = 0; i < n; i++)
		z.p [i] = -p [i];

	return z;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator-(const Vector<TYPE> &x) const
{	if (!n || !x.n)
 Error("Vector: undefined Vector in operator -(vector)",__FILE__,__LINE__); 

	if (n != x.n)
 Error("Vector: incompatible dimension for operation -(vector)",__FILE__,__LINE__); 

	Vector	z(n);

	for (Integer i = 0; i < n; i++)
		z.p [i] = p [i] - x.p [i];

	return z;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator-=(const Vector<TYPE> &x)
{	if (!n || !x.n)
Error("Vector: undefined Vector in operator -=(vector)",__FILE__,__LINE__); 

	if (n != x.n)
Error("Vector: incompatible dimension in operation -=(vector)",__FILE__,__LINE__); 

	for (Integer i = 0; i < n; i++)
		p [i] -= x.p [i];

	return *this;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator* (const TYPE &x) const
{	if (!n)
Error("Vector: undefined Vector in operator *(number)",__FILE__,__LINE__); 

	Vector	z (n);

	for (Integer i = 0; i < n; i++)
		z.p[i] = p [i] * x;

	return z;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator/ (const TYPE &x) const
{       if (!n)
Error("Vector: undefined Vector in operator/",__FILE__,__LINE__); 
 
        Vector  z (n);
 
        for (Integer i = 0; i < n; i++)
                z.p[i] = p [i] / x;
 
        return z;
}

template<class TYPE>
TYPE Vector<TYPE>::operator* (const Vector<TYPE> &x) const
{	if (!n || !x.n)
Error("Vector: undefined Vector in operator *(vector)",__FILE__,__LINE__);

	if (n != x.n)
Error("Vector: incompatible dimension in operator *(vector)",__FILE__,__LINE__);  

	TYPE	z;

	z = p [0] * x.p [0];
	for (Integer i = 1; i < n; i++)
		z += p [i] * x.p [i];

	return z;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::operator*(const Matrix<TYPE> &x) const
{	if (!n)
		Error( "undefined Vector in operator ", __FILE__,__LINE__);

	if (!x.row || !x.col)
                Error( "undefined Vector", __FILE__,__LINE__);

	if (n != x.row)
                Error( "incompatible dimension", __FILE__,__LINE__);

	TYPE	a;
	Vector	z(x.col);

	for (Integer i = 0; i < x.col; i++)
	{	a = p [0] * x.p [0] [i];
		for (Integer j = 1; j < n; j++)
			a += p [j] * x.p [j] [i];
		z.p [i] = a;
	}

	return z;
}

template<class TYPE>
Vector<TYPE> &Vector<TYPE>::operator*= (const TYPE &x)
{	if (!n)
Error("Vector: undefined Vector in operator*=",__FILE__,__LINE__); 

	TYPE	y = x;

	for (Integer i = 0; i < n; i++)
		p [i] *= y;

	return *this;
}

template<class TYPE>

Vector<TYPE> &Vector<TYPE>::operator*=(const Matrix<TYPE> &x)
{	*this = *this * x;

	return *this;
}

template<class TYPE>
Integer Vector<TYPE>::operator== (const Vector<TYPE> &x) const
{	if (!n || !x.n)
  Error("Vector: undefined Vector in operator==",__FILE__,__LINE__);

	for (Integer i = 0; i < n; i++)
		if (p [i] != x.p [i])
			return 0;

	return 1;
}

template<class TYPE>
Integer Vector<TYPE>::operator!= (const Vector<TYPE> &x) const
{	if (!n || !x.n)
  Error("Vector: undefined Vector in operator !=", __FILE__, __LINE__);

	for (Integer i = 0; i < n; i++)
		if (p [i] != x.p [i])
			return 1;

	return 0;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::part (const Integer i1,const Integer i2) const
{	if (!n)
		cerr << "Vector: undefined Vector";

	if (i1 < 0 || i1 > i2 || i2 >= n)
		cerr << "Vector: invalid index";

	Vector	z (i2 - i1 + 1);

	for (Integer i = 0; i < z.n; i++)
		z [i] = p [i1 + i];

	return z;
}

template<class TYPE>
Vector<TYPE> Vector<TYPE>::unit (const Integer n,const Integer i)
{	if (n <= 0)
		cerr << "Vector: invalid dimension";

	if (i < 0 || i >= n)
		cerr << "Vector: invalid index";

	Vector<TYPE>	z (n);

	z.p [i] = 1;

	return z;
}



//template<class TYPE>
//Vector<TYPE> operator* (const TYPE &x,const Vector<TYPE> &y)
//{	if (!y.size())
//		cerr << "Vector: undefined Vector";
//
//	Vector<TYPE>	z (y.size());
//
//	for (Integer i = 0; i < y.size); i++)
//		z.p [i] = x * y.p [i];
//
//	return z;
//}


template<class TYPE>
Double Vector<TYPE> :: norm_2 ()
{	if (!(*this).size())
  Error("Vector: undefined Vector in function norm_2()", __FILE__, __LINE__);

	TYPE z = (*this)*(*this);
	return sqrt (z);
}

template<class TYPE>
void  Vector<TYPE>:: add (const TYPE & y, Integer pos)
{
   Integer i;

   TYPE * help=new TYPE[n+1];

   for (i=0; i < pos; i++) help[i]=p[i];
   help[pos]=y;
   for (i=pos+1; i<n+1; i++) help[i]=p[i-1];

   delete [] p;
   p=help;
   n++;
}

template<class TYPE>
void  Vector<TYPE>:: add (const Vector<TYPE> & y, Integer pos)
{
   Integer i;
   
   Integer l=y.size();
   TYPE * help=new TYPE[n+l];
   for (i=0; i < pos; i++) help[i]=p[i];
   for (i=pos; i < pos+l; i++) help[i]=y[i-pos];
   for (i=pos+l; i < n+l; i++) help[i]=p[i-l];
  
   delete [] p; 
   p=help;
   n+=l;
}

template<class TYPE>
void  Vector<TYPE>:: cut (const Integer pos)
{
   if (pos<0 || pos >=n) Error("Invalid index for cut");
   if (!n) Error("Using .cut() to undefined vector");
 
   Integer i;
 
   TYPE * help=new TYPE[n-1];
   for (i=0; i < pos; i++) help[i]=p[i];
   for (i=pos+1; i < n; i++) help[i-1]=p[i];
 
   delete [] p; 
   p=help;
   n=n-1;
}

template<class TYPE>
void  Vector<TYPE>:: cut (const Integer pos1, const Integer pos2)
{
   Integer i;
 
   Integer l=pos2-pos1+1;
   TYPE * help=new TYPE[n-l];
   for (i=0; i < pos1; i++) help[i]=p[i];
   for (i=n-1; i > pos2 ; i--) help[i-pos2+pos1-1]=p[i];

   delete [] p; 
   p=help;
   
   n=n-l;
}

template<class TYPE>
Integer  Vector<TYPE>:: memory() const
{
 return n*sizeof(TYPE);
}

template<class T>
void swap(Vector<T> & a, Vector<T> & b) 
{
 Integer tmp_size;

 tmp_size=a.size;
 a.size=b.size;
 b.size=a.size;

 T * tmp_p;

 tmp_p=a.p;
 a.p=b.p;
 b.p=a.p;
}

template<class T> void swap(T & a, T & b)
{
   T tmp=a;
   a=b;
   b=tmp;
}

template void swap<Integer>(Integer & ,Integer &);
template void swap<Double>(Double & ,Double &);

template<class S> void sort(S* v, Integer n)
{
  for (Integer gap=n/2; gap >0; gap/=2)
    for (Integer i=gap; i<n; i++)
       for (Integer j=i-gap; j>=0; j-=gap)
          if (v[j+gap]<v[j]) swap<S>(v[j],v[j+gap]);

}
template void sort<Integer>(Integer * ,Integer);
template void sort<Double>(Double * ,Integer);

// Overloading << for Vector

template<class S>
std::ostream & operator << ( std::ostream & out, const Vector<S> & vc)
{
if (vc.size()==0) out << "vector is undefined" ;
for (Integer i=0; i < vc.size(); i++)
 out << vc[i] << " " << std::endl;
 return out;
}

template std::ostream & operator<<<Integer> (std::ostream & , const Vector<Integer> &);
template std::ostream & operator<<<Double> (std::ostream & , const Vector<Double> &);

//Vector<Double> Vector<Double>::null=Vector<Double>();
//Vector<Integer> Vector<Integer>::null=Vector<Integer>();

}
