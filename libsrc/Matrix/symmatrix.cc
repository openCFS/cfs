#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>

#include "symmatrix.hh"

namespace CoupledField
{      

template<class TYPE>
SymMatrix<TYPE>::SymMatrix ()
{	size = 0;
	p = NULL;
}

template<class TYPE>
SymMatrix<TYPE>::SymMatrix (const Integer i)
{
    	if (i <= 0) Error("invalid dimension",__FILE__,__LINE__);

	size = i;
	p = new TYPE* [size];

       if (!p) Error("out of memory",__FILE__,__LINE__); 

       p[0]=new TYPE[size*(size+1)/2];

      if (!p[0]) Error("out of memory",__FILE__,__LINE__); 

      Integer k;
      for (k=1; k < size; k++) p[k]=p[k-1]+k;
//      for (k=0; k < size*(size+1)/2; k++)  p[0][k]=0;
}

template<class TYPE>
SymMatrix<TYPE>::SymMatrix (const SymMatrix<TYPE> &x)
{	if (!x.size)  Error("undefined SymMatrix");

	size = x.size;
        
        if (!p) { delete p[0];
                  delete [] p; }

	p = new TYPE * [size];
	if (!p) Error("out of memory", __FILE__, __LINE__);

        p[0]=new TYPE[size*(size+1)/2];
        if (!p[0]) Error("out of memory", __FILE__, __LINE__);

        Integer k;

        for (k=0; k < size*(size+1)/2; k++)  p[0][k]=x.p[0][k];
        for (k=1; k < size; k++) p[k]=p[k-1]+k;        
}

template<class TYPE>
void SymMatrix<TYPE> :: Resize(const Integer isize, const Integer isizefikc)
{
   
   Integer k;

   if (p) {delete p[0];
           delete [] p;}

   size=isize;

   p = new TYPE* [size];
   p[0]=new TYPE[size*(size+1)/2];
   
   if (!p && !p[0]) Error("Out of memory in Resize() function");

  for (k=1; k < size; k++) p[k]=p[k-1]+k;
}

template<class TYPE>
TYPE * SymMatrix<TYPE>::operator[] (const Integer i) const
{ 
        if (!size) Error("undefined SymMatrix",__FILE__,__LINE__);
        if (i < 0 || i >= size) Error("invalid index",__FILE__,__LINE__);
 
        return p[i];
}

template<class TYPE>
SymMatrix<TYPE> &SymMatrix<TYPE>::operator=(const SymMatrix<TYPE> &x)
{
       if (!x.size) Error("undefined SymMatrix",__FILE__,__LINE__);
 
        if (this == &x)  return *this;
 
        Integer k;

        if (size!=x.size )
        { 
           if (p)   {delete p[0];
                     delete [] p;}

           size=x.size; 

           p = new TYPE* [size];
           p[0]=new TYPE[size*(size+1)/2];

           for (k=1; k < size; k++) p[k]=p[k-1]+k;
           for ( k = 0; k < size*(size+1)/2; k++) p[0][k]=0;
 
           if (!p) Error("out of memory",__FILE__,__LINE__);
           if (!p[0]) Error("out of memory",__FILE__,__LINE__);
        }
         
        for (k = 0; k < size*(size+1)/2; k++)
           p[0][k]=x.p[0][k];

        return *this;
}

template<class TYPE>
SymMatrix<TYPE> SymMatrix<TYPE>::operator+ () const
{
        if (!size) Error("undefined SymMatrix",__FILE__,__LINE__);
        return *this;
}

template<class TYPE>
SymMatrix<TYPE> SymMatrix<TYPE>::operator+(const SymMatrix<TYPE> &x) const
{
        if (!size)
            Error("undefined SymMatrix",__FILE__,__LINE__);
 
        if (size != x.size)
             Error("incompatible dimension",__FILE__,__LINE__);

        SymMatrix<TYPE> z(size);

         Integer k;
         for ( k = 0; k < size*(size+1)/2; k++)
                z [0][k] = x.p [0][k]+p[0][k];

        return z;
}

template<class TYPE>
SymMatrix<TYPE> &SymMatrix<TYPE>::operator+=(const SymMatrix<TYPE> &x)
{
        if (!size || !x.size )
            Error("Don't use += to undefuned matrix",__FILE__,__LINE__);
 
         if (size != x.size )
             Error("incompatible dimension for +",__FILE__,__LINE__); 
 
        Integer k;
        for ( k = 0; k < size*(size+1)/2; k++)
                p [0][k] += x.p [0][k];
 
        return *this;
}

template<class TYPE>
SymMatrix<TYPE> SymMatrix<TYPE>::operator-() const
{
        if (!size )
            Error("undefined SymMatrix",__FILE__,__LINE__);
 
        SymMatrix<TYPE> z(size);
 
         Integer k;
         for ( k = 0; k < size*(size+1)/2; k++)
              if (p[0][k]!=0)  z [0][k] = -p [0][k];
 
        return z;
}

template<class TYPE>
SymMatrix<TYPE> SymMatrix<TYPE>::operator-(const SymMatrix<TYPE> &x) const
{
        if (!size || !x.size)
            Error("undefined SymMatrix",__FILE__,__LINE__);
 
        if (size != x.size )
             Error("incompatible dimension",__FILE__,__LINE__);
 
        SymMatrix<TYPE> z(size);
 
         Integer k;
         for ( k = 0; k < size*(size+1)/2; k++)
                z [0][k] = -x.p [0][k]+p[0][k];
 
        return z;
}

template<class TYPE>
SymMatrix<TYPE> &SymMatrix<TYPE>::operator-=(const SymMatrix<TYPE> &x)
{
        if (!size || !x.size)
            Error("Don't use += to undefuned matrix",__FILE__,__LINE__);
 
         if (size != x.size)
             Error("incompatible dimension for +",__FILE__,__LINE__);
 
        Integer k;
        for ( k = 0; k < size*(size+1)/2; k++)
                p [0][k] -= x.p [0][k];
 
        return *this;
}

template<class TYPE>
Integer SymMatrix<TYPE>::operator== (const SymMatrix<TYPE> &x) const
{
         if (!size || !x.size)
              Error("undefuned matrix",__FILE__,__LINE__);

         Integer k;

         for (k = 0; k < size*(size+1)/2; k++)
              if (p [0][k] != x.p[0][k]) return 0;
 
        return 1;
}

template<class TYPE>
Integer SymMatrix<TYPE>::operator!= (const SymMatrix<TYPE> &x) const
{
         if (!size || !x.size)
              Error("undefuned matrix",__FILE__,__LINE__); 

         Integer k;
         for (k = 0; k < size*(size+1)/2; k++)
              if (p [0][k] != x.p[0][k]) return 1;
 
        return 0;
}

template<class TYPE>
SymMatrix<TYPE> SymMatrix<TYPE>::operator* (const TYPE &x) const 
{       if (!size) Error("undefined SymMatrix",__FILE__,__LINE__);
 
        Integer k;

        SymMatrix<TYPE> z(size);

        for ( k = 0; k < size*(size+1)/2; k++)
                z[0][k] = p[0][k]*x;
 
        return z;
}

template<class TYPE>
Vector<TYPE> SymMatrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
        if (!size) Error("undefined SymMatrix",__FILE__,__LINE__);
 
        if (!x.n) Error("undefined Vector",__FILE__,__LINE__);
 
        if (size != x.n) Error("incompatible dimension",__FILE__,__LINE__);
 
        Vector<TYPE> z(size);
 
        Integer k,kk;
        for ( k = 0; k < size; k++)
        {
           for ( kk = 0; kk < k+1; kk++)
                z.p[k] += p[k][kk]*x[kk];
           for ( kk=k+1; kk < size; kk++)
                z.p[k] += p[kk][k]*x[kk];
        } 
        return z;
}

template<class TYPE>
Matrix<TYPE> SymMatrix<TYPE>::operator*(const Matrix<TYPE> &x) const
{
         if (!size || !x.col || !x.row)
         Error("undefined SymMatrix or Matrix",__FILE__,__LINE__);
 
        if (size != x.row)
        Error("incompatible dimension",__FILE__,__LINE__);
 
        TYPE    a;
        Matrix<TYPE>  z (size, x.col);

        Integer i,j,k; 
        for (i = 0; i < size; i++)
                for (j = 0; j < x.col; j++)
                {       a = p [i] [0] * x.p [0] [j];
                        for ( k = 1; k < i+1; k++)
                                a += p [i] [k] * x.p [k] [j];
                        for ( k=i+1; k < size; k++)
                                a+=  p [k] [i] * x.p [k] [j];
                        z.p [i] [j] = a;
                }
 
        return z;
}

template<class TYPE>
Matrix<TYPE> SymMatrix<TYPE>::operator*(const SymMatrix<TYPE> &x) const
{
         if (!size || !x.size)
         Error("undefined SymMatrix ",__FILE__,__LINE__);
 
        if (size != x.size)
        Error("incompatible dimension",__FILE__,__LINE__);
 
        TYPE    a;
        Matrix<TYPE>  z (size, size);
 
        Integer i,j,k;
        for (i = 0; i < size; i++)
                for (j = 0; j < size; j++)
                {
                        a = p [i] [0] * x.p [j] [0];
                        for ( k = 1; k < i+1; k++)
                                if (k<j)
                                a += p [i] [k] * x.p [j] [k];
                                else 
                                  a += p [i] [k] * x.p [k] [j];
                        for ( k=i+1; k < size; k++)
                                if (k<j)
                                a+=  p [k] [i] * x.p [j] [k];
                                else
                                a+=  p [k] [i] * x.p [k] [j];
                        z.p [i] [j] = a;
                }
 
        return z;
}

template<class TYPE>
SymMatrix<TYPE> &SymMatrix<TYPE>::operator*= (const TYPE &x)
{
         if (!size ) Error("undefined matrix",__FILE__,__LINE__);

        TYPE y=x;
 
        Integer i;
        for (i = 0; i < size*(size+1)/2; i++)
                p [0][i] *= y;
 
        return *this;
}

template<class TYPE>
void SymMatrix<TYPE>::cut(const Integer ai, const Integer aj)
{
if (ai < 0 || ai > size-1 || aj < 0 || aj > size-1)
         Error("invalid index in function cut",__FILE__,__LINE__);
if (!size) Error("SymMatrix is undefined in function cut");
 
SymMatrix<TYPE> help(size-1);
Integer n=size*(size-1)/2;
 
Integer i, j, k=0;
for (i=0; i<size ; i++)
if (i!=ai)  for (j=0; j<size; j++)
if (j!=aj)
   {
     help.p[0][k]=p[i][j];
     k++;
   }
 
delete  p[0];
delete [] p;
size--; 
 
p = new TYPE* [size];
p[0]=new TYPE[n];
for (k=1; k < size; k++) p[k]=p[k-1]+k;
 
for (i=0; i<n; i++) p[0][i]=help.p[0][i];
}

template<class TYPE>
SymMatrix<TYPE> SymMatrix<TYPE>::part (const Integer i) const
{
         if (!size ) Error("undefined SymMatrix",__FILE__,__LINE__); 
 
        if (i < 0 || i >= size)
                      Error("invalid index",__FILE__,__LINE__);
 
        SymMatrix  z (i + 1);
 
        for (Integer k = 0; k < (z.size+1)*z.size/2 ; k++)
                        z [0] [k] = p [0] [k];
 
        return z;
}

template< class TYPE>
void SymMatrix<TYPE>::precond(Vector<TYPE> & e, const Vector<TYPE> r, const Integer type)
{
#ifdef TRACE
  (*trace) << "entering SymMatrix::precond" << std::endl;
#endif
 
 Integer i,j;
 Integer n=r.size();
 Double omega;
 
 switch(type)
{ case non: e=r; break;
  case Jacobi:
       for (i=0; i<n; i++) e[i]=r[i]/p[i][i];
       break;
  case SSOR:
       omega=1.2; /// omega=1.2
       TYPE sum;
       for (i=0; i<n; i++)
       {
         sum=0;
         for (j=0; j<i; j++) sum+=p[i][j]*e[j];
         e[i]=(r[i]-omega*sum)/p[i][i];
       }
 
       for (i=n-1; i>=0; i--)
       {
         sum=0;
         for (j=i+1; j<n; j++) sum+=p[j][i]*e[j];
         e[i]-=omega*sum/p[i][i];
       }
 
       e*=omega*(2-omega);
       break;

default:
       Error("Wrong type of Matrix preconditioning",__FILE__,__LINE__);
 
} //end of switch

}

template<>
void SymMatrix<Double>::CholerskyDecomposition()
{
#ifdef TRACE
  (*trace) << " entering CholerskyDecomposition\n";
#endif

  // CholerskyDecomposition A=BB^T.
  // We store B matrix as lower triangular matrix in memory of A. It means that A is destroyed (:.

  p[0][0]=std::sqrt(p[0][0]);
  
  Integer i,j,k;
  for (i=1; i<size; i++) 
    {
      Double sum=0;

      for (j=0; j<i; j++)
	{
	  sum=0;
	  for (k=0; k<j; k++)
	    sum+=p[i][k]*p[j][k];

	  p[i][j]=(p[i][j]-sum)/p[j][j];
	}

      sum=0;
       for (j=0; j<i; j++)
	sum+=sqr(p[i][j]);
	
      p[i][i]=std::sqrt(p[i][i]-sum);
    } 
}

template<>
void SymMatrix<Integer>::CholerskyDecomposition()
{
  Error("For Integer matrix is not implemented");
}

template<class S>
std::ostream & operator << (std::ostream & out, const SymMatrix<S> & mat)
{
Integer i,j;

out.setf(std::ios::scientific);

for (i=0; i < mat.getSize(); i++)
{
  for (j=0; j < i+1; j++) out << mat[i][j] << " ";

// for printing matrix in full size, uncomment line below
//  for (jj=i+1; jj < mat.getSize(); jj++) out << mat[jj][i] <<" ";

 out << std::endl;
}

//out.setf(0, std::ios::floatfield);

 return out;
}

//template std::ostream & operator<<<Integer> (std::ostream & , const SymMatrix<Integer> &);
template std::ostream & operator<<<Double> (std::ostream & , const SymMatrix<Double> &);

template<class TYPE>
TYPE Spur (const SymMatrix<TYPE> &x)
{
if (!x.getSize()) Error("undefined SymMatrix",__FILE__,__LINE__);
 
        TYPE    a;
 
        a = x.get() [0] [0];
 
        Integer i;
        for ( i = 1; i < x.getSize(); i++)
                a += x.get() [i] [i];
 
        return a;
}

template Integer Spur<Integer>(const SymMatrix<Integer> &);
template Double Spur<Double>(const SymMatrix<Double> &);

} // end of namespace
