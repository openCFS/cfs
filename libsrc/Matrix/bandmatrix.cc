#include <iostream>
#include <fstream>
#include <time.h>
#include <string>

#include <general_head.hh>
#include <utils_head.hh>

namespace CoupledField
{

template<class TYPE>
BandMatrix<TYPE> :: BandMatrix()
{
 row=0;
 li=0;
 ri=0;

 p=NULL;
}

template<class TYPE>
BandMatrix<TYPE>::BandMatrix (const Integer arow, const Integer alcln, const Integer arcln)
{
#ifdef TRACE
    (*trace) << "entering BandMatrix::BandMatrix" << endl;
#endif
        row = arow;
        li=alcln;
        ri=arcln;

        numentry=(row-2)*(ri+li+1)+ri+li+2;

        p = new TYPE* [row];
        p[0]=new TYPE[numentry];

      if (!p || !p[0]) Error("out of memory",__FILE__,__LINE__);
 
      Integer k;

      for (k=1; k < row; k++) p[k]=p[k-1]+ri+li+1;
}

template<class TYPE>
BandMatrix<TYPE>::BandMatrix (const BandMatrix<TYPE> &x)
{       if (!x.row || !x.li ||!x.ri)  Error("undefined BandMatrix");
 
        row = x.row;
        li = x.li;
        ri = x.ri;
        numentry=x.numentry;
 
        if (!p) { delete p[0];
                  delete [] p; }
 
        p = new TYPE * [row];
        p[0]=new TYPE[numentry];

        if (!p || !p[0]) Error("out of memory", __FILE__, __LINE__);
 
        Integer k;
 
        for (k=0; k < numentry; k++)  p[0][k]=x.p[0][k];
        for (k=1; k < row; k++) p[k]=p[k-1]+li+ri+1;
}

template<class TYPE>
void BandMatrix<TYPE> :: Resize(const Integer irow, const Integer ili, const Integer iri)
{
 
   Integer k;
 
   if (p) {delete p[0];
           delete [] p;}
 
   row=irow;
   li=ili;
   ri=iri;
   numentry=(row-2)*(ri+li+1)+ri+li+2;
 
   p = new TYPE* [row];
   p[0]=new TYPE[numentry];

   if (!p || !p[0]) Error("out of memory",__FILE__,__LINE__);
 
  for (k=1; k < row; k++) p[k]=p[k-1]+ri+li+1;
}

template<class TYPE>
TYPE * BandMatrix<TYPE>::operator[] (const Integer i) const
{
        if (!p) Error("undefined BandMatrix",__FILE__,__LINE__);
        if (i < 0 || i >= row) Error("invalid index",__FILE__,__LINE__);
 
        return p[i];
}

template<class TYPE>
BandMatrix<TYPE> & BandMatrix<TYPE>::operator=(const BandMatrix<TYPE> &x)
{
       if (!x.row) Error("undefined BandMatrix",__FILE__,__LINE__);
 
        if (this == &x)  return *this;
 
        Integer k;
 
        if (row!=x.row )
        {
           if (p)   { 
                      delete [] p[0];
                      delete [] p;}
 
           row=x.row;
           li=x.li;
           ri=x.ri;
           numentry=x.numentry;
 
           p = new TYPE* [row];
           p[0] = new TYPE[numentry];
 
           if (!p || !p[0]) Error("out of memory",__FILE__,__LINE__);

           for (k=1; k < row; k++) p[k]=p[k-1]+ri+li+1;
        }

        if (li!=x.li || ri!=x.ri)
        {
          delete [] p[0];

          li=x.li;
          ri=x.ri;
          numentry=x.numentry;

          p[0]=new TYPE[numentry];
           
           for (k=1; k < row; k++) p[k]=p[k-1]+ri+li+1;
        }
 
        for (k = 0; k < numentry; k++)  p[0][k]=x.p[0][k]; 
 
        return *this;
}

template<class TYPE>
BandMatrix<TYPE> & BandMatrix<TYPE>::operator+ () 
{
        if (!row) Error("undefined BandMatrix",__FILE__,__LINE__);
        return *this;
}

template<class TYPE>
BandMatrix<TYPE> & BandMatrix<TYPE>::operator+=(const BandMatrix<TYPE> &x)
{
  if (row != x.row || ri!=x.ri || li!=x.li) Error("Bad matrix sizes");
  for (Integer i=0; i<numentry; i++)  p[0][i]+=x[0][i];
  return *this;
}

template<class TYPE>
BandMatrix<TYPE> BandMatrix<TYPE>::operator+(const BandMatrix<TYPE> &x) const 
{
        if (!row)
            Error("undefined BandMatrix",__FILE__,__LINE__);
 
        if (row != x.row || li!=x.li || ri!=x.ri)
             Error("incompatible dimension",__FILE__,__LINE__);
        
        BandMatrix<TYPE> z(row, ri, li);
        Integer k;
        for (k=0; k< numentry; k++)  z[0][k]=x.p[0][k]+p[0][k]; 
        return z;
}

template<class TYPE>
BandMatrix<TYPE> BandMatrix<TYPE>::operator-() const
{
        if (!row )
            Error("undefined BandMatrix",__FILE__,__LINE__);
 
        BandMatrix<TYPE> z(row,li,ri);
 
         Integer k;
         for ( k = 0; k < numentry; k++)
              if (p[0][k]!=0)  z[0][k] = -p [0][k];
 
        return z;
}

template<class TYPE>
BandMatrix<TYPE> BandMatrix<TYPE>::operator-(const BandMatrix<TYPE> &x) const
{
        if (!row || !x.row)
            Error("undefined BandMatrix",__FILE__,__LINE__);
 
        if (row != x.row || ri !=x.ri || li!=x.li )
             Error("incompatible dimension",__FILE__,__LINE__);
 
        BandMatrix<TYPE> z(row,li,ri);
 
         Integer k;
         for ( k = 0; k < numentry; k++)
                z [0][k] = -x.p [0][k]+p[0][k];
 
        return z;
}

template<class TYPE>
BandMatrix<TYPE> &BandMatrix<TYPE>::operator-=(const BandMatrix<TYPE> &x)
{
        if (!row || !x.row)
            Error("Don't use -= to undefuned matrix",__FILE__,__LINE__);
 
         if (row != x.row || ri !=x.ri || li!=x.li)
             Error("incompatible dimension for -",__FILE__,__LINE__);
 
        Integer k;
        for ( k = 0; k < numentry; k++)
                p [0][k] -= x.p [0][k];
 
        return *this;
}

template<class TYPE>
Integer BandMatrix<TYPE>::operator== (const BandMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
 
         for (k = 0; k < numentry; k++)
              if (p [0][k] != x.p[0][k]) return 0;
 
        return 1;
}

template<class TYPE>
Integer BandMatrix<TYPE>::operator!= (const BandMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
         for (k = 0; k < numentry; k++)
              if (p [0][k] != x.p[0][k]) return 1;
 
        return 0;
}

template<class TYPE>
BandMatrix<TYPE> BandMatrix<TYPE>::operator* (const TYPE &x) const
{       if (!row) Error("undefined BandMatrix",__FILE__,__LINE__);
 
        Integer k;
 
        BandMatrix<TYPE> z(row,li,ri);
 
        for ( k = 0; k < numentry; k++)  z[0][k] = p[0][k]*x;
 
        return z;
}

template<class TYPE>
BandMatrix<TYPE> &BandMatrix<TYPE>::operator*= (const TYPE &x)
{
         if (!row ) Error("undefined BandMatrix",__FILE__,__LINE__);
 
        Integer i;
        for (i = 0; i < numentry; i++)
                p [0][i] *= x;
 
        return *this;
}

template<class TYPE>
Vector<TYPE> BandMatrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
        if (!row) Error("undefined BandMatrix",__FILE__,__LINE__);
 
        if (!x.n) Error("undefined Vector",__FILE__,__LINE__);
 
        if (row != x.n) Error("incompatible dimension",__FILE__,__LINE__);
 
        Vector<TYPE> z(row);
 
        Integer k,kk;

        for (k=0; k < ri+1 ; k++) z.p[0]+=x[k]*p[0][k];
        for ( k = 1; k < row-1; k++)
            for (kk=-li; kk < ri+1; kk++)
                z.p[k]+=p[k][kk]*x[k+kk];
        for (k=-li; k < 1 ; k++) z.p[row-1]+=x[row-1+k]*p[row-1][k];

        return z;
}

template<class TYPE>
Matrix<TYPE> BandMatrix<TYPE>::operator*(const Matrix<TYPE> &x) const 
{
        if (!row || !x.col || !x.row)
        Error("undefined BandMatrix or Matrix",__FILE__,__LINE__);
 
        if (row != x.row)
        Error("incompatible dimension",__FILE__,__LINE__);
 
        TYPE    a;
        Matrix<TYPE>  z (row, x.col);
        Integer i,j,k;
         
        for (i=0; i < x.col; i++)
            { 
              a=p[0][0]*x.p[0][i];
              for (j=1; j<ri+1; j++)
              a+=p[0][j]*x.p[j][i];
              z.p[0][i]=a;
            }

        for (k = 1; k < row-1; k++)
        for (i = 0; i < x.col; i++)
          {
            a=p[k][-li]*x.p[k-li][i];
           for (j = -li+1 ; j < ri+1; j++)
	      a+=p[k][j]*x.p[k+j][i];
            z.p[k][i]=a;
          }
       
          for (i = 0; i < x.col; i++) 
          for (j = -li; j < 1 ; j++) 
                 z.p[row-1][i]+=x.p[row-1+j][i]*p[row-1][j];
               
        return z;
}

template<class TYPE>
BandMatrix<TYPE> BandMatrix<TYPE>::part (const Integer i) const
{
         if (!row ) Error("undefined BandMatrix",__FILE__,__LINE__);
 
        if (i < 0 || i >= row)
                      Error("invalid index",__FILE__,__LINE__);
 
        BandMatrix  z(i+1,li,ri);
 
        for (Integer k = 0; k < z.numentry ; k++)
                        z [0][k] = p [0][k];
 
        return z;
}

template<class S>
S Spur(const BandMatrix<S> & x)
{
  S a;
  Integer i;
  for (i=0; i<x.row; i++) a+=x.p[i][0];
  return a;
}

#ifdef __GNUC__
template Integer Spur<Integer>(const BandMatrix<Integer> &);
template Double Spur<Double>(const BandMatrix<Double> &);
#endif

template<class S>
std::ostream & operator << (std::ostream & out, const BandMatrix<S> &mat)
{

Integer i,j;
Integer n=mat.GetRow();
Integer li=mat.GetLi();
Integer ri=mat.GetRi();

for (i=0; i<=ri; i++)  out << mat[0][i] << " "; 
for (i=ri+1; i <n; i++) out << 0 << " ";
                        out << endl;

for (i=1; i < n-1; i++)
{

  for (j=0; j < i-li; j++) out << 0 << " ";
  for (j=i-li; j < i+ri+1 ; j++) out << mat[i][-i+j] << " ";
  for (j=i+ri+1; j < n; j++) out << 0 <<" ";
 
  out << endl;
}

for (i=0; i< n-li-1; i++) out << 0 << " ";
for (i=-li; i<1; i++)  out << mat[n-1][i] << " ";
                       out << endl;
 
 return out;
}

#ifdef __GNUC__ 
template std::ostream & operator<<<Integer> (std::ostream & , const BandMatrix<Integer>&);
template std::ostream & operator<<<Double> (std::ostream & , const BandMatrix<Double>& );
#endif

} // end of namespace


