#include <iostream.h>
#include <fstream.h>
#include <time.h>
#include <string>
 
#include <general_head.hh>
#include <utils_head.hh>

namespace CoupledField
{

template<class TYPE>
SymBandMatrix<TYPE> :: SymBandMatrix()
{
 row=0;
 li=0;
 p=NULL;
 numentry=0;
}

template<class TYPE>
SymBandMatrix<TYPE>:: SymBandMatrix (const Integer arow, const Integer alcln)
{
#ifdef TRACE
    cout << "entering SymBandMatrix::SymBandMatrix" << endl;
#endif
        row = arow;
        li=alcln;
        
        numentry=(row-1)*(li+1)+1;

        p = new TYPE* [row];
 
        if (!p) Error("out of memory",__FILE__,__LINE__);
 
        p[0]=new TYPE[numentry];

        if (!p[0]) Error("out of memory",__FILE__,__LINE__);
 
        Integer k;

        for (k=1; k < row; k++) p[k]=p[k-1]+li+1;
}

template<class TYPE>
SymBandMatrix<TYPE>::SymBandMatrix (const SymBandMatrix<TYPE> &x)
{       if (!x.row )  Error("undefined SymBandMatrix");
 
        row = x.row;
        li = x.li;
        numentry=x.numentry;
 
        if (!p) { delete p[0];
                  delete [] p; }
 
        p = new TYPE * [row];
        if (!p) Error("out of memory", __FILE__, __LINE__);
 
        p[0]=new TYPE[numentry];
        if (!p[0]) Error("out of memory", __FILE__, __LINE__);
 
        Integer k;
 
        for (k=0; k < numentry; k++)  p[0][k]=x.p[0][k];
        for (k=1; k < row; k++) p[k]=p[k-1]+li+1;
}

template<class TYPE>
void SymBandMatrix<TYPE> :: Resize(const Integer irow, const Integer ili)
{
   Integer k;
 
   if (p) {delete p[0];
           delete [] p;}
 
   row=irow;
   li=ili;
  
   numentry=(row-1)*(li+1)+1;
 
   p = new TYPE* [row];
   p[0]=new TYPE[numentry];

   if (!p) Error("out of memory",__FILE__,__LINE__);
   if (!p[0])  Error("out of memory",__FILE__,__LINE__);
 
  for (k=1; k < row; k++) p[k]=p[k-1]+li+1;
}

template<class TYPE>
void SymBandMatrix<TYPE>::cut(const Integer ai, const Integer aj)
{
if (ai < 0 || ai > row-1 || aj < 0 || aj > li)
         Error("invalid index in function cut",__FILE__,__LINE__);
if (!row) Error("SymBandMatrix is undefined in function cut");

Error("Not implemented cut for SymBandMatrix");
}

template<class TYPE>
TYPE * SymBandMatrix<TYPE>::operator[] (const Integer i) const
{
        if (!p) Error("undefined SymBandMatrix",__FILE__,__LINE__);
        if (i < 0 || i >= row) Error("invalid index",__FILE__,__LINE__);
 
        return p[i];
}

template<class TYPE>
TYPE & SymBandMatrix<TYPE>::At (const Integer i, const Integer j) 
{
    if (!p) Error("undefined SymBandMatrix",__FILE__,__LINE__);
     if (i < 0 || i >= row) Error("invalid index for row",__FILE__,__LINE__); 
     if (j < 0 || j >= row) Error("invalid index for col",__FILE__,__LINE__); 
      if (j< (i-li) || j > (i+li))
 Error("index in operator[][] is out of band",__FILE__,__LINE__);     
      else return p[i][-abs(i-j)];
}

template<class TYPE>
SymBandMatrix<TYPE> & SymBandMatrix<TYPE>::operator=(const SymBandMatrix<TYPE> &x)
{
 cout << "entering Copy" << endl;
       if (!x.row) Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        if (this == &x)  return *this;
 
        Integer k;
 
        if (row!=x.row )
        {
           if (p)   { 
                      delete [] p[0];
                      delete [] p;}
 
           row=x.row;
           li=x.li;
           
           numentry=x.numentry;
 
           p = new TYPE* [row];
           p[0] = new TYPE[numentry];
 
           if (!p) Error("out of memory",__FILE__,__LINE__);
           if (!p[0]) Error("out of memory",__FILE__,__LINE__);

           for (k=1; k < row; k++) p[k]=p[k-1]+li+1;
        }

        if (li!=x.li)
        {
          delete [] p[0];

          li=x.li;
          numentry=x.numentry;

          p[0]=new TYPE[numentry];
           
          for (k=1; k < row; k++) p[k]=p[k-1]+li+1;
        }
 
        for (k = 0; k < numentry; k++)
           p[0][k]=x.p[0][k];
 
        return *this;
}

template<class TYPE>
SymBandMatrix<TYPE> & SymBandMatrix<TYPE>::operator+ () 
{
        if (!row) Error("undefined SymBandMatrix",__FILE__,__LINE__);
        return *this;
}

template<class TYPE>
SymBandMatrix<TYPE> & SymBandMatrix<TYPE>::operator+=(const SymBandMatrix<TYPE> &x)
{
  if (row != x.row || li!=x.li) Error("Bad matrix sizes");
  for (Integer i=0; i<numentry; i++)
       p[0][i]+=x[0][i];
  return *this;
}

template<class TYPE>
SymBandMatrix<TYPE> SymBandMatrix<TYPE>::operator+(const SymBandMatrix<TYPE> &x) const 
{
        if (!row)
            Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        if (row != x.row || li!=x.li )
             Error("incompatible dimension",__FILE__,__LINE__);
        
        SymBandMatrix<TYPE> z(row, li);
        Integer k;
        for (k=0; k< numentry; k++)
         z[0][k]=x.p[0][k]+p[0][k];
        return z;
}

template<class TYPE>
SymBandMatrix<TYPE> SymBandMatrix<TYPE>::operator-() const
{
        if (!row )
            Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        SymBandMatrix<TYPE> z(row,li);
 
         Integer k;
         for ( k = 0; k < numentry; k++)
              if (p[0][k]!=0)  z[0][k] = -p [0][k];
 
        return z;
}

template<class TYPE>
SymBandMatrix<TYPE> SymBandMatrix<TYPE>::operator-(const SymBandMatrix<TYPE> &x) const
{
        if (!row || !x.row)
            Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        if (row != x.row || li!=x.li )
             Error("incompatible dimension",__FILE__,__LINE__);
 
        SymBandMatrix<TYPE> z(row,li);
 
         Integer k;
         for ( k = 0; k < numentry; k++)
                z [0][k] = -x.p [0][k]+p[0][k];
 
        return z;
}

template<class TYPE>
SymBandMatrix<TYPE> &SymBandMatrix<TYPE>::operator-=(const SymBandMatrix<TYPE> &x)
{
        if (!row || !x.row)
            Error("Don't use -= to undefuned matrix",__FILE__,__LINE__);
 
         if (row != x.row || li!=x.li)
             Error("incompatible dimension for -",__FILE__,__LINE__);
 
        Integer k;
        for ( k = 0; k < numentry; k++)
                p [0][k] -= x.p [0][k];
 
        return *this;
}

template<class TYPE>
Integer SymBandMatrix<TYPE>::operator== (const SymBandMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
 
         for (k = 0; k < numentry; k++)
              if (p [0][k] != x.p[0][k]) return 0;
 
        return 1;
}

template<class TYPE>
Integer SymBandMatrix<TYPE>::operator!= (const SymBandMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
         for (k = 0; k < numentry; k++)
              if (p [0][k] != x.p[0][k]) return 1;
 
        return 0;
}

template<class TYPE>
SymBandMatrix<TYPE> SymBandMatrix<TYPE>::operator* (const TYPE &x) const
{       if (!row) Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        Integer k;
 
        SymBandMatrix<TYPE> z(row,li);
 
        for ( k = 0; k < numentry; k++)
                z[0][k] = p[0][k]*x;
 
        return z;
}

template<class TYPE>
SymBandMatrix<TYPE> & SymBandMatrix<TYPE>::operator*= (const TYPE &x)
{
         if (!row ) Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        Integer i;
        for (i = 0; i < numentry; i++)
                p [0][i] *= x;
 
        return *this;
}

template<class TYPE>
Vector<TYPE> SymBandMatrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
        if (!row) Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        if (!x.n) Error("undefined Vector",__FILE__,__LINE__);
 
        if (row != x.n) Error("incompatible dimension",__FILE__,__LINE__);
 
        Vector<TYPE> z(row);
 
        Integer k,kk;

       for ( k = 0; k < row; k++)
         {
           if (k!=0)  for (kk=-li; kk < 0; kk++)
                        z.p[k]+=p[k][kk]*x[k+kk];

                        z.p[k]+=p[k][0]*x[k];

           if (k!=row-1)  for (kk=1; kk < li+1; kk++)
                            z.p[k]+=p[kk+k][-kk]*x[k+kk];
         }
        
        return z;
}

template<class TYPE>
Matrix<TYPE> SymBandMatrix<TYPE>::operator*(const Matrix<TYPE> &x) const 
{
        if (!row || !x.col || !x.row)
        Error("undefined SymBandMatrix or Matrix",__FILE__,__LINE__);
 
        if (row != x.row)
        Error("incompatible dimension",__FILE__,__LINE__);
 
        Matrix<TYPE>  z (row, x.col);
        
        Integer i,j,k;
        
        for (k = 0; k < row; k++)
        for (i = 0; i < x.col; i++)
          {
           if (k!=0) { for (j = -li ; j < 0; j++)
	               z.p[k][i]+=p[k][j]*x.p[k+j][i];}

                      z.p[k][i]+=p[k][0]*x.p[k][i];

           if (k!=row-1) { for (j=1; j < 1+li; j++)
	               z.p[k][i]+=p[j+k][-j]*x.p[k+j][i];}
           }
                       
        return z;
}

template<class TYPE>
SymBandMatrix<TYPE> SymBandMatrix<TYPE>::part (const Integer i) const
{
         if (!row ) Error("undefined SymBandMatrix",__FILE__,__LINE__);
 
        if (i < 0 || i >= row)
                      Error("invalid index",__FILE__,__LINE__);
 
        SymBandMatrix  z(i+1,li);
 
        for (Integer k = 0; k < z.numentry ; k++)
                        z [0][k] = p [0][k];
 
        return z;
}

template< class TYPE>
void SymBandMatrix<TYPE>::precond(Vector<TYPE> & e, const Vector<TYPE> r, const Integer type)
{
#ifdef TRACE
  (*trace) << "entering SymMatrix::precond" << endl;
#endif
 
 Integer i,j;
 Integer n=r.size();
 Double omega;
 
 switch(type)
{ case 0: e=r; break;
  case 1:
       for (i=0; i<n; i++) e[i]=r[i]/p[i][0];  ///// May be (*this)
       break;
  case 2:
       omega=1.2; /// omega=1.2
       TYPE sum;
       Error("Not implemented",__FILE__,__LINE__);
       for (i=0; i<n; i++)
       {
         sum=0;
         for (j=0; j<i; j++) sum+=p[j][i]*e[j];
         e[i]=(r[i]-omega*sum)/p[i][0];
       }
 
       for (i=n-1; i>=0; i--)
       {
         sum=0;
         for (j=i+1; j<n; j++) sum+=p[i][j-i]*e[j];
         e[i]-=omega*sum/p[i][0];
       }
 
       e*=omega*(2-omega);
       break;
 
default:
       Error("Wrong type of Matrix preconditioning",__FILE__,__LINE__);
 
} //end of switch
}
template<class S>
S Spur(const SymBandMatrix<S> & x)
{
  S a;
  Integer i;
  for (i=0; i<x.row; i++) a+=x.p[i][0];
  return a;
}

template Integer Spur<Integer>(const SymBandMatrix<Integer> &);
template Double Spur<Double>(const SymBandMatrix<Double> &);

template<class S>
ostream & operator << (ostream & out, const SymBandMatrix<S> &mat)
{

out.setf(ios::scientific);

if (mat.row==0) out << "Null matrix" << endl;
Integer i,j;

for (i=0; i < mat.row; i++)
{
if (i!=0) { 
              for (j=0; j < i-mat.li; j++) out << 0 << " ";
              for (j=-mat.li; j < 0 ; j++) out << mat[i][j] << " ";
          }
             
              out << mat[i][0] <<" ";
 
if ( i!=mat.row-1 ) 
            {
              for (j=1; j < mat.li+1; j++) out << mat[j+i][-j] << " ";
              for (j=i+mat.li+1; j<mat.row; j++) out << 0 << " ";
             } 

               out << endl;
}

 out.setf(0, ios::floatfield);
 return out;

}
 
template ostream & operator<<<Integer> (ostream & , const SymBandMatrix<Integer>&);
template ostream & operator<<<Double> (ostream & , const SymBandMatrix<Double>& );


} // end of namespace


