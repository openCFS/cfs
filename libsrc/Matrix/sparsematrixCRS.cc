#include <iostream>
#include <fstream>
#include <time.h>
#include <string>

//#include <general_head.hh> 
//#include <utils_head.hh>

#include "sparsematrixCRS.hh"

namespace CoupledField
{

template<class TYPE>
SparseMatrix<TYPE> :: SparseMatrix()
{
 numentry = 0;
 row = 0;
 col = 0;
 p = NULL;
 pc = NULL;
 pf = NULL;
}

template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix (const Integer anumentry, const Integer arow, const Integer acol)
{
#ifdef TRACE
    (*trace) << "entering SparseMatrix::SparseMatrix" << std::endl;
#endif
        numentry = anumentry;
        row = arow;
        col = acol;

        p = new TYPE[numentry];
        pc = new Integer[numentry];
        pf = new Integer[row+1];
        
       if (!p || !pc || !pf) Error("out of memory",__FILE__,__LINE__);
}

template<class TYPE>
void SparseMatrix<TYPE>::Init ()
{
  Integer i;
  for (i=0; i<numentry; i++) { p[i]=0; pc[i]=0;}
}

template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix (const SparseMatrix<TYPE> &x)
{      
        if (!x.row )
        Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        row = x.row;
        col = x.col;
        numentry=x.numentry;

        if (!p) { delete [] p;}
        if (!pc) { delete [] pc;}
        if (!pf) {delete [] pf;}
           
        p = new TYPE[numentry];
        pc = new Integer[numentry];
        pf = new Integer[row+1];
        
        if (!p || !pc || !pf) Error("out of memory", __FILE__, __LINE__);
  
        Integer k;
        for (k=0; k < numentry; k++) p[k]=x.p[k];
	for (k=0; k < numentry; k++) pc[k]=x.pc[k];
        for (k=0; k < row+1; k++) pf[k]=x.pf[k];
}

template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix (const Matrix<TYPE> & x)
{
  row=x.row;
  col=x.col;
  numentry=0;

  pf=new Integer[x.row+1];
 
  SparseMatrix<TYPE> help(x.row*x.col,x.row,x.col);

  Integer k,kk;

  for (k=0; k<row; k++)
    { pf[k]=numentry;
    for (kk=0; kk<col; kk++)
         if (x.p[k][kk]!=0) 
           { help.p[numentry]=x.p[k][kk]; 
	     help.pc[numentry]=kk;
             numentry++;
           };
    }
 
  pf[row]=numentry;

  p=new TYPE[numentry];
  pc=new Integer[numentry];
  
  for (k=0; k<numentry; k++)
    {
     p[k]=help.p[k];
     pc[k]=help.pc[k];
    }
}

template<class TYPE>
void SparseMatrix<TYPE> :: Resize(const Integer anumentry, const Integer arow, const Integer acol)
{
  if (p) delete [] p;
  if (pc) delete [] pc;
  if (pf) delete [] pf;
 
   row=arow;
   col=acol;
   numentry=anumentry;
 
   p = new TYPE [numentry];
   pc=new Integer[numentry];
   pf=new Integer[row];

   if (!p || !pc |!pf) Error("out of memory",__FILE__,__LINE__);
}

template<class TYPE>
TYPE * SparseMatrix<TYPE>::operator[] (const Integer i) const
{
        if (!p) Error("undefined SparseMatrix",__FILE__,__LINE__);
        if (i < 0 || i >= row) Error("invalid index",__FILE__,__LINE__);
 
        return p+pf[i];
}

template<class TYPE>
SparseMatrix<TYPE> & SparseMatrix<TYPE>::operator=(const SparseMatrix<TYPE> &x)
{

  if (!x.row) Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        if (this == &x)  return *this;
 
        Integer k;
 
        if (numentry!=x.numentry )
	  { if (!p) delete [] p;
	    if (!pc) delete [] pc;

            numentry=x.numentry;
             
            p=new TYPE[numentry];
            pc=new Integer[numentry];            
                       
           if (!p || !pc) Error("out of memory",__FILE__,__LINE__);
          }

        if (row!=x.row )
        {
          delete [] pf;

          row=x.row;

          pf=new Integer[row+1];
          if (!pf) Error("out of memory",__FILE__,__LINE__);
           
        }
 
        if (col!=x.col) col=x.col;

        for (k = 0; k < numentry; k++)
	  { p[k]=x.p[k]; pc[k]=x.pc[k];}
        for (k = 0; k < row+1; k++)
	   pf[k]=x.pf[k];
 
        return *this;
}

template<class TYPE>
SparseMatrix<TYPE> & SparseMatrix<TYPE>::operator+ () 
{
        if (!row) Error("undefined SparseMatrix",__FILE__,__LINE__);
        return *this;
}

template<class TYPE>
SparseMatrix<TYPE> & SparseMatrix<TYPE>::operator+=(const SparseMatrix<TYPE> &x)
{
  if (row != x.row || col!=x.col ) Error("Bad matrix sizes",__FILE__,__LINE__);

  SparseMatrix<TYPE> z(*this);
  *this=z+x;
  return *this;
}

template<class TYPE>
SparseMatrix<TYPE> SparseMatrix<TYPE>::operator+(const SparseMatrix<TYPE> &x) const 
{
        if (!row)
            Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        if (row != x.row || col!=x.col )
             Error("incompatible dimension for +",__FILE__,__LINE__);
        
        SparseMatrix<TYPE> z(x.numentry + numentry,x.row,x.col);
        z.Init();

        Integer i,j,l,n,apos,bpos;

        TYPE a;      
        z.pf[0]=0; l=0;
        
        for (i=0; i<row; i++)
	  {
           apos=0; bpos=0;
            
            if ( x.pf[i+1] == x.pf[i] || pf[i+1]==pf[i]) {

               if (pf[i+1]!=pf[i]) 
                     { n=pf[i+1]-pf[i];
		     for (j=0; j<n; j++) { 
                                            z.p[l]+=p[j+pf[i]];
                                            z.pc[l]=pc[j+pf[i]];
                                            l++;}
                     }

                if (x.pf[i+1]!= x.pf[i]) 
                     { n=x.pf[i+1]-x.pf[i];
		     for (j=0; j<n; j++) { 
                                           z.p[l]+=x.p[j+x.pf[i]];
                                           z.pc[l]=x.pc[j+x.pf[i]];
                                           l++;}
                     }
 
	       z.pf[i+1]=l;
               continue;}

	    n=x.pf[i+1]-x.pf[i]+pf[i+1]-pf[i];
       
            for (j=0; j<n; j++)
	      {
     if (x.pc[x.pf[i]+apos]<pc[pf[i]+bpos])
       { if (j!=0 && x.pc[x.pf[i]+apos]==z.pc[l-1])
                    { z.p[l-1]+=x.p[x.pf[i]+apos];continue;}
                else {
                      z.pc[l]=x.pc[apos+x.pf[i]];
                      z.p[l]+=x.p[apos+x.pf[i]];
                      l++;}
	              apos++;
                    }
      else {
        if (x.pc[x.pf[i]+apos] == pc[pf[i]+bpos]) 
	  {       z.pc[l]=pc[bpos+pf[i]];
                  a=x.p[apos+x.pf[i]]+p[bpos+pf[i]];
                  z.p[l]=a;
		 
                  apos++; bpos++; l++; j++;
                  }
	else {
	      if (j!=0 && pc[pf[i]+bpos]==z.pc[l-1])
                     { z.p[l-1]+=p[pf[i]+bpos]; continue;}
	             else {
                           z.pc[l]=pc[bpos+pf[i]];
                           z.p[l]+=p[bpos+pf[i]];
                           l++;}
	      bpos++;
	     }
	     }
         } // end for j
	      z.pf[i+1]=l;
          }
      
      SparseMatrix<TYPE> result(l,x.row,x.col);
      for (i=0; i<l; i++)
	   { result.p[i]=z.p[i];
             result.pc[i]=z.pc[i];
           }
      for (i=0; i < x.row+1; i++)
	     result.pf[i]=z.pf[i];
 
       return result;
}
template<class TYPE>
SparseMatrix<TYPE> & SparseMatrix<TYPE>::operator-=(const SparseMatrix<TYPE> &x)
{
  if (row != x.row || col!=x.col ) Error("Bad matrix sizes",__FILE__,__LINE__);

  SparseMatrix<TYPE> z(*this);
  *this=z-x;
   return *this;
}

template<class TYPE>
SparseMatrix<TYPE> SparseMatrix<TYPE>::operator-(const SparseMatrix<TYPE> &x) const 
{
        if (!row)
            Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        if (row != x.row || col!=x.col )
             Error("incompatible dimension for +",__FILE__,__LINE__);
        
        SparseMatrix<TYPE> z(x.numentry + numentry,x.row,x.col);
        z.Init();

        Integer i,j,l,n,apos,bpos;

        TYPE a;      
        z.pf[0]=0; l=0;
        
        for (i=0; i<row; i++)
	  {
           apos=0; bpos=0;
            
            if ( x.pf[i+1] == x.pf[i] || pf[i+1]==pf[i]) {

               if (pf[i+1]!=pf[i]) 
                     { n=pf[i+1]-pf[i];
		     for (j=0; j<n; j++) { 
                                            z.p[l]+=p[j+pf[i]];
                                            z.pc[l]=pc[j+pf[i]];
                                            l++;}
                     }

                if (x.pf[i+1]!= x.pf[i]) 
                     { n=x.pf[i+1]-x.pf[i];
		     for (j=0; j<n; j++) { 
                                           z.p[l]-=x.p[j+x.pf[i]];
                                           z.pc[l]=x.pc[j+x.pf[i]];
                                           l++;}
                     }
 
	       z.pf[i+1]=l;
               continue;
                                                        }

	    n=x.pf[i+1]-x.pf[i]+pf[i+1]-pf[i];
       
            for (j=0; j<n; j++)
	      {
     if (x.pc[x.pf[i]+apos]<pc[pf[i]+bpos])
       { if (j!=0 && x.pc[x.pf[i]+apos]==z.pc[l-1])
                    { z.p[l-1]-=x.p[x.pf[i]+apos];continue;}
                else {
                      z.pc[l]=x.pc[apos+x.pf[i]];
                      z.p[l]-=x.p[apos+x.pf[i]];
                      if (z.p[l] != 0) l++;}
	              apos++;
                    }
      else {
        if (x.pc[x.pf[i]+apos] == pc[pf[i]+bpos]) 
	  {       z.pc[l]=pc[bpos+pf[i]];
                  a=-x.p[apos+x.pf[i]]+p[bpos+pf[i]];
                  z.p[l]=a;
		 
                  apos++; bpos++;
                  if (z.p[l] != 0) l++; j++;
                  }
	else {
	      if (j!=0 && pc[pf[i]+bpos]==z.pc[l-1])
                     { z.p[l-1]+=p[pf[i]+bpos]; continue;}
	             else {
		           z.pc[l]=pc[bpos+pf[i]];
                           z.p[l]+=p[bpos+pf[i]];
                           if (z.p[l] != 0) l++;}
	      bpos++;
	     }
	     }
         } // end for j
	      z.pf[i+1]=l;
          }
       
      SparseMatrix<TYPE> result(l,x.row,x.col);
      for (i=0; i<l; i++)
	   { result.p[i]=z.p[i];
             result.pc[i]=z.pc[i];
           }
      for (i=0; i < x.row+1; i++)
	     result.pf[i]=z.pf[i];
 
       return result;
}

template<class TYPE>
SparseMatrix<TYPE> SparseMatrix<TYPE>::operator-() const
{
    SparseMatrix<TYPE> z(numentry,row,col);

	Integer k;
         for ( k = 0; k < numentry; k++)
	   { z.p[k] = -p[k]; z.pc[k]=pc[k];}

         for (k=0; k< row+1; k++) z.pf[k]=pf[k]; 
        return z;
}

template<class TYPE>
Integer SparseMatrix<TYPE>::operator== (const SparseMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
 
         for (k = 0; k < numentry; k++)
              if (p[k]!= x.p[k] || pc[k]!=x.pc[k]) return 0;
 
        return 1;
}

template<class TYPE>
Integer SparseMatrix<TYPE>::operator!= (const SparseMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
         for (k = 0; k < numentry; k++)
              if (p[k] != x.p[k] || pc[k]!=x.pc[k]) return 1;
 
        return 0;
}

template<class TYPE>
SparseMatrix<TYPE> SparseMatrix<TYPE>::operator* (const TYPE &x) const
{       if (!row) Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        Integer k;
 
        SparseMatrix<TYPE> z(numentry,row,col);
 
        for ( k = 0; k < numentry; k++)
                { z.p[k] = p[k]*x;  z.pc[k]= pc[k];} 
        for (k=0; k< row + 1; k++) z.pf[k]=pf[k];
 
        return z;
}

template<class TYPE>
SparseMatrix<TYPE> &SparseMatrix<TYPE>::operator*= (const TYPE &x)
{
         if (!row ) Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        Integer i;
        for (i = 0; i < numentry; i++)  p[i] *= x;
 
        return *this;
}

template<class TYPE>
Vector<TYPE> SparseMatrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
        if (!row) Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        if (!x.n) Error("undefined Vector",__FILE__,__LINE__);
 
        if (col != x.n) Error("incompatible dimension",__FILE__,__LINE__);
 
        Vector<TYPE> z(row); // should be initialized by zero
 
        Integer k,kk;

        for (k=0; k < row ; k++) 
            for (kk=pf[k]; kk<pf[k+1]; kk ++)
                 z[k]+=p[kk]*x[pc[kk]];

        return z;
}

template<class TYPE>
Matrix<TYPE> SparseMatrix<TYPE>::operator*(const Matrix<TYPE> &x) const
{
        if (!numentry || !x.col || !x.row)
        Error("undefined SparseMatrix or Matrix",__FILE__,__LINE__);
 
        if (col != x.row)
        Error("incompatible dimension for multiplycation",__FILE__,__LINE__);
 
        Matrix<TYPE>  z (row, x.col); // should be initialized by zero
        Integer i,j,k;
      
	for (i=0; i<row; i++)
	  for (j=0; j<x.col; j++)
            for (k=pf[i]; k < pf[i+1]; k++)
	      z[i][j]+=p[k]*x.p[pc[k]][j];
   
         return z;
}

template<class TYPE>
SparseMatrix<TYPE> SparseMatrix<TYPE>::part (const Integer index) const
{
         if (!row ) Error("undefined SparseMatrix",__FILE__,__LINE__);
 
        if (index < 0 || index >= row)
                      Error("invalid index",__FILE__,__LINE__);
         
        SparseMatrix  help(numentry,index,index);
 
        Integer i,k,l;
        
        l=0;
        help.pf[0]=0;
        for (i=0; i<index; i++)
	  { for (k=pf[i]; k < pf[i+1]; k++)
	    { if (pc[k]<index)
                              { help.p[l]=p[k];
      	                        help.pc[l]=pc[k];
                                l+=1;}
	    else k=pf[i+1];                  
            }
	  help.pf[i+1]=l;
          }

        SparseMatrix z(l,index,index);
  
        for (i=0; i<l; i++)
           { z.p[i]=help.p[i];
	   z.pc[i]=help.pc[i];}
        for (i=0; i<index+1; i++) z.pf[i]=help.pf[i];
         
        return z;
}

template< class TYPE>
void SparseMatrix<TYPE>::precond(Vector<TYPE> & e, const Vector<TYPE> r, enum precond type)
{
#ifdef TRACE
  (*trace) << "entering SparseMatrix::precond" << std::endl;
#endif
  if (!e.size()) Error("Define vector before use of precond");

  Integer n=r.size(); 
  if (row!=r.size()) Error("Incompatible dimension in precond of SparseMatrix"); 
 
 Integer i,j;
 Double omega;
 Vector<Double> diag(row);
 
  switch(type)
{
  case non: e=r; break;

  case Jacobi:
       for (i=0; i< row; i++) 
      {  for  (j=pf[i]; j < pf[i+1]; j++)
          diag[i]+=p[j]*p[j];
        //  idiag[i]=sqrt(diag[i]);
     
          e[i]=r[i]/diag[i];  ///// May be (*this)
      }
       break;

  case SSOR:
       omega=1.2; /// omega=1.2
       TYPE sum;
 
          for (i=0; i< row; i++)
      {  for  (j=pf[i]; j < pf[i+1]; j++)
          { if (pc[j]==i) diag[i]+=p[j]*p[j];
            else diag[i]+=omega*omega*p[j]*p[j];
          }
 //         diag[i]=sqrt(diag[i]);  // ####################################
      }

       for (i=0; i<n; i++)
       {
         sum=0;
         for (j=pf[i]; j<pf[i+1]; j++) 
         { if (pc[j]<i) sum+=p[j]*e[pc[j]]; }
         e[i]=(r[i]-omega*sum)/diag[i];
       }
 
       for (i=n-1; i>=0; i--)
       {
         sum=0;
         for (j=pf[i]; j<pf[i+1]; j++)
         { if (pc[j]>i) sum+=p[j]*e[pc[j]];}
         e[i]-=omega*sum/diag[i];
       }
 
       e*=omega*(2-omega);
       break;
default:
       Error("Wrong type of precondition",__FILE__,__LINE__);
 
} //end of switch
}

template<class TYPE>
TYPE Spur (const SparseMatrix<TYPE> &x)
{
if (!x.row) Error("undefined SparseMatrix",__FILE__,__LINE__);
 
if (x.row != x.col)
  Error("incompatible dimension for calculation trace",__FILE__,__LINE__);

        TYPE a=0;
        Integer i,j,n;
        for ( i = 0; i < x.row; i++)
	  {   n=x.pf[i+1]-x.pf[i];
                   for (j=0; j<n; j++)
                   if (x.pc[j+x.pf[i]]==i) a+=x.p[j+x.pf[i]];}
 
        return a;
}

template Double Spur(const SparseMatrix<Double> &);
template Integer Spur(const SparseMatrix<Integer> &);

template<class S>
std::ostream & operator << (std::ostream & out, const SparseMatrix<S> &m)
{

Integer i,j,k;

for (i=0; i<m.row; i++)
  {
          // If we have zero row
       if (m.pf[i]==m.pf[i+1]) 
       { for (k=0; k<m.col; k++) out << 0 << " ";
         out << std::endl;
         continue;
       }
 
       for (k=0; k<m.pc[m.pf[i]]; k++) out << 0 << " ";

       for (j=m.pf[i]; j<m.pf[i+1] ; j++)
       {
                  out << m.p[j] << " ";
        for (k=m.pc[j]+1; k<m.pc[j+1] && j+1 < m.numentry; k++)
                  out << 0 << " ";   
       }

       for (k=m.pc[j-1]+1; k<m.col; k++) out << 0 << " ";

   out << std::endl;

  }
  return out;
}
 
template std::ostream & operator<<<Integer> (std::ostream & , const SparseMatrix<Integer>&);
template std::ostream & operator<<<Double> (std::ostream & , const SparseMatrix<Double>& );


} // end of namespace


