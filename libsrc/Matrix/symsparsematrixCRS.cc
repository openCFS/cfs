#include <iostream>
#include <fstream>
//#include <time.h>
//#include <string>

#include <general_head.hh> 
#include <utils_head.hh>

namespace CoupledField
{

template<class TYPE>
SymSparseMatrix<TYPE> :: SymSparseMatrix()
{
 numentry = 0;
 row = 0;

 p = NULL;
 pc = NULL;
 pf = NULL;
}

template<class TYPE>
SymSparseMatrix<TYPE>:: SymSparseMatrix (const Integer anumentry, const Integer arow)
{
#ifdef TRACE
    (*trace) << "entering SymSparseMatrix:: SymSparseMatrix" << std::endl;
#endif
        numentry = anumentry;
        row = arow;
        
        p = new TYPE[numentry];
        pc = new Integer[numentry];
        pf = new Integer[row+1];
        
       if (!p || !pc || !pf) Error("out of memory",__FILE__,__LINE__);
}

template<class TYPE>
void SymSparseMatrix<TYPE>::Init ()
{
  Integer i;
  for (i=0; i<numentry; i++) { p[i]=0; pc[i]=0;}
}

template<class TYPE>
SymSparseMatrix<TYPE>:: SymSparseMatrix (const SymSparseMatrix<TYPE> &x)
{      
        if (!x.row )
        Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        row = x.row;
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
SymSparseMatrix<TYPE>:: SymSparseMatrix (const SymMatrix<TYPE> & x)
{
  row=x.size;
  numentry=0;

  pf=new Integer[x.size+1];
 
  SymSparseMatrix<TYPE> help(x.size*(x.size+1)/2, x.size);

  Integer k,kk;

  for (k=0; k<row; k++)
    { pf[k]=numentry;
    for (kk=0; kk<k+1; kk++)
         if (x.p[k][kk]!=0) 
           { help.p[numentry]=x.p[k][kk]; 
	     help.pc[numentry]=kk;
             numentry++;
           }
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
void SymSparseMatrix<TYPE> :: Resize(const Integer anumentry, const Integer arow)
{
  if (p) delete [] p;
  if (pc) delete [] pc;
  if (pf) delete [] pf;
 
   row=arow;
   numentry=anumentry;
 
   p = new TYPE [numentry];
   pc=new Integer[numentry];
   pf=new Integer[row];

   if (!p || !pc |!pf) Error("out of memory",__FILE__,__LINE__);
}

template<class TYPE>
TYPE * SymSparseMatrix<TYPE>::operator[] (const Integer i) const
{
        if (!p) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
        if (i < 0 || i >= row) Error("invalid index",__FILE__,__LINE__);
 
        return p+pf[i];
}

template<class TYPE>
SymSparseMatrix<TYPE> &SymSparseMatrix<TYPE>::operator=(const SymSparseMatrix<TYPE> &x)
{

  if (!x.row) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
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
 
        for (k = 0; k < numentry; k++)  { p[k]=x.p[k]; pc[k]=x.pc[k];} 
        for (k = 0; k < row+1; k++)  pf[k]=x.pf[k];
 
        return *this;
}

template<class TYPE>
SymSparseMatrix<TYPE> &SymSparseMatrix<TYPE>::operator+ () 
{
        if (!row) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
        return *this;
}

template<class TYPE>
SymSparseMatrix<TYPE> &SymSparseMatrix<TYPE>::operator+=( const SymSparseMatrix<TYPE> &x)
{
  if (row != x.row ) Error("Bad matrix sizes",__FILE__,__LINE__);

  SymSparseMatrix<TYPE> z(*this);
  *this=z+x;
  return *this;
}

template<class TYPE>
SymSparseMatrix<TYPE> SymSparseMatrix<TYPE>::operator+(const SymSparseMatrix<TYPE> &x) const 
{
        if (!row)
            Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        if (row != x.row )
             Error("incompatible dimension for +",__FILE__,__LINE__);
        
        SymSparseMatrix<TYPE> z(x.numentry + numentry,x.row);

        Integer i,j,l,n,apos,bpos;

        z.pf[0]=0; l=0;
        for (i=0; i<x.numentry+numentry; i++) z.p[i]=0;
        
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
                  z.p[l]=x.p[apos+x.pf[i]]+p[bpos+pf[i]];
		 
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
      
      SymSparseMatrix<TYPE> result(l,x.row);
      for (i=0; i<l; i++)
	   { result.p[i]=z.p[i];
             result.pc[i]=z.pc[i];
           }
      for (i=0; i < x.row+1; i++)
	     result.pf[i]=z.pf[i];
 
       return result;
}
template<class TYPE>
SymSparseMatrix<TYPE> &SymSparseMatrix<TYPE>::operator-=(const SymSparseMatrix<TYPE> &x)
{
  if (row != x.row ) Error("Bad matrix sizes",__FILE__,__LINE__);

  SymSparseMatrix<TYPE> z(*this);
  *this=z-x;
   return *this;
}

template<class TYPE>
SymSparseMatrix<TYPE> SymSparseMatrix<TYPE>::operator-(const SymSparseMatrix<TYPE> &x) const 
{
        if (!row)
            Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        if (row != x.row )
             Error("incompatible dimension for +",__FILE__,__LINE__);
        
        SymSparseMatrix<TYPE> z(x.numentry + numentry,x.row);
        z.Init();       
 
        Integer i,j,l,n,apos,bpos;

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
                  z.p[l]=-x.p[apos+x.pf[i]]+p[bpos+pf[i]];
		 
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
       
      SymSparseMatrix<TYPE> result(l,x.row);
      for (i=0; i<l; i++)
	   { result.p[i]=z.p[i];
             result.pc[i]=z.pc[i];
           }
      for (i=0; i < x.row+1; i++)
	     result.pf[i]=z.pf[i];
 
       return result;
}

template<class TYPE>
SymSparseMatrix<TYPE> SymSparseMatrix<TYPE>::operator-() const
{
    SymSparseMatrix<TYPE> z(numentry,row);

	Integer k;
         for ( k = 0; k < numentry; k++)
	   { z.p[k] = -p[k]; z.pc[k]=pc[k];}
         for (k=0; k< row+1; k++) z.pf[k]=pf[k]; 
        return z;
}

template<class TYPE>
Integer SymSparseMatrix<TYPE>::operator== (const SymSparseMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
 
         for (k = 0; k < numentry; k++)
              if (p[k]!= x.p[k] || pc[k]!=x.pc[k]) return 0;
 
        return 1;
}

template<class TYPE>
Integer SymSparseMatrix<TYPE>::operator!= (const SymSparseMatrix<TYPE> &x) const
{
         if (!row || !x.row)
              Error("undefined matrix",__FILE__,__LINE__);
 
         Integer k;
         for (k = 0; k < numentry; k++)
              if (p[k] != x.p[k] || pc[k]!=x.pc[k]) return 1;
 
        return 0;
}

template<class TYPE>
SymSparseMatrix<TYPE> SymSparseMatrix<TYPE>::operator* (const TYPE &x) const
{       if (!row) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        Integer k;
 
        SymSparseMatrix<TYPE> z(numentry,row);
 
        for ( k = 0; k < numentry; k++)
                { z.p[k] = p[k]*x;
                  z.pc[k]= pc[k];}
        for (k=0; k< row + 1; k++) z.pf[k]=pf[k];
 
        return z;
}

template<class TYPE>
SymSparseMatrix<TYPE> &SymSparseMatrix<TYPE>::operator*= (const TYPE &x)
{
         if (!row ) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        Integer i;
        for (i = 0; i < numentry; i++)  p[i] *= x; 
 
        return *this;
}

template<class TYPE>
Vector<TYPE> SymSparseMatrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
        if (!row) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        if (!x.n) Error("undefined Vector",__FILE__,__LINE__);
 
        if (row != x.n) Error("incompatible dimension",__FILE__,__LINE__);
 
        Vector<TYPE> z(row); // should be initialized by zero
 
        Integer i, j, k, kk;

        for (k=0; k < row ; k++) 
        {
          for (kk=pf[k]; kk<pf[k+1]; kk ++)
	    z[k]+=p[kk]*x[pc[kk]];

          for (j=k+1; j < row; j++)
            { 
              for (i=pf[j]; i<pf[j+1]; i++) 
             {  if (pc[i] > k) continue;
                if (pc[i] == k)  
                      { z[k]+=p[i]*x[j]; continue; }
             }
            }
        }
        return z;
}

template<class TYPE>
Matrix<TYPE> SymSparseMatrix<TYPE>::operator*(const Matrix<TYPE> &x) const
{
        if (!numentry || !x.row)
        Error("undefined SymSparseMatrix or Matrix",__FILE__,__LINE__);
 
        if (row != x.row)
        Error("incompatible dimension for multiplication",__FILE__,__LINE__);
 
        Matrix<TYPE>  z (row, x.col); // should be initialized by zero
        Integer i,j,k,kk,l;
      
	for (k=0; k<row; k++)
	  for (l=0; l<x.col; l++)
            { for (kk=pf[k]; kk<pf[k+1]; kk ++)
                  z[k][l]+=p[kk]*x[pc[kk]][l];
 
              for (j=k+1; j < row; j++)
                {
                 for (i=pf[j]; i<pf[j+1]; i++)
                  {  if (pc[i] > k) continue;
                     if (pc[i] == k)
                      { z[k][l]+=p[i]*x[j][l]; continue; }
                  }
                }
            }

         return z;
}

template<class TYPE>
SymSparseMatrix<TYPE>SymSparseMatrix<TYPE>::part (const Integer index) const
{
         if (!row ) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        if (index < 0 || index >= row)
                      Error("invalid index",__FILE__,__LINE__);
         
        SymSparseMatrix  help(numentry,index);
 
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

        SymSparseMatrix z(l,index);
  
        for (i=0; i<l; i++)
           { z.p[i]=help.p[i]; z.pc[i]=help.pc[i];} 

        for (i=0; i<index+1; i++) z.pf[i]=help.pf[i];
         
        return z;
}


template<class TYPE>
TYPE Spur (const SymSparseMatrix<TYPE> &x)
{
if (!x.row) Error("undefined SymSparseMatrix",__FILE__,__LINE__);
 
        TYPE a=0;
        Integer i,j,n;
        for ( i = 0; i < x.row; i++)
	  {   n=x.pf[i+1]-x.pf[i];
                   for (j=0; j<n; j++)
                   if (x.pc[j+x.pf[i]]==i) a+=x.p[j+x.pf[i]];}
 
        return a;
}

template Double Spur(const SymSparseMatrix<Double> &);
template Integer Spur(const SymSparseMatrix<Integer> &);

template<class S>
std::ostream & operator << (std::ostream & out, const SymSparseMatrix<S> &m)
{

Integer i,j,k;

for (i=0; i<m.row; i++)
  {
          // If we have zero row
       if (m.pf[i]==m.pf[i+1]) 
       { for (k=0; k<i+1; k++) out << 0 << " ";
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

       for (k=m.pc[j-1]+1; k<i+1; k++) out << 0 << " ";

   out << std::endl;

  }
  return out;
}
 
template std::ostream & operator<<<Integer> (std::ostream & , const SymSparseMatrix<Integer>&);
template std::ostream & operator<<<Double> (std::ostream & , const SymSparseMatrix<Double>& );


} // end of namespace


