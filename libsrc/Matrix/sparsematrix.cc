#include <iostream>
#include <fstream>
#include <string>

#include "sparsematrix.hh"

namespace CoupledField
{

template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix()
{
 numrows=0;
 numcol=0;
 ptRow=0;
}

template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix(const Integer anumrows, const Integer anumcol)
{
 if (anumrows <0 || anumcol <0) Error("Invalid index in constructor of SparseMatrix", __FILE__,__LINE__);

 numrows=anumrows;
 numcol=anumcol;
 ptRow=new ElemSparseMatrix<TYPE> * [numrows];

 if (!ptRow)
  Error("Not enought memory for constructor of SparseMatrix",__FILE__,__LINE__);
 
 Integer i;
 for (i=0; i<numrows; i++) ptRow[i]=0;
}

template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix(const SparseMatrix<TYPE> & x)
{
 if (!x.numrows)  Error("undefined SparseMatrix in Copy-Initializer Constructor",__FILE__,__LINE__);

 if (!ptRow) delete [] ptRow;

 numrows=x.numrows;
 numcol=x.numcol;

  ptRow=new ElemSparseMatrix<TYPE> * [numrows];

 if (!ptRow)
 Error("Not enought memory for constructor of SparseMatrix",__FILE__,__LINE__);

 Integer i;
 for (i=0; i<numrows; i++) ptRow[i]=0;

 Copy(x);
}

template<class TYPE>
void SparseMatrix<TYPE>::Copy(const SparseMatrix<TYPE> & x)
{
Integer irow;
for (irow=0; irow<numrows; irow++)
 {
   ElemSparseMatrix<TYPE> * ptElem, * ptElemX;
   Integer first=1;
   for (ptElemX=x.ptRow[irow]; ptElemX!=0; ptElemX=ptElemX->next)
    {
     ElemSparseMatrix<TYPE> * newElem=new ElemSparseMatrix<TYPE>;

     if (!newElem) Error("Not enought memory", __FILE__,__LINE__);
  
     newElem->column=ptElemX->column;
     newElem->value=ptElemX->value;
     newElem->next=0;
     if (first==1)
         { ptElem=ptRow[irow]=newElem;
           first=0;
         }
     else
        {
          ptElem->next=newElem;
          ptElem=ptElem->next;
        }
    }
 }
}

//! Constructor(const Matrix<TYPE> & x)
template<class TYPE>
SparseMatrix<TYPE>::SparseMatrix(const Matrix<TYPE> & x)
{
 Error("Not implemented yet");
}

template<class TYPE>
void SparseMatrix<TYPE>::Resize(const Integer anumrows, const Integer anumcol)
{
 if (!ptRow) delete [] ptRow;

 numrows=anumrows;
 numcol=anumcol;

 ptRow=new ElemSparseMatrix<TYPE> * [numrows];

 if (!ptRow)
 Error("Not enought memory for constructor of SparseMatrix",__FILE__,__LINE__);

 Integer i;
 for (i=0; i<numrows; i++) ptRow[i]=0;
 
}

template<class TYPE>
TYPE & SparseMatrix<TYPE>::operator()(Integer row,Integer col)
{
if (row<0 || row >=numrows || col <0 || col>=numcol) 
  Error("Invalid index in SparsMatrix", __FILE__, __LINE__);

ElemSparseMatrix<TYPE> * ptElem=ptRow[row];
if (ptElem==0 || col<ptElem->column)
   return InsertElement(ptElem,row,col,1);
if (ptElem!=0 && col==ptElem->column)
   return ptElem->value;
for (;ptElem->next !=0; ptElem=ptElem->next)
{
if (col==ptElem->next->column)
   return ptElem->next->value;
else if (col<ptElem->next->column) break;
}
return InsertElement(ptElem,row,col,0);
}

template<class TYPE>
TYPE & SparseMatrix<TYPE>::InsertElement(ElemSparseMatrix<TYPE> * ptElem,const Integer row, const Integer col, const Integer first)
{
ElemSparseMatrix<TYPE> * newElem=new ElemSparseMatrix<TYPE>;

if (!newElem) Error("Not enought memory");

newElem->column=col;
newElem->value=0.;

if (first==1)
{
 newElem->next=ptElem;
 ptRow[row]=newElem; 
}
else
{
newElem->next=ptElem->next;
ptElem->next=newElem;
}

return newElem->value;
}

template<class TYPE>
Vector<TYPE>    SparseMatrix<TYPE>::operator*       (const Vector<TYPE> & vec) const
{
 if (!numrows) Error("undefined SparseMatrix",__FILE__,__LINE__);

 if (numcol!=vec.size()) Error("incompatible dimensions",__FILE__,__LINE__);

 Vector<TYPE> result;

}

template<class TYPE>
void SparseMatrix<TYPE>::DeleteElement(const Integer irow, const Integer icol)
{
if (ptRow==0) return;
if (irow < 0 || irow >= numrows) return;
if (icol < 0 || icol >= numcol) return;

ElemSparseMatrix<TYPE> * ptElem,*tmp;

ptElem=ptRow[irow];
if (ptElem==0) return;
if (icol == ptElem->column)
 {
 ptRow[irow]=ptElem->next;
 delete ptElem;
 return;
 }
for(;ptElem->next !=0; ptElem=ptElem->next)
 {
if (icol==ptElem->next->column)
  {
 tmp=ptElem->next;
 ptElem->next=ptElem->next->next;
 delete tmp;
 return;
  }
 }
}

template<class TYPE>
SparseMatrix<TYPE> & SparseMatrix<TYPE>::operator=(const SparseMatrix<TYPE> & x)
{

if (!x.numrows) Error("undefined SparseMatrix",__FILE__,__LINE__);

if (this==&x) return *this;

if (x.numrows!=numrows)
{
if (!ptRow) delete [] ptRow;

 numrows=x.numrows;

  ptRow=new ElemSparseMatrix<TYPE> * [numrows];

 if (!ptRow)
 Error("Not enought memory for constructor of SparseMatrix",__FILE__,__LINE__);
}

 numcol=x.numcol;

 Integer i;
 for (i=0; i<numrows; i++) ptRow[i]=0;

 Copy(x);
}

template<class TYPE>
SparseMatrix<TYPE>  &  SparseMatrix<TYPE>::operator+ () 
{
 if (!numrows)  Error("undefined Matrix",__FILE__,__LINE__);
 return *this;
}

template<class TYPE>
SparseMatrix<TYPE>  & SparseMatrix<TYPE>::operator+= (const SparseMatrix<TYPE> &x) 
{
 if (!numrows || !x.numrows)
   Error("Don't use += to undefined matrix",__FILE__,__LINE__);

  if (numrows != x.numrows || numcol != x.numcol)
             Error("incompatible dimension for +",__FILE__,__LINE__);

  Integer irow;
  for (irow=0; irow <numrows; irow++)
  {
   ElemSparseMatrix<TYPE> *ptElemX,*ptElem;

   ptElem=ptRow[irow];
   ptElemX=x.ptRow[irow];

   Integer first=1;
  
   std::cout << " row " << irow << std::endl;
 
   if (ptElemX!=0)
   {
     if (ptElem==0)
        for (; ptElemX!=0; ptElemX=ptElemX->next)
          (*this)(irow,ptElemX->column)=ptElemX->value;
     else
      for (; ptElemX!=0; ptElemX=ptElemX->next)
       {
       std::cout << ptElem->value <<"(" << ptElem->column << ")" << " "; 
       (*this)(irow,ptElemX->column)+=ptElemX->value;
     //   std::cout << "we are here" << std::endl;
     //  	while ((ptElem->column < ptElemX->column) && ptElem->next!=0)
     //            { ptElem=ptElem->next; first=0;}
     //   if (ptElem->column == ptElemX->column) 
     //          ptElem->value+= ptElemX->value;
     //   if (ptElem->column > ptElemX->column)
     //        InsertElement(ptElemX,irow,ptElemX->column,first);
       }
   }
  }
}     

template<class TYPE>
SparseMatrix<TYPE>  SparseMatrix<TYPE>::operator+ (const SparseMatrix<TYPE> &x) const
{
 if (!numrows || !x.numrows)
   Error("Don't use += to undefined matrix",__FILE__,__LINE__);

  if (numrows != x.numrows || numcol != x.numcol)
             Error("incompatible dimension for +",__FILE__,__LINE__);

  SparseMatrix<TYPE> result;

  result=(*this);
  result+=x;

  return result;
}

template<class TYPE>
std::ostream & operator << (std::ostream & out, const SparseMatrix<TYPE> &mat)
{
out.setf(std::ios::scientific);

Integer i,ii;
ElemSparseMatrix<TYPE> * ptElem;

for (i=0; i < mat.numrows; i++)
{
 ptElem=mat.ptRow[i];
 if (ptElem == 0) for (ii=0; ii<mat.numcol; ii++) out << 0.0 << " ";
 else
    { for (ii=0; ii<ptElem->column; ii++)  out << 0.0 << " ";
      out << ptElem->value <<" ";
      for (; ptElem->next!=0; ptElem=ptElem->next)
       {
        for (ii=ptElem->column+1; ii<ptElem->next->column; ii++)
          out << 0.0 << " ";
        out << ptElem->next->value<<" ";
       }
      for (ii=ptElem->column+1; ii <mat.numcol; ii++)
        out << 0.0 << " ";
    }
  out << std::endl;
}

out.setf(0, std::ios::floatfield);

return out;

}

template<class TYPE>
void SparseMatrix<TYPE>::test()
{
 Integer i;
 for (i=0; i< numrows; i++)
 {
   ElemSparseMatrix<TYPE> * ptElem=ptRow[i];
   if (ptElem==0) std::cout << "zero" << std::endl;
   else
 {  for (; ptElem!=0; ptElem=ptElem->next)
 { std::cout << ptElem->value <<"(" << ptElem->column << ")" << " "; }
  std::cout << std::endl;
 }

 }
}

template std::ostream & operator<<<Integer> (std::ostream & , const SparseMatrix<Integer> &);
template std::ostream & operator<<<Double> (std::ostream & , const SparseMatrix<Double> &);

} // end of namespace
