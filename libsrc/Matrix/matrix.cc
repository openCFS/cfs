#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <iomanip>

#include "matrix.hh"

namespace CoupledField
{      

template<class TYPE>
Matrix<TYPE>::Matrix ()
{	row = 0;
	col = 0;
	p = NULL;
}



template<class TYPE>
Matrix<TYPE>::Matrix (const Integer i, const Integer j)
{
    	if (i <= 0 || j <= 0) Error("invalid dimension",__FILE__,__LINE__);

	row = i;
	col = j;
	p = new TYPE* [row];

       if (!p) Error("out of memory",__FILE__,__LINE__); 

       p[0]=new TYPE[col*row];

      if (!p[0]) Error("out of memory",__FILE__,__LINE__); 

      for (Integer k=1; k < row; k++) p[k]=p[k-1]+col;
}

template<class TYPE>
Matrix<TYPE>::Matrix (const Integer i,const Vector<TYPE> * const x)
{	
    if (i <= 0) Error("invalid dimension",__FILE__,__LINE__);

	row = i;
	col = x [0].size();

	if (!col) Error("invalid dimension",__FILE__,__LINE__);

	p = new TYPE *[row];
	if (!p) Error("out of memory",__FILE__,__LINE__); 

        Integer k,kk;

        for (k=1; k < row; k++)

        { if (x[k].size()!=col)  Error(" Not all vectors for initialization has the same size",__FILE__,__LINE__);
        }

	for (k=0; k < row; k++)
         for (kk=0; kk<col; kk++)
               p[k][kk]=x[k][kk];
}

template<class TYPE>
Matrix<TYPE>::Matrix (const Matrix<TYPE> &x)
{	
  if (!x.row || !x.col)  Error("undefined Matrix",__FILE__,__LINE__);

	row = x.row;
	col = x.col;
        
//         if (p) 
// 	  { 
// 	    	    delete p[0];
// 	    	    delete [] p; 
// 	  }

	p = new TYPE * [row];
	if (!p) Error("out of memory", __FILE__, __LINE__);

        p[0]=new TYPE[row*col];
        if (!p[0]) Error("out of memory", __FILE__, __LINE__);

        Integer k;

        for (k=0; k < row*col; k++)  p[0][k]=x.p[0][k];
        for (k=1; k < row; k++) p[k]=p[k-1]+col;        
}

template<class TYPE>
void Matrix<TYPE> :: Resize(const Integer irow, const Integer icol)
{
   
   Integer k;

   if (p) {delete [] p[0];
           delete [] p;}

   row=irow; col=icol;

   p = new TYPE* [row];
   p[0]=new TYPE[row*col];
   if (!p[0] || !p) Error("out of memory", __FILE__, __LINE__);

  for (k=1; k < row; k++) p[k]=p[k-1]+col;
 
  for ( k = 0; k < row*col; k++) p[0][k]=0;
}


template<class TYPE>
void Matrix<TYPE> :: Resize(const Integer col)
{
  Resize(col,col);  
}


template<class TYPE>
TYPE * Matrix<TYPE>::operator[] (const Integer i) const
{ 
        if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
        if (i < 0 || i >= row) Error("invalid index",__FILE__,__LINE__);
 
        return p[i];
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator=(const Matrix<TYPE> &x)
{
       if (!x.row || !x.col) Error("undefined Matrix",__FILE__,__LINE__);
 
        if (this == &x)  return *this;
 
        Integer k;

        if (row!=x.row || col !=x.col )
        { 
           if (p)   {delete p[0];
                     delete [] p;}

           row=x.row; col=x.col; 

           p = new TYPE* [row];
           p[0]=new TYPE[col*row];

           for (k=1; k < row; k++) p[k]=p[k-1]+col;
           for ( k = 0; k < row*col; k++) p[0][k]=0;
 
           if (!p) Error("out of memory",__FILE__,__LINE__);
           if (!p[0]) Error("out of memory",__FILE__,__LINE__);
        }
         
        for (k = 0; k < row*col; k++)
           p[0][k]=x.p[0][k];

        return *this;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator+ () const
{
        if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
        return *this;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator+(const Matrix<TYPE> &x) const
{
        if (!row || !col || !x.row || !x.col)
            Error("undefined Matrix",__FILE__,__LINE__);
 
        if (row != x.row || col != x.col)
             Error("incompatible dimension",__FILE__,__LINE__);

        Matrix<TYPE> z(row,col);

         Integer k;
         for ( k = 0; k < row*col; k++)
                z [0][k] = x.p [0][k]+p[0][k];

        return z;
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator+=(const Matrix<TYPE> &x)
{
        if (!row || !col || !x.row || !x.col)
            Error("Don't use += to undefuned matrix",__FILE__,__LINE__);
 
         if (row != x.row || col != x.col)
             Error("incompatible dimension for +",__FILE__,__LINE__); 
 
        Integer k;
        for ( k = 0; k < row*col; k++)
                p [0][k] += x.p [0][k];
 
        return *this;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator-() const
{
        if (!row || !col )
            Error("undefined Matrix",__FILE__,__LINE__);
 
        Matrix<TYPE> z(row,col);
 
         Integer k;
         for ( k = 0; k < row*col; k++)
              if (p[0][k]!=0)  z [0][k] = -p [0][k];
 
        return z;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator-(const Matrix<TYPE> &x) const
{
        if (!row || !col || !x.row || !x.col)
            Error("undefined Matrix",__FILE__,__LINE__);
 
        if (row != x.row || col != x.col)
             Error("incompatible dimension",__FILE__,__LINE__);
 
        Matrix<TYPE> z(row,col);
 
         Integer k;
         for ( k = 0; k < row*col; k++)
                z [0][k] = -x.p [0][k]+p[0][k];
 
        return z;
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator-=(const Matrix<TYPE> &x)
{
        if (!row || !col || !x.row || !x.col)
            Error("Don't use += to undefuned matrix",__FILE__,__LINE__);
 
         if (row != x.row || col != x.col)
             Error("incompatible dimension for +",__FILE__,__LINE__);
 
        Integer k;
        for ( k = 0; k < row*col; k++)
                p [0][k] -= x.p [0][k];
 
        return *this;
}

template<class TYPE>
Integer Matrix<TYPE>::operator== (const Matrix<TYPE> &x) const
{
         if (!col || !row || !x.col || !x.row)
              Error("undefuned matrix",__FILE__,__LINE__);

         Integer k;

         for (k = 0; k < row*col; k++)
              if (p [0][k] != x.p[0][k]) return 0;
 
        return 1;
}

template<class TYPE>
Integer Matrix<TYPE>::operator!= (const Matrix<TYPE> &x) const
{
         if (!col || !row || !x.col || !x.row)
              Error("undefuned matrix",__FILE__,__LINE__); 

         Integer k;
         for (k = 0; k < row*col; k++)
              if (p [0][k] != x.p[0][k]) return 1;
 
        return 0;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator* (const TYPE &x) const 
{       if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
 
        Integer k;

        Matrix<TYPE> z(row,col);

        for ( k = 0; k < row*col; k++)
                z [0][k] = p[0][k]*x;
 
        return z;
}

template<class TYPE>
Vector<TYPE> Matrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
        if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
 
        if (!x.n) Error("undefined Vector",__FILE__,__LINE__);
 
        if (col != x.n) Error("incompatible dimension",__FILE__,__LINE__);
 
        Vector<TYPE> z(row);
 
        Integer k,kk;
        for ( k = 0; k < row; k++)
           for ( kk = 0; kk < col; kk++)
                z.p[k] += p[k][kk]*x[kk];
 
        return z;
}


template<class TYPE>
std::vector<TYPE> Matrix<TYPE>::operator*(const std::vector<TYPE> &x) const
{
  if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
 
  if (!x.size()) Error("undefined Vector",__FILE__,__LINE__);
 
  if (col != x.size()) Error("incompatible dimension",__FILE__,__LINE__);
 
  std::vector<TYPE> z(row);
 
  Integer k,kk;
  for ( k = 0; k < row; k++)
    for ( kk = 0; kk < col; kk++)
      z[k] += p[k][kk]*x[kk];
  
  return z;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator*(const Matrix<TYPE> &x) const
{
         if (!col || !row || !x.col || !x.row)
         Error("undefined Matrix",__FILE__,__LINE__);
 
        if (col != x.row)
        Error("incompatible dimension",__FILE__,__LINE__);
 
        TYPE    a;
        Matrix  z (row, x.col);

        Integer i,j; 
        for (i = 0; i < row; i++)
                for (j = 0; j < x.col; j++)
                {       a = p [i] [0] * x.p [0] [j];
                        for (Integer k = 1; k < col; k++)
                                a += p [i] [k] * x.p [k] [j];
                        z.p [i] [j] = a;
                }
 
        return z;
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator*= (const TYPE &x)
{
         if (!col || !row ) Error("undefined matrix",__FILE__,__LINE__);

        TYPE y=x;
 
        Integer i;
        for (i = 0; i < row*col; i++)
                p [0][i] *= y;
 
        return *this;
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator*=(const Matrix<TYPE> &x)
{       *this = *this * x;
 
        return *this;
}

template<class TYPE>
void Matrix<TYPE>::add_row (const Vector<TYPE> &x, const Integer pos)
{
      if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__); 

      TYPE ** help=new TYPE*[row+1];
      help[0]=new TYPE[col*(row+1)];
      Integer k;
      for (k=1; k < row+1; k++) help[k]=help[k-1]+col; 

      Integer i,ii;
      for (i=0; i < col; i++)
    {
      for (ii=0; ii < pos; ii++) help[ii][i]=p[ii][i];
                                 help[pos][i]=x[i];
      for (ii=pos+1; ii < row+1; ii++) help[ii][i]=p[ii-1][i];
    }
      row++;

      (*this).Resize(row,col);
 
      p=help;
      p[0]=help[0];
}

template<class TYPE>
void Matrix<TYPE>::add_col (const Vector<TYPE> &x, const Integer pos)
{
      if (!row || !col) Error("undefined Matrix",__FILE__,__LINE__);
 
      TYPE ** help=new TYPE*[row];
      help[0]=new TYPE[(col+1)*row];
      Integer k;
      for (k=1; k < row; k++) help[k]=help[k-1]+col+1;
 
      Integer i,ii;
      for (i=0; i < row; i++)
    {
      for (ii=0; ii < pos; ii++) help[i][ii]=p[i][ii];
                                 help[i][pos]=x[i];
      for (ii=pos+1; ii < col+1; ii++) help[i][ii]=p[i][ii-1];
    }

      col++;
 
      (*this).Resize(row,col);
 
      p=help;
      p[0]=help[0]; 
}

template<class TYPE>
void Matrix<TYPE>::cut(const Integer ai, const Integer aj)
{
if (ai < 0 || ai > row-1 || aj < 0 || aj > col-1)
         Error("invalid index in function cut",__FILE__,__LINE__);
if (!row) Error("Matrix is undefined in function cut",__FILE__,__LINE__);

Integer n=(row-1)*(col-1);
Matrix<TYPE> help(row-1,col-1);

Integer i, j, k=0;
for (i=0; i<row ; i++)
if (i!=ai)  for (j=0; j<col; j++)
if (j!=aj)
   {
     help.p[0][k]=p[i][j];
     k++;
   }

delete  p[0];
delete [] p;
row--; col--;
 
p = new TYPE* [row];
p[0]=new TYPE[n];
for (k=1; k < row; k++) p[k]=p[k-1]+col;

for (i=0; i<n; i++) p[0][i]=help.p[0][i];
}
 

// does not work in this way, because z is allocated INSIDE the method, and, thererfore, "dies" 
// a horrible death outside !!!!!!
//
// template<class TYPE>
// Matrix<TYPE> Matrix<TYPE>::part (const Integer i1,const Integer i2, const Integer j1,
//                                  const Integer j2) const
// {
//          if (!col || !row ) Error("undefined Matrix",__FILE__,__LINE__); 
 
//         if (i1 < 0 || i1 > i2 || i2 >= row)
//                       Error("invalid index",__FILE__,__LINE__); 
//         if (j1 < 0 || j1 > j2 || j2 >= col)
//                     Error("invalid index",__FILE__,__LINE__);   
//         Matrix  z (i2 - i1 + 1, j2 - j1 + 1);
 
//         for (Integer i = 0; i < z.row; i++)
//                 for (Integer j = 0; j < z.col; j++)
//                         z [i] [j] = p [i1 + i] [j1 + j];
 
//         return z;
// }

template< class TYPE>
void Matrix<TYPE>::precond(Vector<TYPE> & e, const Vector<TYPE> r, enum precond type)
{

#ifdef TRACE
  (*trace) << "entering Matrix::precond" << std::endl;
#endif
 
 if (!e.size()) Error("Define vector before use of precond");
 
 Integer i,j;
 Integer n=r.size();
 Double omega;

 switch(type)
{ case non: e=r; break;
  case Jacobi: 
       for (i=0; i<n; i++) e[i]=r[i]/p[i][i];  ///// May be (*this)
       break;
  case SSOR:
       omega=1.2; /// omega=1.2
       TYPE sum;

       for (i=0; i<n; i++)
       {
         sum=0;
         for (j=0; j<i; j++) sum+=p[i][j]*e[j];
         e[i]=(TYPE)((r[i]-omega*sum)/p[i][i]);
       }

       for (i=n-1; i>=0; i--)
       {
         sum=0;
         for (j=i+1; j<n; j++) sum+=p[i][j]*e[j];
         e[i]-= (TYPE)(omega*sum/p[i][i]);
       }
       
       e*= (TYPE)(omega*(2-omega));
       break;
default:
       Error("Wrong type of precondition",__FILE__,__LINE__);
              
} //end of switch
}

template<class S>
std::ostream & operator << (std::ostream & out, const Matrix<S> &mat)
{
out.setf(std::ios::scientific);

for (Integer i=0; i < mat.size_row(); i++)
{for (Integer j=0; j < mat.size_col(); j++)
 out << mat[i][j] << " ";
 out << std::endl;}

//out.setf(0, std::ios::floatfield);

 return out;
}

template std::ostream & operator<<<Integer> (std::ostream & , const Matrix<Integer> &);
template std::ostream & operator<<<Double> (std::ostream & , const Matrix<Double> &);

template<class TYPE>
TYPE Spur (const Matrix<TYPE> &x)
{
if (!x.size_row() || !x.size_col()) Error("undefined Matrix",__FILE__,__LINE__);
 
if (x.size_row() != x.size_col())
  Error("incompatible dimension",__FILE__,__LINE__);

        TYPE    a;
 
        a = x.get() [0] [0];
 
        Integer i;
        for ( i = 1; i < x.size_row(); i++)
                a += x.get() [i] [i];
 
        return a;
}


template<class TYPE>
void Matrix<TYPE>::Transpose (Matrix<TYPE> &transposedMat)
{
  transposedMat.Resize(col, row);
  
  Integer i,j;
   
  for (i = 0; i < col; i++)
    for (j = 0; j < row; j++)
      transposedMat[i] [j] = p[j] [i];
}


template<class TYPE>
void Matrix<TYPE>::DyadicMult(std::vector<TYPE> v1)
{
  DyadicMult(v1, v1);
}




template<class TYPE>
void Matrix<TYPE>::DyadicMult(std::vector<TYPE> v1, std::vector<TYPE> v2)
{
  Integer row = v1.size();
  Integer col = v2.size();
  
  this->Resize(row,col);

  for(Integer actRow=0; actRow<row; actRow++)
    for(Integer actCol=0; actCol < col; actCol++)
      p[actRow][actCol] = v1[actRow] * v2[actCol];
}

template<class TYPE>
void Matrix<TYPE>::Invert (Matrix <Double> & inv) const
{       
  if (!row || !col) Error("Undefined Matrix!",__FILE__,__LINE__);
  if (row != col ) Error("No quadratic matrix!",__FILE__,__LINE__);


  switch (row)
    {
    case 1: 
      inv.Resize(1);
      inv[0][0] = 1/p[0][0];
      break;
    case 2:
      inv.Resize(2,2);
      inv[0][0] = p[1][1];
      inv[0][1] = - p[0][1];
      inv[1][0] = - p[1][0];
      inv[1][1] = p[0][0];
      inv *= 1/this->Det();
      break;

    case 3:
      // see Stöcker: "Taschenbuch Mathematischer Formeln und Moderner Verfahren" p.418
      inv.Resize(3,3);
      for(Integer i=0; i<3; i++)
	for(Integer j=0; j<3; j++)
	  inv[j][i] = Adjunct(i,j);      
      
      inv *= 1/this->Det();      
      break;
      
    default: 
      Error("Inversion not implemented fo dimension larger than 2!",__FILE__,__LINE__);
    }
}


template<class TYPE>
 Double Matrix<TYPE>::Adjunct (Integer i, Integer j) const
{       
  if (!row || !col) Error("Undefined Matrix!",__FILE__,__LINE__);
  if (row != col ) Error("No quadratic matrix!",__FILE__,__LINE__);
  if (row != 3 ) Error("Matrix::Adjunct only implemented for matrix size 3!",__FILE__,__LINE__);  
  
  std::vector<Integer> iVec(2);
  std::vector<Integer> jVec(2);
  Integer runningIndexI = 0;
  Integer runningIndexJ = 0;
  
  for (Integer actI = 0; actI<=2; actI++)
    {
      if (actI != i)
	{	  
	  iVec[runningIndexI] = actI;
	  runningIndexI++;
	}
      
      if (actI != j)
	{
	  jVec[runningIndexJ] = actI;
	  runningIndexJ++;
	}
    }

  Double adj = pow(-1,i+j) * 
    ( p[iVec[0]][jVec[0]] * p[iVec[1]][jVec[1]] - 
      p[iVec[0]][jVec[1]] * p[iVec[1]][jVec[0]]);


  return adj;
  
}

template Integer Spur<Integer>(const Matrix<Integer> &);
template Double Spur<Double>(const Matrix<Double> &);

// template<class TYPE>
// Matrix<TYPE> Trans (const Matrix<TYPE> &x)
// {       if (!x.size_row() || !x.size_col())
//             Error("undefined Matrix",__FILE__,__LINE__);
 
//         Matrix<TYPE>    z (x.size_col(), x.size_row());
 
//         Integer i,j;
//         for (i = 0; i < x.size_col(); i++)
//                 for (j = 0; j < x.size_row(); j++)
//                         z[i] [j] = x [j] [i];
//         return z;
// }
 
// template Matrix<Integer> Trans<Integer>(const Matrix<Integer> &);
// template Matrix<Double> Trans<Double>(const Matrix<Double> &);

} // end of namespace
