#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <iomanip>

#include "matrix.hh"
#include <Utils/vector.hh>

namespace CoupledField
{      

template<class TYPE>
Matrix<TYPE>::Matrix ()
{
  ENTER_FCN("Matrix::Matrix");	
  size_row_ = 0;
  size_col_ = 0;
  data_ = NULL;
}


template<class TYPE>
Matrix<TYPE>::Matrix (const Integer nRows, const Integer nCols)
{
  ENTER_FCN("Matrix::Matrix");
#ifdef CHECK_INDEX 
  if (nRows <= 0 || nCols <= 0) Error("invalid dimension",__FILE__,__LINE__);
#endif

  size_row_ = nRows;
  size_col_ = nCols;
  data_ = new TYPE* [size_row_];

  data_[0]=new TYPE[size_col_*size_row_];

  for (Integer k=1; k < size_row_; k++) 
    data_[k]=data_[k-1]+size_col_;
}


template<class TYPE>
Matrix<TYPE>::Matrix (const Integer nRows,const Vector<TYPE> * const x)
{
  ENTER_FCN("Matrix::Matrix");
#ifdef CHECK_INDEX
  if (nRows <= 0) Error("invalid dimension",__FILE__,__LINE__);
#endif
  
  size_row_ = nRows;
  size_col_ = x[0].size_;

#ifdef CHECK_INDEX 
  if (size_col_ == 0) Error("invalid dimension",__FILE__,__LINE__);
#endif 

  data_ = new TYPE *[size_row_];

  Integer k,kk;
#ifdef CHECK_INDEX
  for (k=1; k < size_row_; k++)
    
    { if (x[k].size_!=size_col_)  Error(" Not all vectors for initialization have the same size",__FILE__,__LINE__);
    }
#endif
  
  for (k=0; k < size_row_; k++)
    for (kk=0; kk<size_col_; kk++)
      data_[k][kk]=x[k][kk];
}

template<class TYPE>
Matrix<TYPE>::Matrix (const Matrix<TYPE> &x)
{
  ENTER_FCN("Matrix::Matrix");

#ifdef CHECK_INITIALIZED
 if (x.size_row_ == 0 || x.size_col_ == 0)  Error("undefined Matrix",__FILE__,__LINE__);
#endif

 
 size_row_ = x.size_row_;
 size_col_ = x.size_col_;
 
 data_ = new TYPE * [size_row_];
 data_[0]=new TYPE[size_row_ * size_col_];
 
 
 Integer k;
 
 for (k=0; k < size_row_*size_col_; k++)  
   data_[0][k]=x.data_[0][k];
 for (k=1; k < size_row_; k++) 
   data_[k]=data_[k-1]+size_col_;        
}

template<class TYPE>
Matrix<TYPE>::~Matrix ()
{
  ENTER_FCN("Matrix::~Matrix");
  if (data_ != NULL)
    {
      delete[] data_[0];
      delete[] data_;
    }
}

template<class TYPE>
Boolean Matrix<TYPE>::IsComplex()
{
  return FALSE;
}

template<>
Boolean Matrix<Complex>::IsComplex()
{
  return TRUE;
}


template<class TYPE>
void Matrix<TYPE> :: Resize(const Integer nRows, const Integer nCols)
{
  ENTER_FCN("Matrix::Resize");
  
  Integer k;
  
  if (nRows != size_row_ || nCols != size_col_)
    {
      
      if (data_ != NULL) {
	delete [] data_[0];
	delete [] data_;
      }
      
      size_row_ = nRows; 
      size_col_ = nCols;
      
      data_ = new TYPE* [size_row_];
      data_[0]=new TYPE[size_row_*size_col_];
      
      for (k=1; k < size_row_; k++) 
	data_[k]=data_[k-1]+size_col_;
    }
  
  // initialize values to 0
  for ( k = 0; k < size_row_ * size_col_; k++) 
    data_[0][k]=0;
}


template<class TYPE>
void Matrix<TYPE>::Resize(const Integer col)
{
  ENTER_FCN("Matrix::Resize");
  Resize(col,col);  
}




template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator=(const Matrix<TYPE> &x)
{
  ENTER_IFCN("Matrix::operator=");

#ifdef CHECK_INITIALIZED
  if (x.size_row_ == 0 || x.size_col_ == 0) Error("undefined Matrix",__FILE__,__LINE__);
#endif  

  if (this == &x)  return *this;
  
  Integer k;
  
  if (size_row_ != x.size_row_ || size_col_ != x.size_col_ )
    {
      
      if (data_)
	{
	  delete[] data_[0];
	  delete[] data_;
	}
      
      size_row_ = x.size_row_; 
      size_col_ = x.size_col_; 
      
      data_ = new TYPE* [size_row_];
      data_[0]=new TYPE[size_row_*size_col_];
  
      for (k=1; k < size_row_; k++) 
	data_[k]=data_[k-1]+size_col_;
    }
  
  for ( k = 0; k < size_row_ * size_col_; k++) 
    data_[0][k]=x.data_[0][k];
  
  return *this;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator+ () const
{
  ENTER_IFCN("Matrix::operator+");
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) Error("undefined Matrix",__FILE__,__LINE__);
#endif

  return *this;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator+(const Matrix<TYPE> &x) const
{
  ENTER_IFCN("Matrix::operator+");
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) Error("undefined Matrix",__FILE__,__LINE__);
#endif 

#ifdef CHECK_INDEX
  if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
    Error("incompatible dimension",__FILE__,__LINE__);
#endif
  
  Matrix<TYPE> z(size_row_,size_col_);
  
  Integer k;
  for ( k = 0; k < size_row_*size_col_; k++)
    z [0][k] = x.data_ [0][k]+data_[0][k];
  
  return z;
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator+=(const Matrix<TYPE> &x)
{
  ENTER_IFCN("Matrix::operator+=");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 || 
      x.size_row_ == 0 || x.size_col_ == 0)
    Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX  
  if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
    Error("incompatible dimension for +",__FILE__,__LINE__); 
#endif
  Integer k;
  
  for ( k = 0; k < size_row_ * size_col_; k++)
    data_[0][k] += x.data_[0][k];
  
  return *this;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator-() const
{
  ENTER_IFCN("Matrix::operator-");
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 )
    Error("undefined Matrix",__FILE__,__LINE__);
#endif

  Matrix<TYPE> z(size_row_,size_col_);
  
  Integer k;
  for ( k = 0; k < size_row_*size_col_; k++)
    z [0][k] = -data_[0][k];
  
  return z;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator-(const Matrix<TYPE> &x) const
{
  ENTER_IFCN("Matrix::operator-");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 || 
      x.size_row_ == 0 || x.size_col_ == 0)
    Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX  
  if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
    Error("incompatible dimension for +",__FILE__,__LINE__); 
#endif
  
  Matrix<TYPE> z(size_row_,size_col_);
  
  Integer k;
  for ( k = 0; k < size_row_*size_col_; k++)
    z[0][k] = -x.data_[0][k]+data_[0][k];
  
  return z;
}

template<class TYPE>
Matrix<TYPE> & Matrix<TYPE>::operator-=(const Matrix<TYPE> &x)
{
  ENTER_IFCN("Matrix::operator-=");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 || 
      x.size_row_ == 0 || x.size_col_ == 0)
    Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX  
  if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
    Error("incompatible dimension for +",__FILE__,__LINE__); 
#endif
  Integer k;
  for ( k = 0; k < size_row_ * size_col_; k++)
    data_ [0][k] -= x.data_ [0][k];
  
  return *this;
}

template<class TYPE>
Integer Matrix<TYPE>::operator== (const Matrix<TYPE> &x) const
{
  ENTER_IFCN("Matrix::operator==");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 || 
      x.size_row_ == 0 || x.size_col_ == 0)
    Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
  Integer k;
  
  for (k = 0; k < size_row_*size_col_; k++)
    if (data_ [0][k] != x.data_[0][k]) return FALSE;
  
  return TRUE;
}

template<class TYPE>
Integer Matrix<TYPE>::operator!= (const Matrix<TYPE> &x) const
{
  ENTER_IFCN("Matrix::operator!=");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 || 
      x.size_row_ == 0 || x.size_col_ == 0)
    Error("undefined Matrix",__FILE__,__LINE__);
#endif 
  
  Integer k;
  for (k = 0; k < size_row_*size_col_; k++)
    if (data_ [0][k] != x.data_[0][k]) return FALSE;
  
  return TRUE;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator* (const TYPE &x) const 
{ 
  ENTER_IFCN("Matrix::operator*");
  
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
  Integer k;
  
  Matrix<TYPE> z(size_row_,size_col_);
  
  for ( k = 0; k < size_row_*size_col_; k++)
    z [0][k] = data_[0][k]*x;
  
  return z;
}

template<class TYPE>
Vector<TYPE> Matrix<TYPE>::operator*(const Vector<TYPE> &x) const
{
  ENTER_IFCN("Matrix::operator*");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
  if (x.size_ == 0) Error("undefined Vector",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (size_col_ != x.size_) Error("incompatible dimension",__FILE__,__LINE__);
#endif
  
  Vector<TYPE> z(size_row_);
  
  Integer k,kk;
  for ( k = 0; k < size_row_; k++)
    for ( kk = 0; kk < size_col_; kk++)
      z.data_[k] += data_[k][kk]*x.data_[kk];
  
  return z;
}

template<class TYPE>
Matrix<TYPE> Matrix<TYPE>::operator*(const Matrix<TYPE> &x) const
{
  ENTER_IFCN("Matrix::operator*");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0 || 
      x.size_row_ == 0 || x.size_col_ == 0)
    Error("undefined Matrix",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX  
  if (size_col_ != x.size_row_)
    Error("incompatible dimension",__FILE__,__LINE__);
#endif
 
  TYPE    a;
  Matrix  z (size_row_, x.size_col_);
  
  Integer i,j; 
  for (i = 0; i < size_row_; i++)
    for (j = 0; j < x.size_col_; j++)
      {       
	a = data_ [i] [0] * x.data_ [0] [j];
	for (Integer k = 1; k < size_col_; k++)
	  a += data_ [i] [k] * x.data_ [k] [j];
	z.data_ [i] [j] = a;
      }
  
  return z;
}


template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator*= (const TYPE &x)
{
  ENTER_IFCN("Matrix::operator*=");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
  TYPE y=x;
  
  Integer i;
  for (i = 0; i < size_row_*size_col_; i++)
    data_ [0][i] *= y;
  
  return *this;
}

  // Perform a matrix-vector multiplication rvec = this*mvec
template<class TYPE>
void Matrix<TYPE>::Mult(CFSVector & mvec, CFSVector & rvec)
{
  ENTER_IFCN("Matrix::Mult");
  Vector<TYPE> & mvec1 = dynamic_cast<Vector<TYPE>& >(mvec);
  Vector<TYPE> & rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
  Integer size_mvec = mvec1.GetSize();
  Integer size_rvec = rvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
  if (size_mvec == 0) 
    Error("undefined Vector",__FILE__,__LINE__);
  if (size_rvec == 0) 
    Error("undefined Vector",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (size_col_ != size_mvec) Error("incompatible dimension",__FILE__,__LINE__);
  if (size_row_ != size_rvec) Error("incompatible dimension",__FILE__,__LINE__);
#endif
   
  Integer k,kk;
  for ( k = 0; k < size_row_; k++)
    for ( kk = 0; kk < size_col_; kk++)
        rvec1[k] += data_[k][kk]*mvec1[kk];

}

// Perform a matrix-matrix multiplication rMat = this*mMat
template<class TYPE>
void Matrix<TYPE>::Mult(CFSMatrix & mMat, CFSMatrix & rMat)
{
  ENTER_IFCN("Matrix::Mult");
  Matrix<TYPE> & mMat1 = dynamic_cast<Matrix<TYPE> & >(mMat);
  Matrix<TYPE> & rMat1 = dynamic_cast<Matrix<TYPE>& >(rMat);
  
  Integer size_mMatRow = mMat1.GetSizeRow();
  Integer size_mMatCol = mMat1.GetSizeCol();

  Integer size_rMatRow = rMat1.GetSizeRow();
  Integer size_rMatCol = rMat1.GetSizeCol();
 
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
  if (size_mMatRow == 0 || size_mMatCol==0) 
    Error("undefined Matrix",__FILE__,__LINE__);
  if (size_rMatRow == 0||size_rMatCol==0) 
    Error("undefined Matrix",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (size_col_ != size_mMatRow) Error("incompatible dimension while matrix-matrix multiplication",__FILE__,__LINE__);
  if (size_row_ != size_rMatRow) Error("incompatible dimension while matrix-matrix multiplication",__FILE__,__LINE__);
  if (size_mMatCol != size_rMatCol) Error("incompatibel dimension while matrix-matrix multiplication",__FILE__,__LINE__);
#endif
   
  Vector<TYPE> temp(1);

  for (Integer i = 0; i < size_row_; i++)
    for (Integer j = 0; j < size_mMatCol; j++)
      {       
	temp = data_[i][0] * mMat1[0][j];
	for (Integer k = 1; k <size_mMatRow; k++)
	  temp[0] += data_[i][k] * mMat1[k][j];
	rMat1[i][j] = temp[0];
      }  
} 


// Perform a matrix(Double)-vector(Complex) multiplication rvec = this*mvec
// where the matrix is supposed to be of type Double, rvec and mvec are complex valued
template<>
void Matrix<Double>::MatVecMult_DC(Vector<Complex> & mvec, Vector<Complex> & rvec){
 ENTER_IFCN("Matrix::MatVecMult_DC");
  
  Integer size_mvec = mvec.GetSize();
  Integer size_rvec = rvec.GetSize();
 
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
  if (size_mvec == 0) 
    Error("undefined Vector",__FILE__,__LINE__);
  if (size_rvec == 0) 
    Error("undefined Vector",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (size_col_ != size_mvec) Error("incompatible dimension",__FILE__,__LINE__);
  if (size_row_ != size_rvec) Error("incompatible dimension",__FILE__,__LINE__);
#endif
   
  Integer k,kk;
  for ( k = 0; k < size_row_; k++)
    for ( kk = 0; kk < size_col_; kk++)
        rvec[k] += data_[k][kk]*mvec[kk];

}

// Perform a matrix(Double)-vector(Complex) multiplication rvec = this*mvec
// where the matrix is supposed to be of type Double, rvec and mvec are complex valued
template<>
void Matrix<Complex>::MatVecMult_CD(Vector<Double> & mvec, Vector<Complex> & rvec){
 ENTER_IFCN("Matrix::MatVecMult_CD");

  Integer size_mvec = mvec.GetSize();
  Integer size_rvec = rvec.GetSize();
 
#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
  if (size_mvec == 0) 
    Error("undefined Vector",__FILE__,__LINE__);
  if (size_rvec == 0) 
    Error("undefined Vector",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (size_col_ != size_mvec) Error("incompatible dimension",__FILE__,__LINE__);
  if (size_row_ != size_rvec) Error("incompatible dimension",__FILE__,__LINE__);
#endif
   
  Integer k,kk;
  for ( k = 0; k < size_row_; k++)
    for ( kk = 0; kk < size_col_; kk++)
        rvec[k] += data_[k][kk]*mvec[kk];

}


template<class TYPE>
Matrix<TYPE> & Matrix<TYPE>::operator*=(const Matrix<TYPE> &x)
{   
  ENTER_IFCN("Matrix::operator*=");    
  *this = *this * x;
  
  return *this;
}

template<class TYPE>
Matrix<TYPE> &Matrix<TYPE>::operator/= (const TYPE &x)
{
  ENTER_IFCN("Matrix::operator/=");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
#endif

  TYPE y=x;
  
  Integer i;
  for (i = 0; i < size_row_*size_col_; i++)
    data_ [0][i] /= y;
  
  return *this;
}


template<class TYPE>
void Matrix<TYPE>::AddRow (const Vector<TYPE> &x, const Integer pos)
{
  ENTER_IFCN("Matrix::add_row");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
#endif


  TYPE ** help=new TYPE*[size_row_+1];
  help[0]=new TYPE[size_col_ * (size_row_+1)];

  Integer k;
  for (k=1; k < size_row_+1; k++) 
    help[k]=help[k-1]+size_col_; 
  
  Integer i,ii;
  for (i=0; i < size_col_; i++)
    {
      for (ii=0; ii < pos; ii++) 
	help[ii][i]=data_[ii][i];
      help[pos][i]=x[i];
      for (ii=pos+1; ii < size_row_+1; ii++) 
	help[ii][i]=data_[ii-1][i];
    }
  size_row_++;
  
  // Inefficient?
  (*this).Resize(size_row_,size_col_);
  
  data_=help;
  data_[0]=help[0];
}


template<class TYPE>
void Matrix<TYPE>::AddColumn (const Vector<TYPE> &x, const Integer pos)
{
  ENTER_IFCN("Matrix::AddColumn");

#ifdef CHECK_INITIALIZED
  if (size_row_ == 0 || size_col_ == 0) 
    Error("undefined Matrix",__FILE__,__LINE__);
#endif

  TYPE ** help=new TYPE*[size_row_];
  help[0]=new TYPE[(size_col_+1)*size_row_];
  Integer k;
  for (k=1; k < size_row_; k++) 
    help[k]=help[k-1]+size_col_+1;
  
  Integer i,ii;
  for (i=0; i < size_row_; i++)
    {
      for (ii=0; ii < pos; ii++) 
	help[i][ii]=data_[i][ii];
      help[i][pos]=x[i];
      for (ii=pos+1; ii < size_col_+1; ii++) 
	help[i][ii]=data_[i][ii-1];
    }
  
  size_col_++;
  
  (*this).Resize(size_row_,size_col_);
  
  data_=help;
  data_[0]=help[0]; 
}

template<class S>
std::ostream & operator << (std::ostream & out, const Matrix<S> &mat)
{
  ENTER_IFCN("Matrix::operator<<");
  out.setf(std::ios::scientific);
  
  for (Integer i=0; i < mat.GetSizeRow(); i++)
    {for (Integer j=0; j < mat.GetSizeCol(); j++)
      out << mat[i][j] << " ";
    out << std::endl;}
  
  //out.setf(0, std::ios::floatfield);
  
  return out;
}

template std::ostream & operator<<<Integer> (std::ostream & , const Matrix<Integer> &);
template std::ostream & operator<<<Double> (std::ostream & , const Matrix<Double> &);


template<class TYPE>
void Matrix<TYPE>::Transpose (Matrix<TYPE> &transposedMat)
{
ENTER_FCN("Matrix::Transpose");
 transposedMat.Resize(size_col_, size_row_);
 
 Integer i,j;
 
 for (i = 0; i < size_col_; i++)
   for (j = 0; j < size_row_; j++)
     transposedMat.data_[i][j] = data_[j] [i];
}


template<class TYPE>
void Matrix<TYPE>::DyadicMult(CFSVector & v1)
{
  ENTER_FCN("Matrix::DyadicMult");
  DyadicMult(v1, v1);
}




template<class TYPE>
void Matrix<TYPE>::DyadicMult(CFSVector & v1, CFSVector & v2)
{
  ENTER_FCN("Matrix::DyadicMult");
  
  Vector<TYPE> & vec1 = dynamic_cast<Vector<TYPE>& >(v1);
  Vector<TYPE> & vec2 = dynamic_cast<Vector<TYPE>& >(v2);
  
  Integer row = vec1.GetSize();
  Integer col = vec2.GetSize();
  
  this->Resize(row,col);

  for(Integer actRow=0; actRow<row; actRow++)
    for(Integer actCol=0; actCol < col; actCol++)
      data_[actRow][actCol] = vec1[actRow] * vec2[actCol];
}

template<class TYPE>
void Matrix<TYPE>::Invert (Matrix <TYPE> & inv) const
{   
 ENTER_FCN("Matrix::Invert");   

#ifdef CHECK_INITIALIZED
 if (size_row_ == 0 || size_col_ == 0) 
   Error("Undefined Matrix!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
 if (size_row_ != size_col_ ) 
   Error("No quadratic matrix!",__FILE__,__LINE__);
#endif

 TYPE det;
  switch (size_row_)
    {
    case 1: 
      inv.Resize(1);
      inv[0][0] = 1/data_[0][0];
      break;
    case 2:
      inv.Resize(2,2);
      inv[0][0] = data_[1][1];
      inv[0][1] = - data_[0][1];
      inv[1][0] = - data_[1][0];
      inv[1][1] = data_[0][0];
      this->Determinant(det);
      inv *= 1/det;
      break;

    case 3:
      // see Stöcker: "Taschenbuch Mathematischer Formeln und Moderner Verfahren" p.418
      inv.Resize(3,3);
      for(Integer i=0; i<3; i++)
	for(Integer j=0; j<3; j++)
	  inv[j][i] = Adjunct(i,j);      
      this->Determinant(det);
      inv *= 1/det;      
      break;
      
    default: 
      Error("Inversion not implemented fo dimension larger than 2!",__FILE__,__LINE__);
    }
}

template<> void Matrix<Complex>::Invert (Matrix <Complex> & inv) const
{
  ENTER_FCN("Matrix::Adjunct"); 
  Error("Matrix<Complex>::Invert: Not implemented!",__FILE__,__LINE__);
}

template<class TYPE>
TYPE Matrix<TYPE>::Adjunct (Integer i, Integer j) const
{
  ENTER_FCN("Matrix::Adjunct");  

#ifdef CHECK_INITIALIZED
 if (size_row_ == 0 || size_col_ == 0) 
   Error("Undefined Matrix!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
 if (size_row_ != size_col_ ) 
   Error("No quadratic matrix!",__FILE__,__LINE__);
#endif
     
 if (size_row_ != 3 ) Error("Matrix::Adjunct only implemented for matrix size 3!",__FILE__,__LINE__);  
  
  Vector<Integer> iVec(2);
  Vector<Integer> jVec(2);
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

  TYPE adj = (TYPE) pow(-1,i+j) * 
    ( data_[iVec[0]][jVec[0]] * data_[iVec[1]][jVec[1]] - 
      data_[iVec[0]][jVec[1]] * data_[iVec[1]][jVec[0]]);


  return adj;
  
};



// copies a submatrix at the position (row, col) into subMat, 
// the amount of copied elements depends on the size of subMat
template<class TYPE>
void Matrix<TYPE>::GetSubMatrix(Matrix<TYPE>& subMat, Integer startRow, Integer startCol) const
{
ENTER_FCN("Matrix::GetSubMatrix");

#ifdef CHECK_INITIALIZED
  if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
    Error("undefined matrix",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
  if (((subMat.size_row_ + startRow) > size_row_) || ((subMat.size_col_ + startCol) > size_col_) )
    Error("Submatrx to be written is to large! ",__FILE__,__LINE__);
#endif

  for(int actRow=0; actRow < subMat.size_row_; actRow++)
    for(int actCol=0; actCol < subMat.size_col_; actCol++)
      subMat[actRow][actCol] = data_[actRow + startRow][actCol + startCol];  
}



// overwrites the matrix elements at the position (row, col) with subMat
// in a rectangular (submatrix) way
template<class TYPE>
void Matrix<TYPE>::SetSubMatrix(Matrix<TYPE>& subMat, Integer startRow, Integer startCol)
{
  ENTER_FCN("Matrix::SetSubMatrix");
#ifdef CHECK_INITIALIZED
  if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
    Error("undefined matrix",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX
  if ((subMat.size_row_ + startRow > size_row_) || (subMat.size_col_ + startCol > size_col_) )
    Error("Submatrix to be read is to large! ",__FILE__,__LINE__);
#endif
  
  for(int actRow=0; actRow < subMat.size_row_; actRow++)
    for(int actCol=0; actCol < subMat.size_col_; actCol++)
      data_[actRow + startRow][actCol + startCol] = subMat[actRow][actCol];
}



/// converts a matrix into a vector, by appending successively all rows
template<class TYPE>
void Matrix<TYPE>::ConvertToVec_AppendRows(CFSVector & v) const
{
  ENTER_FCN("Matrix::ConvertToVec_AppendRows");
  
  Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

  vec.Resize(size_row_ * size_col_);
  
  for(int i=0; i < size_row_; i++)
    for(int j=0; j < size_col_; j++)
      vec[i*(size_row_-1) + j] = (*this)[i][j];
}

/// converts a matrix into a vector, by appending successively all colums
template<class TYPE>
void Matrix<TYPE>::ConvertToVec_AppendCols(CFSVector &v) const
{
  ENTER_FCN("Matrix::ConvertToVec_AppendCols");

  Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

  vec.Resize(size_row_ * size_col_);
  
  for(int actCol=0; actCol < size_col_; actCol++)
    for(int actRow=0; actRow < size_row_; actRow++)
      {
	vec[actCol*size_row_ + actRow] = data_[actRow][actCol];
      }
}


/// scales the diagonal elements of a  matrix by a factor
template<class TYPE>
void Matrix<TYPE>::ScaleDiagElems(TYPE factor) 
{
ENTER_FCN("Matrix::ScaleDiagElems");

#ifdef CHECK_INITIALIZED
 if (size_row_ == 0 || size_col_ == 0) 
   Error("Undefined Matrix!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX 
  if (size_row_ != size_col_ ) Error("No square- matrix!",__FILE__,__LINE__);
#endif

  Integer i;
  for (i = 0; i < size_row_; i++)
    data_[i][i] *= factor;
  
}


/// gets the diagonal elements of a  matrix in a one column matrix
template<class TYPE>
void Matrix<TYPE>::GetDiagInMatrix(Matrix<TYPE>& columnMat) 
{
ENTER_FCN("Matrix::GetDiagInMatrix");

#ifdef CHECK_INITIALIZED
 if (size_row_ == 0 || size_col_ == 0) 
   Error("Undefined Matrix!",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX 
  if (size_row_ != size_col_ ) Error("No square- matrix!",__FILE__,__LINE__);
#endif

  columnMat.Resize(size_row_, 1);
  
  Integer i;
  for (i = 0; i < size_row_; i++)
    columnMat [i][0] = data_[i][i];
  

}


// Two different methods for getting a double pointer to the data
// in the matrix:
// 1.) Double-Matrices: Here simply the pointer data_ is passed
// 2.) Complex-Matrices: A new array is created, which holds the 
//                       complex values as a serial sequence of
//                       Double values.
template<> Double* Matrix<Integer>::GetDoublePointer()
{
  ENTER_IFCN("Matrix::GetDoublePointer");
  Error("Matrix<Integer>::GetDoublePointer: Function not implemented!",__FILE__,__LINE__);
}

template<>
Double* Matrix<Double>::GetDoublePointer()
{
  ENTER_IFCN("Matrix::GetDoublePointer");
  return data_[0];
}

template<>
Double* Matrix<Complex>::GetDoublePointer()
{
  ENTER_IFCN("Matrix::GetDoublePointer");
  Double * help = new Double[size_col_ * size_row_ * 2];

  for (Integer i=0; i < size_row_*size_col_; i++)
    {
      help[2*i]   = data_[0][i].real();
      help[2*i+1] = data_[0][i].imag();
    }

  return help;
}

template<class TYPE> bool Matrix<TYPE>::IsSymmetric() {
  ENTER_FCN( "Matrix::IsSymmetric" );
  bool amSymm = true;
  for ( Integer i = 1; i < size_row_; i++ ) {
    for ( Integer j = i+1; j < size_col_; j++ ) {
	 if ( data_[i][j] != data_[j][i] ) {
	   amSymm = false;
	   break;
	 }
    }
  }
  return amSymm;
}


// Alternate version of symmetry checker, that will report
// asymmetries to standard output

// template<class TYPE> bool Matrix<TYPE>::IsSymmetric() {
//   ENTER_FCN( "Matrix::IsSymmetric" );
//   bool amSymm = true;
//   for ( Integer i = 1; i < size_row_; i++ ) {
//     for ( Integer j = i+1; j < size_col_; j++ ) {
//       if ( data_[i][j] != data_[j][i] ) {
// 	amSymm = false;
// 	(*debug) << " (" << i << "," << j << "): " << data_[i][j] << " <--> "
// 		 << data_[j][i] << "(abs.diff= " << (data_[i][j]-data_[j][i])
// 		 << ") "
// 		 << "(rel.diff= "
// 		 << (data_[i][j]-data_[j][i])/data_[i][j] << ")\n";
//       }
//     }
//   }
//   return amSymm;
// }

// explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Matrix<Integer>
#pragma instantiate Matrix<Double>
#pragma instantiate Matrix<Complex>
#endif

} // end of namespace
