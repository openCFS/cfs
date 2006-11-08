#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <iomanip>

#include "matrix.hh"
#include "Utils/vector.hh"

#include "Utils/boost-serialization.hh"

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
  Matrix<TYPE>::Matrix (const UInt nRows, const UInt nCols)
  {
    ENTER_FCN("Matrix::Matrix");
#ifdef CHECK_INDEX 
    if (nRows <= 0 || nCols <= 0) Error("invalid dimension",__FILE__,__LINE__);
#endif

    size_row_ = nRows;
    size_col_ = nCols;
    data_ = new TYPE* [size_row_];

    data_[0]=new TYPE[size_col_*size_row_];

    for (UInt k=1; k < size_row_; k++) 
      data_[k]=data_[k-1]+size_col_;
    Init();
  }


  template<class TYPE>
  Matrix<TYPE>::Matrix (const UInt nRows,const Vector<TYPE> * const x)
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

    UInt k,kk;
#ifdef CHECK_INDEX
    for (k=1; k < size_row_; k++)
    
      { if (x[k].size_!=size_col_)  
          Error(" Not all vectors for initialization have the same size",
                __FILE__,__LINE__);
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
//     if (x.size_row_ == 0 || x.size_col_ == 0)  
//       Error("undefined Matrix",__FILE__,__LINE__);
#endif

 
    size_row_ = x.size_row_;
    size_col_ = x.size_col_;

    if (size_row_ > 0 &&
        size_col_ > 0 ) {
      data_ = new TYPE * [size_row_];
      data_[0]=new TYPE[size_row_ * size_col_];
    } else {
      data_ = NULL;
    }
 
 
    UInt k;
 
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
  void Matrix<TYPE> :: Resize(const UInt nRows, const UInt nCols )
  {
    ENTER_FCN("Matrix::Resize");
  
    UInt k;
  
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
  
    // // initialize values to 0
//     if( init == true ) {
//       for ( k = 0; k < size_row_ * size_col_; k++) 
//         data_[0][k]=0;
    //}
  }


  template<class TYPE>
  void Matrix<TYPE>::Resize(const UInt col )
  {
    ENTER_FCN("Matrix::Resize");
    Resize(col,col);  
  }


#ifndef EXPR_TEMPLATES

  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator=(const Matrix<TYPE> &x)
  {
    ENTER_IFCN("Matrix::operator=");

#ifdef CHECK_INITIALIZED
    if (x.size_row_ == 0 || x.size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif  

    if (this == &x)  {
      return *this;
    }
  
    UInt k;
  
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
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif

    return *this;
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
    UInt k;
  
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
  
    UInt k;
    for ( k = 0; k < size_row_*size_col_; k++)
      z [0][k] = -data_[0][k];
  
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
    UInt k;
    for ( k = 0; k < size_row_ * size_col_; k++)
      data_ [0][k] -= x.data_ [0][k];
  
    return *this;
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
  
    UInt i;
    for (i = 0; i < size_row_*size_col_; i++)
      data_ [0][i] *= y;
  
    return *this;
  }

 template<class TYPE>
  Matrix<TYPE> & Matrix<TYPE>::operator*=(const Matrix<TYPE> &x)
  {   
    ENTER_IFCN("Matrix::operator*=");    
 
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
  
    UInt i,j; 
    for (i = 0; i < size_row_; i++)
      for (j = 0; j < x.size_col_; j++)
        {       
          a = data_ [i] [0] * x.data_ [0] [j];
          for (UInt k = 1; k < size_col_; k++)
            a += data_ [i] [k] * x.data_ [k] [j];
          z.data_ [i] [j] = a;
        }
  
    *this = z;
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
  
    UInt i;
    for (i = 0; i < size_row_*size_col_; i++)
      data_ [0][i] /= y;
  
    return *this;
  }


#endif // EXPR_TEMPLATES

  template<class TYPE>
  bool Matrix<TYPE>::operator== (const Matrix<TYPE> &x) const
  {
    ENTER_IFCN("Matrix::operator==");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      Error("undefined Matrix",__FILE__,__LINE__);
#endif
  
    UInt k;
  
    for (k = 0; k < size_row_*size_col_; k++)
      if (data_ [0][k] != x.data_[0][k]) return false;
  
    return true;
  }

  template<class TYPE>
  bool Matrix<TYPE>::operator!= (const Matrix<TYPE> &x) const
  {
    ENTER_IFCN("Matrix::operator!=");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      Error("undefined Matrix",__FILE__,__LINE__);
#endif 
  
    UInt k;
    for (k = 0; k < size_row_*size_col_; k++)
      if (data_ [0][k] != x.data_[0][k]) return false;
  
    return true;
  }

 

  // Perform a matrix-vector multiplication rvec = this*mvec
  template<class TYPE>
  void Matrix<TYPE>::Mult(const CFSVector & mvec, CFSVector & rvec) const
  {
    ENTER_IFCN("Matrix::Mult");
    Vector<TYPE> const & mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> & rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
    UInt size_mvec = mvec1.GetSize();
    UInt size_rvec = rvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
    if (size_mvec == 0) 
      Error("undefined Vector",__FILE__,__LINE__);
    if (size_rvec == 0) 
      Error("undefined Vector",__FILE__,__LINE__);
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mvec) 
      Error("incompatible dimension",__FILE__,__LINE__);
    if (size_row_ != size_rvec) 
      Error("incompatible dimension",__FILE__,__LINE__);
#endif
   
    UInt k,kk;
    rvec1.Init();
    for ( k = 0; k < size_row_; k++)
      for ( kk = 0; kk < size_col_; kk++)
        rvec1[k] += data_[k][kk]*mvec1[kk];

  }

  // // Perform a matrix-matrix multiplication rMat = this*mMat
  // template<class TYPE>
  // void Matrix<TYPE>::Mult(CFSMatrix & mMat, CFSMatrix & rMat)
  // {
  //   ENTER_IFCN("Matrix::Mult");
  //   Matrix<TYPE> & mMat1 = dynamic_cast<Matrix<TYPE> & >(mMat);
  //   Matrix<TYPE> & rMat1 = dynamic_cast<Matrix<TYPE>& >(rMat);
  //   
  //   UInt size_mMatRow = mMat1.GetSizeRow();
  //   UInt size_mMatCol = mMat1.GetSizeCol();
  // 
  //   UInt size_rMatRow = rMat1.GetSizeRow();
  //   UInt size_rMatCol = rMat1.GetSizeCol();
  //  
  // #ifdef CHECK_INITIALIZED
  //   if (size_row_ == 0 || size_col_ == 0) 
  //     Error("undefined Matrix",__FILE__,__LINE__);
  //   if (size_mMatRow == 0 || size_mMatCol==0) 
  //     Error("undefined Matrix",__FILE__,__LINE__);
  //   if (size_rMatRow == 0||size_rMatCol==0) 
  //     Error("undefined Matrix",__FILE__,__LINE__);
  // #endif
  // 
  // #ifdef CHECK_INDEX
  //   if (size_col_ != size_mMatRow) Error("incompatible dimension while matrix-matrix multiplication",__FILE__,__LINE__);
  //   if (size_row_ != size_rMatRow) Error("incompatible dimension while matrix-matrix multiplication",__FILE__,__LINE__);
  //   if (size_mMatCol != size_rMatCol) Error("incompatibel dimension while matrix-matrix multiplication",__FILE__,__LINE__);
  // #endif
  //    
  // //  Vector<TYPE> temp(1);
  // //  for (UInt i = 0; i < size_row_; i++)
  // //    for (UInt j = 0; j < size_mMatCol; j++)
  // //      {       
  // //   temp = data_[i][0] * mMat1[0][j];
  // //   for (UInt k = 1; k <size_mMatRow; k++)
  // //     temp[0] += data_[i][k] * mMat1[k][j];
  // //   rMat1[i][j] = temp[0];
  // //      }  
  // 
  //   for (UInt i = 0; i < size_row_; i++ ) {
  //     for (UInt j = 0; j < size_mMatCol; j++ ) {
  // 
  //       rMat1[i][j] = data_[i][0] * mMat1[0][j];
  // 
  //       for ( UInt k = 1; k < size_mMatRow; k++ ) {
  //         rMat1[i][j] += data_[i][k] * mMat1[k][j];
  //       }
  //     }
  //   }
  // } 




  template<class TYPE>
  void Matrix<TYPE>::AddRow (const Vector<TYPE> &x, const UInt pos)
  {
    ENTER_IFCN("Matrix::add_row");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif


    TYPE ** help=new TYPE*[size_row_+1];
    help[0]=new TYPE[size_col_ * (size_row_+1)];

    UInt k;
    for (k=1; k < size_row_+1; k++) 
      help[k]=help[k-1]+size_col_; 
  
    UInt i,ii;
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
  void Matrix<TYPE>::AddColumn (const Vector<TYPE> &x, const UInt pos)
  {
    ENTER_IFCN("Matrix::AddColumn");

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      Error("undefined Matrix",__FILE__,__LINE__);
#endif

    TYPE ** help=new TYPE*[size_row_];
    help[0]=new TYPE[(size_col_+1)*size_row_];
    UInt k;
    for (k=1; k < size_row_; k++) 
      help[k]=help[k-1]+size_col_+1;
  
    UInt i,ii;
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


  // **************
  //   operator<<
  // **************
  template<class S>
  std::ostream &operator<< ( std::ostream &out, const Matrix<S> &mat ) {

    ENTER_IFCN( "Matrix::operator<<" );

    // Set output format
    out.setf( std::ios::scientific );
    UInt oldPrec = out.precision();
    UInt newPrec = 4;
    out.precision( newPrec );

    for ( UInt i = 0; i < mat.GetSizeRow(); i++ ) {
      for ( UInt j = 0; j < mat.GetSizeCol(); j++ ) {
        out << std::setw( newPrec + 7 ) << mat[i][j] << " ";
      }
      out << std::endl;
    }

    // Restore old settings
    out.precision( oldPrec );

    return out;
  }

  template std::ostream &operator<<<UInt>    ( std::ostream &,
                                               const Matrix<UInt>    & );
  template std::ostream &operator<<<Integer> ( std::ostream &,
                                               const Matrix<Integer> & );
  template std::ostream &operator<<<Double>  ( std::ostream &,
                                               const Matrix<Double>  & );
  template std::ostream &operator<<<Complex> ( std::ostream &,
                                               const Matrix<Complex> & );


  template<class TYPE>
  void Matrix<TYPE>::Transpose (Matrix<TYPE> &transposedMat) const
  {
    ENTER_FCN("Matrix::Transpose");
    transposedMat.Resize(size_col_, size_row_);
 
    UInt i,j;
 
    for (i = 0; i < size_col_; i++)
      for (j = 0; j < size_row_; j++)
        transposedMat.data_[i][j] = data_[j] [i];
  }


  template<class TYPE>
  void Matrix<TYPE>::DirectSolve(CFSVector & x1, CFSVector & b1)
  {
    ENTER_FCN("Matrix::DirectSolve");

    Vector<TYPE> & x = dynamic_cast<Vector<TYPE>& >(x1);
    Vector<TYPE> & b = dynamic_cast<Vector<TYPE>& >(b1);
    
    Integer nmat = size_row_-1;
    Integer i, j, k, k1;
    
    //  the Gauss elimination 
    for (k=0; k<=nmat-1; ++k)
      {
        k1 = k + 1;
        for (i=k1; i<=nmat; ++i)
          {
            if (data_[k][k] != 0.0)
              {
                data_[i][k] /= data_[k][k];
                for (j=k1; j<=nmat; ++j)
                  data_[i][j] -= data_[i][k] * data_[k][j];
              }
            else
              {
                std::cerr<<"\n Matrix::DirectSolve Step " << k <<std::endl;
                std::cerr<<"\n The element at position ("<< k << "," << k << ") of the matrix is zero. "<<std::endl;
                //                std::exit(0);
              }
          }
      }

    // solve Ly = b by forward substitution 

   
    Vector<TYPE> y(b.GetSize());

    for (i=0; i<=nmat; ++i)
      { 
        y[i] = b[i];
        for (j=0; j<=i-1; ++j)
          y[i] -= data_[i][j] * y[j];
      }

    // solve Ux = y backward substitution
    
    for (i=nmat; i>=0; --i)
      {
        x[i] = y[i];
        for (j=nmat; j>=i+1; --j)
          x[i] -= data_[i][j] * x[j];
        x[i] /= data_[i][i];
      }
  }


#ifdef USE_LAPACK
  // Compile OLAS and CFS++ with USE_LAPACK
  template<>
  void Matrix<Complex>::solveWithLapack(Matrix<Complex> & b1,
                                        lapackSysMatType & LAPACK_MATRIX_TYPE)
  {
    ENTER_FCN("Matrix::solveWithLapack");

    if(size_row_!=size_col_)
      Error("Matrix is not quadratic, no solver! ",__FILE__,__LINE__);

    char lp_matType;
    lp_matType='L';
    
    Integer lp_nrRHS, lp_loadDim, lp_loadDimRHS, lp_info, lp_lwork,lp_dim;
    Integer lp_lda, lp_ldb;
    lp_nrRHS=b1.GetSizeCol();
    lp_dim=size_row_;
    lp_lda=size_row_;
    lp_ldb=size_row_;

    Integer *lp_interchanges;
    
    F77complex16 *lp_rhsVecf77;
    F77complex16 *lp_sysVecf77;
    F77complex16 *lp_workf77;
    F77complex16 auxVal2;
    
    Vector<Complex> lp_sysVec, lp_work;
    lp_sysVec.Resize(size_row_*size_row_);
    
      // copy values from system and RHS - Matrix into vector
    for (UInt i=0;i<size_row_;i++)
      for (UInt j=0;j<size_row_;j++){
        lp_sysVec[i+j*size_row_]=data_[i][j];
      }
    
    Vector<Complex> lp_rhsVec;
    lp_rhsVec.Resize(b1.GetSizeRow()*b1.GetSizeCol());
 
    for (UInt i=0;i<b1.GetSizeRow();i++)
      for (UInt j=0;j<b1.GetSizeCol();j++){
        lp_rhsVec[i+j*b1.GetSizeRow()]=b1[i][j];
      }
    std::cout<<lp_rhsVec<<std::endl;

    lp_rhsVecf77 = new F77complex16[size_row_*lp_nrRHS];
    lp_interchanges = new int[size_row_*size_row_];
    lp_workf77 = new F77complex16[size_row_*size_row_];
    lp_sysVecf77 = new F77complex16[size_row_*size_row_];
   

    // Convert CFS++ Vector<Complex> to Vector<F77complex16>
    for ( UInt count = 0; count < size_row_*b1.GetSizeCol(); count++ ) {
      CC2F77( lp_rhsVec[count], auxVal2 );
      lp_rhsVecf77[count] = auxVal2;
      }
    
    for (UInt count = 0; count < size_row_*size_row_; count++ ) {
      CC2F77( lp_sysVec[count], auxVal2 );
      lp_sysVecf77[count] = auxVal2;
    }
    
    for (UInt count=0; count < lp_work.GetSize();count++){
      CC2F77(lp_work[count], auxVal2);
      lp_workf77[count] = auxVal2;
    }
    
    switch (LAPACK_MATRIX_TYPE){
      
    case ZGESV:
      // solves systems with general system matrix
      zgesv_(&lp_dim , &lp_nrRHS, lp_sysVecf77, &lp_lda, 
             lp_interchanges, lp_rhsVecf77, &lp_ldb, &lp_info);

      if ( lp_info != 0 ) {
        Error( "ZGESV reports invalid input parameter", __FILE__, __LINE__ );
      }    
      break;
    case ZSYSV:
      
      lp_lwork=192;
      // solves systems with symmetric system matrix
      zsysv_(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77, 
             &lp_lda, lp_interchanges, lp_rhsVecf77, &lp_ldb,
             lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        Error( "ZSYSV reports invalid input parameter", __FILE__, __LINE__ );
      }
      break;
    case ZHESV:
      lp_lwork=192;
      // solves systems with hermitian system matrix
      zhesv_(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77,
             &lp_lda, lp_interchanges,lp_rhsVecf77, &lp_ldb, 
             lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        Error( "ZHESV reports invalid input parameter", __FILE__, __LINE__ );
      }
      break;

    default:
      std::cout<<LAPACK_MATRIX_TYPE<<std::endl;
      Error( " the chosen routine is not yet implemented in solveWithLapack",
             __FILE__, __LINE__ );


    } // matches switch ()...
      
          
     //reconvert Fortran77 -> CFS ++ datatypes
    for ( UInt count = 0; count < size_row_*b1.GetSizeCol(); count++ ) 
      F772CC( lp_rhsVecf77[count], lp_rhsVec[count] );
    std::cout<<lp_rhsVec<<std::endl;
      
    for ( UInt count = 0; count < size_row_*size_row_; count++ ) 
      F772CC( lp_sysVecf77[count], lp_sysVec[count]);
    std::cout<<lp_sysVec<<std::endl;
      
    for (UInt count=0; count < lp_work.GetSize();count++)
      F772CC(lp_workf77[count], lp_work[count]);
      
    // Writes result into b1
    for (UInt i=0;i<size_row_;i++)
      for (UInt j=0;j<b1.GetSizeCol();j++)
        b1[i][j]=lp_rhsVec[i+j*b1.GetSizeRow()];
  
    delete[] (lp_rhsVecf77);
    delete[] (lp_interchanges);
    delete[] (lp_sysVecf77);
    delete[] (lp_workf77 );
    
  }
#endif

#ifdef USE_LAPACK
  // Compile OLAS and CFS++ with USE_LAPACK
  template<>
  void Matrix<Complex>::eigenvaluesWithLapack(Vector<Double> & lp_w)
  {
    ENTER_FCN("Matrix::eigenvaluesWithLapack");
    // computes all eigenvalues of a complex hermitian matrix

    char lp_jobz='V';
    char lp_uplo='L';
    Integer lp_N=size_row_;                
    Integer lp_lda=size_row_;
      
    // array contains ev in ascending order
    //    Vector<Double> lp_w;
    lp_w.Resize(size_row_);
    lp_w.Init();
    Integer lp_lworkf77=99;
      
    // workspace array - complex 16 array
    Vector<Complex> lp_work;
    lp_work.Resize(lp_lworkf77);
    lp_work.Init();
      
    // workspace array - double precission
    Vector<Double> lp_rwork;
    lp_rwork.Resize(3*size_row_-2);
    lp_rwork.Init();
      
    Integer lp_infof77;
    F77complex16 auxValC;
    F77real8 auxValR;

    F77complex16 * lp_af77 = new F77complex16[size_row_*size_row_];
    F77real8 * lp_wf77 = new F77real8[size_row_];
    F77complex16 * lp_workf77 = new F77complex16[99];
    F77real8 * lp_rworkf77 = new F77real8[3*size_row_-2];
      
    // Convert CFS++ Vector<Complex> to Vector<F77complex16>
    for ( UInt count = 0; count < size_row_; count++ ) 
      for ( UInt countC = 0; countC <size_row_; countC++ ) {
        CC2F77( data_[count][countC], auxValC );
        lp_af77[countC*size_row_+count] = auxValC;
      }
    
    for ( UInt count = 0; count < size_row_; count++ ) {
      CC2F77( lp_w[count], auxValR );
      lp_wf77[count] = auxValR;
    }
    
    for (UInt count=0; count < lp_rwork.GetSize();count++){
      CC2F77(lp_rwork[count], auxValR);
      lp_rworkf77[count] = auxValR;
    }
    
    zheev_( &lp_jobz, &lp_uplo, &lp_N, lp_af77, &lp_lda, lp_wf77, 
            lp_workf77, &lp_lworkf77, lp_rworkf77 ,&lp_infof77); 
    
    // reconvert f772C++
    for (UInt count=0; count < lp_work.GetSize();count++)
      F772CC(lp_workf77[count], lp_work[count]);
    
    for ( UInt count = 0; count < size_row_; count++ ) 
      F772CC( lp_wf77[count], lp_w[count] );
    
    delete lp_workf77;
    delete lp_rworkf77;
    delete lp_af77;
    delete lp_wf77;
    
  }

#endif
  

  template<class TYPE>
  void Matrix<TYPE>::DyadicMult(const CFSVector & v1)
  {
    ENTER_FCN("Matrix::DyadicMult");
    DyadicMult(v1, v1);
  }



  template<class TYPE>
  void Matrix<TYPE>::DyadicMult(const CFSVector & v1, const CFSVector & v2)
  {
    ENTER_FCN("Matrix::DyadicMult");
  
    Vector<TYPE> const & vec1 = dynamic_cast<const Vector<TYPE>& >(v1);
    Vector<TYPE> const & vec2 = dynamic_cast<const Vector<TYPE>& >(v2);
  
    UInt row = vec1.GetSize();
    UInt col = vec2.GetSize();
  
    this->Resize(row,col);

    for(UInt actRow=0; actRow<row; actRow++)
      for(UInt actCol=0; actCol < col; actCol++)
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
        for(UInt i=0; i<3; i++)
          for(UInt j=0; j<3; j++)
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
    ENTER_FCN("Matrix::Invert"); 
    Error("Matrix<Complex>::Invert: Not implemented!",__FILE__,__LINE__);
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::Adjunct (UInt i, UInt j) const
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
    UInt runningIndexI = 0;
    UInt runningIndexJ = 0;
  
    for (UInt actI = 0; actI<=2; actI++)
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
    
  }
  
  template<class TYPE>
  Matrix<Double> Matrix<TYPE>::GetPart(  DataType part ) const {
    Error( "Matrix::GetPart: Only Implemented for Real and Complex matrices!", 
           __FILE__, __LINE__ );
    Matrix<Double> temp;
    return temp;
  }

  
  template<>
  Matrix<Double> Matrix<Double>::GetPart(  DataType part ) const {
    ENTER_FCN( "Matrix<Double>::GetPart" )
    
    if ( part != REAL ) {
      Error("Matrix<Double>::GetPart: Only possible for REAL part.",__FILE__, __LINE__ );
    }
      return *this;
    
  }

  template<>
  Matrix<Double> Matrix<Complex>::GetPart(  DataType part ) const {
    ENTER_FCN( "Matrix<Complex>::GetPart" );
    
    Matrix<Double> ret;
    if ( part == REAL ) {
      ret.Resize( size_row_, size_col_ );
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          ret[iRow][iCol]  = data_[iRow][iCol].real();
        }
      }
    } else if ( part == IMAG ) {
      ret.Resize( size_row_, size_col_ );
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          ret[iRow][iCol]  = data_[iRow][iCol].imag();
        }
      }
    } else {
      Error("Matrix<Complex>::GetPart: Only possible for REAL or IMAG part!" , 
	    __FILE__, __LINE__ );
    }
    
    return ret;
  }
  

  template<class TYPE>
  void Matrix<TYPE>::SetPart( DataType part, const Matrix<Double> & partMatrix ) {
    Error( "Matrix::SetPart: Only Implemented for Real and Complex matrices!", 
           __FILE__, __LINE__ );
  }
  
  template<>
  void Matrix<Double>::SetPart( DataType part, const Matrix<Double> & partMatrix ) {
    
    if ( size_col_ != partMatrix.GetSizeCol() ||
	 size_row_ != partMatrix.GetSizeRow () ) {
      Error( "Matrix<Double>::SetPart: Dimension of matrices do not match!", __FILE__, __LINE__ );
    }
 
    if ( part != REAL ) {
      Error( "Matrix<Double>::SetPart: Only possible for REAL part.", __FILE__, __LINE__ );
    }
    *this = partMatrix;
  }

  template<>
  void Matrix<Complex>::SetPart( DataType part, const Matrix<Double> & partMatrix ) {

    if ( size_col_ != partMatrix.GetSizeCol() ||
         size_row_ != partMatrix.GetSizeRow () ) {
      Error( "Matrix<Complex>::SetPart: Dimension of matrices do not match!", __FILE__, __LINE__ );
    }
        
    if ( part == REAL ) {
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          data_[iRow][iCol]  = Complex( partMatrix[iRow][iCol], 
                                        data_[iRow][iCol].imag() );
        }
      }
    } else if ( part == IMAG ) {
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          data_[iRow][iCol]  = Complex( data_[iRow][iCol].real(),
                                        partMatrix[iRow][iCol] );
        }
      }
    } else {
      Error( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!",
	     __FILE__, __LINE__ );
    }
  }
 

  // copies a submatrix at the position (row, col) into subMat, 
  // the amount of copied elements depends on the size of subMat
  template<class TYPE>
  void Matrix<TYPE>::GetSubMatrix(Matrix<TYPE>& subMat, UInt startRow, UInt startCol) const
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

    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        subMat[actRow][actCol] = data_[actRow + startRow][actCol + startCol];  
  }



  // overwrites the matrix elements at the position (row, col) with subMat
  // in a rectangular (submatrix) way
  template<class TYPE>
  void Matrix<TYPE>::SetSubMatrix(const Matrix<TYPE>& subMat, UInt startRow, UInt startCol)
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
  
    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        data_[actRow + startRow][actCol + startCol] = subMat[actRow][actCol];
  }



  // adds subMat to the matrix elements at the position (row, col)
  // in a rectangular (submatrix) way
  template<class TYPE>
  void Matrix<TYPE>::AddSubMatrix(const Matrix<TYPE>& subMat, UInt startRow, UInt startCol)
  {
    ENTER_FCN("Matrix::AddSubMatrix");
#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      Error("undefined matrix",__FILE__,__LINE__);
#endif
  
#ifdef CHECK_INDEX
    if ((subMat.size_row_ + startRow > size_row_) || (subMat.size_col_ + startCol > size_col_) )
      Error("Submatrix to be read is to large! ",__FILE__,__LINE__);
#endif
  
    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        data_[actRow + startRow][actCol + startCol] += subMat[actRow][actCol];
  }



  /// converts a matrix into a vector, by appending successively all rows
  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_AppendRows(CFSVector & v) const
  {
    ENTER_FCN("Matrix::ConvertToVec_AppendRows");
  
    Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

    vec.Resize(size_row_ * size_col_);
  
    for( UInt i=0; i < size_row_; i++)
      for( UInt j=0; j < size_col_; j++)
        vec[i*(size_col_) + j] = (*this)[i][j];
  }

  /// converts a matrix into a vector, by appending successively all colums
  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_AppendCols(CFSVector &v) const
  {
    ENTER_FCN("Matrix::ConvertToVec_AppendCols");

    Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

    vec.Resize(size_row_ * size_col_);
  
    for( UInt actCol=0; actCol < size_col_; actCol++)
      for( UInt actRow=0; actRow < size_row_; actRow++)
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

    UInt i;
    for (i = 0; i < size_row_; i++)
      data_[i][i] *= factor;
  
  }


  /// gets the diagonal elements of a  matrix in a one column matrix
  template<class TYPE>
  void Matrix<TYPE>::GetDiagInMatrix(Matrix<TYPE>& columnMat) const
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
    columnMat.Init();
  
    UInt i;
    for (i = 0; i < size_row_; i++)
      columnMat [i][0] = data_[i][i];
  

  }


  template<class TYPE> bool Matrix<TYPE>::IsSymmetric() const {
    ENTER_FCN( "Matrix::IsSymmetric" );
    bool amSymm = true;
    for ( UInt i = 1; i < size_row_; i++ ) {
      for ( UInt j = i+1; j < size_col_; j++ ) {
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
  //   for ( UInt i = 1; i < size_row_; i++ ) {
  //     for ( UInt j = i+1; j < size_col_; j++ ) {
  //       if ( data_[i][j] != data_[j][i] ) {
  //      amSymm = false;
  //      (*debug) << " (" << i << "," << j << "): " << data_[i][j] << " <--> "
  //               << data_[j][i] << "(abs.diff= " << (data_[i][j]-data_[j][i])
  //               << ") "
  //               << "(rel.diff= "
  //               << (data_[i][j]-data_[j][i])/data_[i][j] << ")\n";
  //       }
  //     }
  //   }
  //   return amSymm;
  // }

#ifdef __GNUC__
  template class Matrix<Double>;
  template class Matrix<Integer>;
  template class Matrix<UInt>;
  template class Matrix<Complex>;
  template class Matrix<bool>;
#endif


// explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Matrix<UInt>
#pragma instantiate Matrix<Integer>
#pragma instantiate Matrix<Double>
#pragma instantiate Matrix<Complex>
#pragme instantiate Matrix<bool>
#endif

} // end of namespace

#include "Utils/boost-serialization.hh"
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Double>, "CoupledField_Matrix_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Complex>, "CoupledField_Matrix_Complex")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Integer>, "CoupledField_Matrix_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::UInt>, "CoupledField_Matrix_UInt")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<bool>, "CoupledField_Matrix_Bool")
    

