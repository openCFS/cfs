// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "MatVec/denseMatrix.hh"
#include "MatVec/matrixLapackSupport.hh"
#include "MatVec/opdefs.hh"
#include "MatVec/vector.hh"
#include "Utils/StdVector.hh"
#include "Utils/tools.hh"
#include "boost/lexical_cast.hpp"
#include "boost/tokenizer.hpp"
#include "boost/type_traits/is_same.hpp"
#include "def_build_type_options.hh"
#include "matrix.hh"

using boost::tokenizer;
using boost::char_separator;

namespace CoupledField {

class SingleVector;
}  // namespace CoupledField

namespace CoupledField
{      

  template<class TYPE>
  Matrix<TYPE>::Matrix () :
    DenseMatrix(),
    size_row_(0),
    size_col_(0),
    data_(NULL)
  { }


  template<class TYPE>
  Matrix<TYPE>::Matrix (const UInt nRows, const UInt nCols) :
    DenseMatrix(),
    size_row_(nRows),
    size_col_(nCols),
    data_(new TYPE* [size_row_])
  {
#ifdef CHECK_INDEX 
    if (nRows <= 0 || nCols <= 0) EXCEPTION("invalid dimension");
#endif

    data_[0]=new TYPE[size_col_*size_row_];

    for (UInt k=1; k < size_row_; k++) 
      data_[k]=data_[k-1]+size_col_;

    Init();
  }


  template<class TYPE>
  Matrix<TYPE>::Matrix (const UInt nRows,const Vector<TYPE> * const x) :
    DenseMatrix(),
    size_row_(nRows),
    size_col_(x[0].size_),
    data_(new TYPE* [size_row_])
  {
#ifdef CHECK_INDEX
    if (nRows <= 0) EXCEPTION("invalid dimension");
    if (size_col_ == 0) EXCEPTION("invalid dimension");
#endif 

    UInt k,kk;
#ifdef CHECK_INDEX
    for (k=1; k < size_row_; k++)
    
      { if (x[k].size_!=size_col_)  
          EXCEPTION(" Not all vectors for initialization have the same size" );
      }
#endif
  
    for (k=0; k < size_row_; k++)
      for (kk=0; kk<size_col_; kk++)
        data_[k][kk]=x[k][kk];
  }

  template<class TYPE>
  Matrix<TYPE>::Matrix (const Matrix<TYPE> &x) :
    DenseMatrix(),
    size_row_(x.size_row_),
    size_col_(x.size_col_),
    data_(NULL)
  {
    // We are able to copy empty matrices!
    if (size_row_ > 0 && size_col_ > 0 )
    {
      data_ = new TYPE * [size_row_];
      data_[0]=new TYPE[size_row_ * size_col_];
    }
 
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k]=x.data_[0][k];
    
    for(UInt k = 1; k < size_row_; ++k)
      data_[k]=data_[k-1]+size_col_;        
  }

  template<class TYPE>
  Matrix<TYPE>::~Matrix ()
  {
    if (data_ != NULL)
    {
      delete[] data_[0];
      delete[] data_;
      data_= NULL;
    }
  }


  template<class TYPE>
  std::string Matrix<TYPE>::ToXMLFormat(const std::string& name, int n_offset) const
  {
    std::string offset(std::max(n_offset, 0), ' ');

    std::ostringstream os;

    bool is_complex = boost::is_same<TYPE, std::complex<double> >::value;

    os << "<" << name << " dim1=\"" << size_row_ << "\" dim2=\"" << size_col_ << "\">";
    os << std::endl << offset << "  " << (is_complex ? "<complex>" : "<real>");

    for(unsigned int r = 0; r < size_row_; ++r)
    {
      os << std::endl << offset << "    ";
      for(unsigned int c = 0; c < size_col_; ++c)
      {
        os << std::scientific;
        os.precision(6);
        os.width(13);
        os << (IsNoise(data_[r][c]) ? 0.0 : data_[r][c]);
        if(c < size_col_ - 1) os << " ";
      }
    }
    os << std::endl << offset << "  " << (is_complex ? "</complex>" : "</real>");
    os << std::endl << offset << "</" << name << ">";

    return os.str();
  }

  template<class TYPE>
  std::string Matrix<TYPE>::ToString(const int level, const bool newline) const
  {
    std::ostringstream os;

    switch(level)
    {
      case -1:
        for(UInt j = 0; j < size_row_; j++)
          for(UInt i = 0; i < size_col_; i++)
            os << data_[j][i] << (j == size_row_ - 1 && i == size_col_ - 1 ? "" : " "); // space not for last element

        os << "]";
      break;

      case 0:
        for(UInt j = 0; j < size_row_; j++)
        {
          os << j << " : [";

          for(UInt i = 0; i < size_col_; i++)
            os << data_[j][i] << (i < size_col_-1 ? " " : "");

          os << "]";

          if(newline) os << std::endl;
          else os << " ";
        }
        break;

      case 1:
        for(UInt j = 0; j < size_row_; j++)
        {
          //os << j << ":"
          os << "[";

          for(UInt i = 0; i < size_col_; i++)
            os << data_[j][i] << (i < size_col_-1 ? " " : "");

          os << "]\n";
        }
        break;

      case 2:
        os << "[";
        for(UInt j = 0; j < size_row_; j++)
        {
          for(UInt i = 0; i < size_col_; i++)
            os << data_[j][i] << (i < size_col_-1 ? " " : "");
          if(j < size_row_-1)
            os << "; ";
        }
        os << "]";
        break;


      default:
        os << "size_row=" << size_row_ << " size_col=" << size_col_;
        if(size_row_ > 0 && size_col_ > 0)
        {
          // the min/max for complex is the real part and we cannot compare complex numbers anyway
          Double min = static_cast<Complex>(data_[0][0]).real();
          Double max = static_cast<Complex>(data_[0][0]).real();

          for(UInt j = 0; j < size_row_; j++)
            for(UInt i = 0; i < size_col_; i++)
            {
              min = std::min(min, static_cast<Complex>(data_[j][i]).real());
              max = std::max(max, static_cast<Complex>(data_[j][i]).real());
            }
          os << " min=" << min << " max=" << max;
        }
    }

    return os.str();
  }

  template<class TYPE>
  unsigned int Matrix<TYPE>::ParseLineHelper(const std::string& input, StdVector<TYPE>& out)
  {
    // the rows are in the form "1 2 3] : " or "1.2 2.5 3.5] : " or "(1,1) (0,-2)] : "

    // out will be resized in all but the first calls. So Push_back is not that expensive!
    char_separator<char> sep(" ]:"); // spaces, closing and colon
    tokenizer<char_separator<char> > tokens(input, sep);

    for(tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it)
    {
      try
      {
        out.Push_back(boost::lexical_cast<TYPE>(*it));
      }
      catch(boost::bad_lexical_cast &)
      {
        EXCEPTION("Cannot cast to value while parsing matrix: '" << *it << "'");
      }
    }

    if(out.GetSize() == 0) EXCEPTION("cannot interpret matrix row as data: '" << input << "'");

    return out.GetSize();
  }


  template<class TYPE>
  void Matrix<TYPE>::Parse(const std::string& input)
  {
    // we don't know the number or rows and cols in advance.
    StdVector<std::string> rows;
    char_separator<char> sep("["); // ignore the row number and colon
    tokenizer<char_separator<char> > tokens(input, sep);

    for(tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it)
      rows.Push_back(*it);

    // the rows are in the form "1 2 3] : " or "1.2 2.5 3.5] : " or "(1,1) (0,-2)] : "
    if(rows.GetSize() == 0) EXCEPTION("cannot interpret data as matrix with data: '" << input << "'");
    StdVector<TYPE> row;
    // check first row
    unsigned int n_cols = ParseLineHelper(rows[0], row); // row is now set to the proper size

    // we have the information, store the data!
    this->Resize(rows.GetSize(), n_cols);

    for(unsigned int r = 0; r < rows.GetSize(); r++)
    {
      row.Clear();
      row.Reserve(n_cols); // I know no Clear() Reserve() combination, but it's not too expensive

      ParseLineHelper(rows[r], row); // row is now set to the proper size
      if(row.GetSize() != n_cols)
        EXCEPTION("matrix has inconsistent number or columns within row " << (r+1) << ": '" << input << "'");

      for(unsigned int c = 0; c < row.GetSize(); c++)
        data_[r][c] = row[c];
    }
  }


  template<class TYPE>
  void Matrix<TYPE>::Resize(const UInt nRows, const UInt nCols )
  {
    if(nRows != size_row_ || nCols != size_col_)
    {
      // set the size to requested values
      size_row_ = nRows; 
      size_col_ = nCols;
      
      if (data_ != NULL)
      {
        delete [] data_[0];
        delete [] data_;
      }

      data_    = new TYPE*[size_row_];
      data_[0] = new TYPE [size_row_*size_col_];

      for (UInt k = 1; k < size_row_; ++k) 
        data_[k] = data_[k-1] + size_col_;

    }
  }


  template<class TYPE>
  void Matrix<TYPE>::Resize(const UInt col )
  {
    Resize(col,col);  
  }

  template<class TYPE>
  void Matrix<TYPE>::Resize(const Matrix<TYPE>& other)
  {
    Resize(other.size_row_, other.size_col_);  
  }

#ifndef EXPR_TEMPLATES

  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator=(const Matrix<TYPE> &x)
  {
    // allows to copy an empty matrix. !

    // Note! it shall be possible to copy an empty matrix!
//#ifdef CHECK_INITIALIZED
//    if (x.size_row_ == 0 || x.size_col_ == 0) 
//      EXCEPTION("undefined Matrix");
//#endif  

    if (this == &x)
      return *this;

    // set the size in any case to the size of the assigned matrix
    Resize(x.size_row_, x.size_col_);
  
    // copy the entries
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] = x.data_[0][k];
  
    return *this;
  }

  template<class TYPE>
  Matrix<TYPE> Matrix<TYPE>::operator+ () const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif

    return *this;
  }

 

  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator+=(const Matrix<TYPE> &x)
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif
  
#ifdef CHECK_INDEX  
    if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
      EXCEPTION("incompatible dimension for +"); 
#endif
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] += x.data_[0][k];
  
    return *this;
  }

  template<class TYPE>
  Matrix<TYPE> Matrix<TYPE>::operator-() const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 )
      EXCEPTION("undefined Matrix");
#endif

    Matrix<TYPE> z(size_row_,size_col_);
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        z[0][k] = -data_[0][k];
    
    return z;
  }

 

  template<class TYPE>
  Matrix<TYPE> & Matrix<TYPE>::operator-=(const Matrix<TYPE> &x)
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif
  
#ifdef CHECK_INDEX  
    if (size_row_ != x.size_row_ || size_col_ != x.size_col_)
      EXCEPTION("incompatible dimension for +"); 
#endif
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] -= x.data_[0][k];
    
    return *this;
  }


  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator/= (const TYPE x)
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] /= x;

    return *this;
  }


#endif // EXPR_TEMPLATES
  
  template<class TYPE>
  Matrix<TYPE> &Matrix<TYPE>::operator*= (const TYPE x)
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
#endif

    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] *= x;

    return *this;
  }

 template<class TYPE>
  Matrix<TYPE> & Matrix<TYPE>::operator*=(const Matrix<TYPE> &x)
  {   
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX  
    if (size_col_ != x.size_row_)
      EXCEPTION("incompatible dimension");
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
  bool Matrix<TYPE>::operator== (const Matrix<TYPE> &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if(data_ [0][k] != x.data_[0][k]) return false;
  
    return true;
  }

  template<class TYPE>
  bool Matrix<TYPE>::operator!= (const Matrix<TYPE> &x) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || 
        x.size_row_ == 0 || x.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif 
  
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if (data_ [0][k] != x.data_[0][k]) return true;
  
    return false;
  }
  
  /** Adds the multiple of another matrix */
  template<class TYPE>
  void Matrix<TYPE>::Add(const TYPE factor, const Matrix<TYPE>& mat)
  {
    const Matrix<TYPE>& other_mat = dynamic_cast<const Matrix<TYPE> & >(mat);
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0 || other_mat.size_row_ == 0 || other_mat.size_col_ == 0)
      EXCEPTION("undefined Matrix");
#endif

#ifdef CHECK_INDEX
    if (size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("matrices do not match");
#endif
    
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] += factor * other_mat.data_[0][k];
  }

  /** Adds the multiple of the transpose of another matrix */
   template<class TYPE>
   void Matrix<TYPE>::AddT(const TYPE factor, const Matrix<TYPE>& mat)
   {
     const Matrix<TYPE>& other_mat = dynamic_cast<const Matrix<TYPE> & >(mat);
 #ifdef CHECK_INITIALIZED
     if (size_row_ == 0 || size_col_ == 0 || other_mat.size_row_ == 0 || other_mat.size_col_ == 0)
       EXCEPTION("undefined Matrix");
 #endif

 #ifdef CHECK_INDEX
     if (size_row_ != other_mat.size_col_ || size_col_ != other_mat.size_row_)
       EXCEPTION("matrices do not match");
 #endif

     for(UInt k = 0, s = size_row_ ; k < s; ++k)
       for (UInt kk =0, ss = size_col_; kk < ss; ++kk)
         data_[k][kk] += factor * other_mat.data_[kk][k];
   }


  /** Assigns a multiple of another matrix */
  template<class TYPE>
  void Matrix<TYPE>::Assign(const Matrix<TYPE>& other_mat, TYPE factor)
  {
#ifdef CHECK_INDEX
    if(size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_) 
      EXCEPTION("matrices do not match");
#endif
    
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] = factor * other_mat.data_[0][k];
  }
  
  // Perform a matrix-vector multiplication rvec = this*mvec
  template<class TYPE>
  void Matrix<TYPE>::Mult(const SingleVector & mvec, SingleVector & rvec) const
  {
    Vector<TYPE> const & mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> & rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
#if defined CHECK_INITIALIZED || defined CHECK_INDEX
    UInt size_mvec = mvec1.GetSize();
    UInt size_rvec = rvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mvec == 0) 
      EXCEPTION("undefined Vector");
    if (size_rvec == 0) 
      EXCEPTION("undefined Vector");
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mvec) 
      EXCEPTION("incompatible dimension");
    if (size_row_ != size_rvec) 
      EXCEPTION("incompatible dimension");
#endif

#endif
    
    UInt k,kk;
    rvec1.Init();
    for ( k = 0; k < size_row_; k++)
      for ( kk = 0; kk < size_col_; kk++)
        rvec1[k] += data_[k][kk]*mvec1[kk];
  }


  template<class TYPE>
  TYPE Matrix<TYPE>::FrobeniusProduct(const Matrix<TYPE>& other_mat) const
  {
#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("incompatible dimension");
#endif

    TYPE result(0);

    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += OpType<TYPE>::dotProduct(data_[0][k], other_mat.data_[0][k]);

    return result;
  }

  template<>
  Matrix<double> Matrix<double>::EntryMult(const Matrix<double>& other_mat) const
  {
#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(size_row_ != other_mat.size_row_ || size_col_ != other_mat.size_col_)
      EXCEPTION("incompatible dimension");
#endif

    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      data_[0][k] *= other_mat.data_[0][k];

    return *this;
  }


  // Perform a matrix-vector multiplication rvec = this*mvec via the inner product
  template<class TYPE>
  void Matrix<TYPE>::MultInner(const SingleVector & mvec, SingleVector & rvec) const
  {
    Vector<TYPE> const & mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> & rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
#if defined CHECK_INITIALIZED || defined CHECK_INDEX
    UInt size_mvec = mvec1.GetSize();
    UInt size_rvec = rvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mvec == 0) 
      EXCEPTION("undefined Vector");
    if (size_rvec == 0) 
      EXCEPTION("undefined Vector");
#endif

#ifdef CHECK_INDEX
    if (size_col_ != size_mvec) 
      EXCEPTION("incompatible dimension");
    if (size_row_ != size_rvec) 
      EXCEPTION("incompatible dimension");
#endif

#endif

    UInt k,kk;
    rvec1.Init();
    for ( k = 0; k < size_row_; k++)
      for ( kk = 0; kk < size_col_; kk++)
        rvec1[k] += OpType<TYPE>::dotProduct(data_[k][kk],mvec1[kk]);
  }
  
  // Perform a matrix-vector multiplication rvec = transpose(this)*mvec
  template<class TYPE>
  void Matrix<TYPE>::MultT(const SingleVector &mvec, SingleVector &rvec) const
  {
    Vector<TYPE> const &mvec1 = dynamic_cast<const Vector<TYPE>& >(mvec);
    Vector<TYPE> &rvec1 = dynamic_cast<Vector<TYPE>& >(rvec);
  
#if defined CHECK_INITIALIZED || defined CHECK_INDEX
    UInt size_mvec = mvec1.GetSize();
 
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("undefined Matrix");
    if (size_mvec == 0) 
      EXCEPTION("undefined input Vector");
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_mvec) 
      EXCEPTION("incompatible dimension");
#endif

#endif
    
    // overwrite output vector with 0.0s and set correct length
    rvec1.Resize(size_col_, 0);
    for ( UInt k = 0; k < size_row_; ++k)
      for ( UInt kk = 0; kk < size_col_; ++kk)
        rvec1[kk] += data_[k][kk]*mvec1[k];
  }
  
  template<>
  void Matrix<int>::MultT(const SingleVector &mvec, SingleVector &rvec) const
  {
    EXCEPTION("undefined")
  }
  template<>
  void Matrix<unsigned int>::MultT(const SingleVector &mvec, SingleVector &rvec) const
  {
    EXCEPTION("undefined")
  }
  
  // **************
  //   operator<<
  // **************
  template<class S>
  std::ostream &operator<< ( std::ostream &out, const Matrix<S> &mat ) {


    // Set output format
    out.setf( std::ios::scientific );
    const UInt oldPrec = out.precision();
    const UInt newPrec = 4;
    out.precision( newPrec );

    const unsigned int mrows(mat.GetNumRows());
    const unsigned int mcols(mat.GetNumCols());
    
    for ( UInt i = 0; i < mrows; ++i ) {
      for ( UInt j = 0; j < mcols; ++j ) {
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
    transposedMat.Resize(size_col_, size_row_);
 
    UInt i,j;
 
    for (i = 0; i < size_col_; i++)
      for (j = 0; j < size_row_; j++)
        transposedMat.data_[i][j] = data_[j] [i];
  }


  template<class TYPE>
  void Matrix<TYPE>::DirectSolve( SingleVector & x1, const SingleVector & b1 ) const
  {

    Vector<TYPE> & x = dynamic_cast<Vector<TYPE>& >(x1);
    const Vector<TYPE> & b = dynamic_cast<const Vector<TYPE>& >(b1);
    
    Integer nmat = size_row_-1;
    Integer i, j, k, k1;
    
    //  the Gauss elimination 
    for(k = 0; k < nmat; ++k)
    {
      k1 = k + 1;
      for(i = k1; i <= nmat; ++i)
      {
        if (data_[k][k] != 0.0)
        {
          data_[i][k] /= data_[k][k];
          for (j=k1; j<=nmat; ++j)
            data_[i][j] -= data_[i][k] * data_[k][j];
        }
        else
        {
          //std::cerr<<"\n Matrix::DirectSolve Step " << k <<std::endl;
         //std::cerr<<"\n The element at position ("<< k << "," << k << ") of the matrix is zero. "<<std::endl;
          //                std::exit(0);
        }
      }
    }

    // solve Ly = b by forward substitution 

   
    Vector<TYPE> y(b.GetSize());

    for (i=0; i<=nmat; ++i)
      { 
        y[i] = b[i];
        for(j = 0; j < i; ++j)
          y[i] -= data_[i][j] * y[j];
      }

    // solve Ux = y backward substitution
    
    for(i = nmat; i >= 0; --i)
    {
      x[i] = y[i];
      for(j = nmat; j > i; --j)
        x[i] -= data_[i][j] * x[j];
      x[i] /= data_[i][i];
    }
  }
  
  template<>
  void Matrix<int>::DirectSolve( SingleVector & x1, const SingleVector & b1 ) const
  {
    EXCEPTION("undefined");
  }
  template<>
  void Matrix<unsigned int>::DirectSolve( SingleVector & x1, const SingleVector & b1 ) const
  {
    EXCEPTION("undefined");
  }


#ifdef USE_LAPACK
  // Compile OLAS and CFS++ with USE_LAPACK
  template<>
  void Matrix<Complex>::solveWithLapack(Matrix<Complex> & b1,
                                        lapackSysMatType & LAPACK_MATRIX_TYPE)
  {

    if(size_row_!=size_col_)
      EXCEPTION("Matrix is not quadratic, no solver!" );

    char lp_matType;
    lp_matType='L';
    
    Integer lp_nrRHS, lp_info, lp_lwork,lp_dim;
    Integer lp_lda, lp_ldb;
    const unsigned int bcols(b1.GetNumCols());
    lp_nrRHS=bcols;
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
    for (UInt i=0;i<size_row_;++i)
      for (UInt j=0;j<size_row_;++j){
        lp_sysVec[i+j*size_row_]=data_[i][j];
      }
    
    Vector<Complex> lp_rhsVec;
    const unsigned int brows(b1.GetNumRows());
    
    lp_rhsVec.Resize(brows*bcols);
 
    for (UInt i=0;i<brows;++i)
      for (UInt j=0;j<bcols;++j){
        lp_rhsVec[i+j*brows]=b1[i][j];
      }
    
    lp_rhsVecf77 = new F77complex16[size_row_*lp_nrRHS];
    lp_interchanges = new int[size_row_*size_row_];
    lp_workf77 = new F77complex16[size_row_*size_row_];
    lp_sysVecf77 = new F77complex16[size_row_*size_row_];
   

    // Convert CFS++ Vector<Complex> to Vector<F77complex16>
    for ( UInt count = 0; count < size_row_*bcols; ++count ) {
      CC2F77( lp_rhsVec[count], auxVal2 );
      lp_rhsVecf77[count] = auxVal2;
      }
    
    for (UInt count = 0; count < size_row_*size_row_; ++count ) {
      CC2F77( lp_sysVec[count], auxVal2 );
      lp_sysVecf77[count] = auxVal2;
    }
    
    for (UInt count=0, ss = lp_work.GetSize(); count < ss; ++count ) {
      CC2F77(lp_work[count], auxVal2);
      lp_workf77[count] = auxVal2;
    }
    
    switch (LAPACK_MATRIX_TYPE){
      
    case ZGESV:
      // solves systems with general system matrix
      zgesv_(&lp_dim , &lp_nrRHS, lp_sysVecf77, &lp_lda, 
             lp_interchanges, lp_rhsVecf77, &lp_ldb, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZGESV reports invalid input parameter" );
      }    
      break;
    case ZSYSV:
      
      lp_lwork=192;
      // solves systems with symmetric system matrix
      zsysv_(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77, 
             &lp_lda, lp_interchanges, lp_rhsVecf77, &lp_ldb,
             lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZSYSV reports invalid input parameter" );
      }
      break;
    case ZHESV:
      lp_lwork=192;
      // solves systems with hermitian system matrix
      zhesv_(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77,
             &lp_lda, lp_interchanges,lp_rhsVecf77, &lp_ldb, 
             lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZHESV reports invalid input parameter" );
      }
      break;

    default:
      std::cout<<LAPACK_MATRIX_TYPE<<std::endl;
      EXCEPTION( " the chosen routine is not yet implemented in solveWithLapack" );


    } // matches switch ()...
      
          
     //reconvert Fortran77 -> CFS ++ datatypes
    for ( UInt count = 0; count < size_row_*bcols; ++count ) 
      F772CC( lp_rhsVecf77[count], lp_rhsVec[count] );
          
    for ( UInt count = 0, ss = size_row_*size_row_; count < ss; ++count ) 
      F772CC( lp_sysVecf77[count], lp_sysVec[count]);
          
    for (UInt count=0, ss = lp_work.GetSize(); count < ss; ++count )
      F772CC(lp_workf77[count], lp_work[count]);
      
    // Writes result into b1
    for (UInt i=0;i<size_row_;++i)
      for (UInt j=0;j<bcols;++j)
        b1[i][j]=lp_rhsVec[i+j*brows];
  
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
  void Matrix<TYPE>::DyadicMult(const SingleVector & v1, const SingleVector & v2)
  {
    Vector<TYPE> const & vec1 = dynamic_cast<const Vector<TYPE>& >(v1);
    Vector<TYPE> const & vec2 = dynamic_cast<const Vector<TYPE>& >(v2);
  
    const UInt row = vec1.GetSize();
    const UInt col = vec2.GetSize();
  
    Resize(row, col);

    for(UInt actRow=0; actRow<row; actRow++)
      for(UInt actCol=0; actCol < col; actCol++)
        data_[actRow][actCol] = vec1[actRow] * vec2[actCol];
  }

  template<class TYPE>
  void Matrix<TYPE>::Invert (Matrix <TYPE> & inv) const
  {   

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION( "Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_col_ ) 
      EXCEPTION( "No quadratic matrix!" );
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
        // see Stoecker: "Taschenbuch Mathematischer Formeln und Moderner Verfahren" p.418
        inv.Resize(3,3);

        // === Old, recursive version ===
        // see Str: "Taschenbuch Mathematischer Formeln und Moderner Verfahren" p.418
        //for(UInt i=0; i<3; i++)
        //  for(UInt j=0; j<3; j++)
        //    inv[j][i] = Adjunct(i,j);      

        // === New, explicit version (from Wikipedia) ===
        inv[0][0] = data_[1][1] * data_[2][2] - data_[1][2] * data_[2][1];
        inv[0][1] = data_[0][2] * data_[2][1] - data_[0][1] * data_[2][2];
        inv[0][2] = data_[0][1] * data_[1][2] - data_[0][2] * data_[1][1];
        inv[1][0] = data_[1][2] * data_[2][0] - data_[1][0] * data_[2][2];
        inv[1][1] = data_[0][0] * data_[2][2] - data_[0][2] * data_[2][0];
        inv[1][2] = data_[0][2] * data_[1][0] - data_[0][0] * data_[1][2];
        inv[2][0] = data_[1][0] * data_[2][1] - data_[1][1] * data_[2][0];
        inv[2][1] = data_[0][1] * data_[2][0] - data_[0][0] * data_[2][1];
        inv[2][2] = data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];

        this->Determinant(det);
        inv *= 1/det;      
        break;
      
      default: 
        Double eps = 1e-40;
        TYPE pivot;
        TYPE pinv;

        //just copy the matrix
        inv.Resize(size_row_);
        for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
          inv.data_[0][k] = data_[0][k];

        //compute the inverse
        for ( UInt k=0; k<size_row_; k++) {
          pivot = inv[k][k];
          pinv  = 1/pivot;
          // std::abs not possible for uint
          if ( ( pivot >  0 ? pivot : -pivot) > eps ) {
            for (UInt j=0; j<size_row_; j++) {
              if ( j != k ) {
                inv[k][j] = -inv[k][j] * pinv;
                for ( UInt i=0; i<size_row_; i++ ) {
                  if ( i != k ) {
                    inv[i][j] = inv[i][j] + inv[i][k] * inv[k][j];
                  }
                }
              }
            }
            for ( UInt i=0;  i<size_row_; i++ ) {
              inv[i][k] = inv[i][k] * pinv;
            }
            inv[k][k] = pinv;
            //	    std::cout << "inv current:\n" << inv << std::endl;
          }
          else {
            EXCEPTION("Get divison by zero in matrix inversion" );
          }
        }
      }  
  }

  template<> void Matrix<Complex>::Invert (Matrix <Complex> & inv) const
  {
    EXCEPTION("Matrix<Complex>::Invert: Not implemented!" );
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::Adjunct (UInt i, UInt j) const
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX
    if (size_row_ != size_col_ ) 
      EXCEPTION("No quadratic matrix!" );
#endif
     
    assert(size_row_ == 3);
  
    StdVector<Integer> iVec(2);
    StdVector<Integer> jVec(2);
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
  Matrix<Double> Matrix<TYPE>::GetPart(  Global::ComplexPart part ) const {
    EXCEPTION( "Matrix::GetPart: Only Implemented for Real and Complex matrices!" );
    Matrix<Double> temp;
    return temp;
  }

  
  template<>
  Matrix<Double> Matrix<Double>::GetPart(  Global::ComplexPart part ) const {
    
    if ( part != Global::REAL ) {
      EXCEPTION("Matrix<Double>::GetPart: Only possible for REAL part." );
    }
      return *this;
    
  }

  template<>
  Matrix<Double> Matrix<Complex>::GetPart(Global::ComplexPart part) const
  {  
    Matrix<Double> ret(size_row_, size_col_);
    switch(part)
    {
    case Global::REAL:
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          ret[iRow][iCol]  = data_[iRow][iCol].real();
        }
      }
      break;
    case Global::IMAG:
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          ret[iRow][iCol]  = data_[iRow][iCol].imag();
        }
      }
      break;
    default:
      EXCEPTION("Matrix<Complex>::GetPart: Only possible for REAL or IMAG part!" );
    }
    
    return ret;
  }
  

  template<class TYPE>
  void Matrix<TYPE>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix ) {
    EXCEPTION( "Matrix::SetPart: Only Implemented for Real and Complex matrices!" );
  }
  
  template<>
  void Matrix<Double>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix )
  {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());
    assert((part == Global::REAL));

    *this = partMatrix;
  }

  template<>
  void Matrix<Complex>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix )
  {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());
    
    switch(part)
    {
    case Global::REAL:
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          data_[iRow][iCol]  = Complex( partMatrix[iRow][iCol],
              data_[iRow][iCol].imag() );
        }
      }
      break;
    case Global::IMAG:
      for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
        for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
          data_[iRow][iCol]  = Complex( data_[iRow][iCol].real(),
              partMatrix[iRow][iCol] );
        }
      }
      break;
    default:
      EXCEPTION( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!" );
    }
  }

  template<class TYPE>
  void Matrix<TYPE>::GetRow(Vector<TYPE>& vec, UInt row) const
  {
    assert(size_row_ > row);
    vec.Resize(size_col_);

    for(unsigned int i = 0; i < size_col_; i++)
      vec[i] = (*this)[row][i]; // do it faster if you like
  }

  template<class TYPE>
  void Matrix<TYPE>::GetCol(Vector<TYPE>& vec, UInt col) const
  {
    assert(size_col_ > col);
    vec.Resize(size_row_);

    for(unsigned int i = 0; i < size_row_; i++)
      vec[i] = (*this)[i][col]; // do it faster if you like
  }


  // copies a submatrix at the position (row, col) into subMat, 
  // the amount of copied elements depends on the size of subMat
  template<class TYPE>
  void Matrix<TYPE>::GetSubMatrix(Matrix<TYPE>& subMat, UInt startRow, UInt startCol) const
  {

#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      EXCEPTION( "undefined matrix" );
#endif

#ifdef CHECK_INDEX
    if (((subMat.size_row_ + startRow) > size_row_) || ((subMat.size_col_ + startCol) > size_col_) )
      EXCEPTION( "Submatrx to be written is to large! " );
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
#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      EXCEPTION( "undefined matrix" );
#endif
  
#ifdef CHECK_INDEX
    if ((subMat.size_row_ + startRow > size_row_) || (subMat.size_col_ + startCol > size_col_) )
      EXCEPTION("Submatrix to be read is to large!" );
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
#ifdef CHECK_INITIALIZED
    if (subMat.size_row_ == 0 || subMat.size_col_ == 0 || size_col_ == 0 || size_row_ == 0 ) 
      EXCEPTION("undefined matrix" );
#endif
  
#ifdef CHECK_INDEX
    if ((subMat.size_row_ + startRow > size_row_) || (subMat.size_col_ + startCol > size_col_) )
      EXCEPTION("Submatrix to be read is to large!");
#endif
  
    for( UInt actRow=0; actRow < subMat.size_row_; actRow++)
      for( UInt actCol=0; actCol < subMat.size_col_; actCol++)
        data_[actRow + startRow][actCol + startCol] += subMat[actRow][actCol];
  }



  /// converts a matrix into a vector, by appending successively all rows
  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_AppendRows(SingleVector & v) const
  {
    Vector<TYPE> & vec = dynamic_cast<Vector<TYPE>&>(v);

    vec.Resize(size_row_ * size_col_);
  
    for( UInt i=0; i < size_row_; i++)
      for( UInt j=0; j < size_col_; j++)
        vec[i*(size_col_) + j] = data_[i][j];
  }

  /// converts a matrix into a vector, by appending successively all colums
  template<class TYPE>
  void Matrix<TYPE>::ConvertToVec_AppendCols(SingleVector &v) const
  {

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
#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX 
    if (size_row_ != size_col_ ) 
      EXCEPTION("No square- matrix!" );
#endif

    for (UInt i = 0; i < size_row_; ++i)
      data_[i][i] *= factor;
  }


  /// gets the diagonal elements of a  matrix in a one column matrix
  template<class TYPE>
  void Matrix<TYPE>::GetDiagInMatrix(Matrix<TYPE>& columnMat) const
  {

#ifdef CHECK_INITIALIZED
    if (size_row_ == 0 || size_col_ == 0) 
      EXCEPTION("Undefined Matrix!" );
#endif

#ifdef CHECK_INDEX 
    if (size_row_ != size_col_ ) 
      EXCEPTION( "No square- matrix!" );
#endif

    columnMat.Resize(size_row_, 1);
    columnMat.Init();
  
    UInt i;
    for (i = 0; i < size_row_; i++)
      columnMat [i][0] = data_[i][i];

  }

  template<class TYPE>
  bool Matrix<TYPE>::IsSymmetric() const
  {
    for(UInt i = 1; i < size_row_; ++i)
      for(UInt j = i+1; j < size_col_; ++j)
        if(data_[i][j] != data_[j][i]) return false;
    
    return true;
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::NormL2() const
  {
    TYPE result(0);
    
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += data_[0][k] * data_[0][k];

    return static_cast<TYPE>(std::sqrt(result)); // for compilers
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::DiffNormL2(const Matrix<TYPE>& other) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ != other.size_row_ || size_col_ != other.size_col_)
      EXCEPTION("Incompatible matrices");
#endif

    TYPE result(0);
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
    {
      TYPE tmp = data_[0][k] - other.data_[0][k];
      tmp *= tmp;
      result += tmp;
    }

    return static_cast<TYPE>(std::sqrt(result)); // for compilers
  }

  template<class TYPE>
  TYPE Matrix<TYPE>::DiffNormL1(const Matrix<TYPE>& other) const
  {
#ifdef CHECK_INITIALIZED
    if (size_row_ != other.size_row_ || size_col_ != other.size_col_)
      EXCEPTION("Incompatible matrices");
#endif

    TYPE result(0);
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      result += static_cast<TYPE>(Abs(data_[0][k] - other.data_[0][k]));

    return result;
  }

  template<class TYPE>
  bool Matrix<TYPE>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if(std::isnan(data_[0][k])) return true;

    return false;
  }

  template<>
  bool Matrix<Complex>::ContainsNaN() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
    {
      if(std::isnan(data_[0][k].real())) return true;
      if(std::isnan(data_[0][k].imag())) return true;
    }
    return false;
  }


  template<class TYPE>
  bool Matrix<TYPE>::ContainsInf() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
      if(std::isinf(data_[0][k])) return true;

    return false;
  }

  template<>
  bool Matrix<Complex>::ContainsInf() const
  {
    for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
    {
      if(std::isinf(data_[0][k].real())) return true;
      if(std::isinf(data_[0][k].imag())) return true;
    }
    return false;
  }


  // Alternate version of symmetry checker, that will report
  // asymmetries to standard output

  // template<class TYPE> bool Matrix<TYPE>::IsSymmetric() {
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
#endif


// explicit template instantiation for SGI compiler
#ifdef __sgi
#pragma instantiate Matrix<UInt>
#pragma instantiate Matrix<Integer>
#pragma instantiate Matrix<Double>
#pragma instantiate Matrix<Complex>
#endif

} // end of namespace

