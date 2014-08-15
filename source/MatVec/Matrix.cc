#include "Matrix.hh"
#include "Vector.hh"
#include "opdefs.hh"

#include <string>
#include <cmath>
#include <def_build_type_options.hh>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_same.hpp>

#include "Utils/boost-serialization.hh"
#include "Utils/tools.hh"

#include "BLASLAPACKInterface.hh"

using namespace std;

using boost::tokenizer;
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
    std::string offset(n_offset, ' ');

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

    os << std::scientific << std::setprecision(7);
    
    switch(level)
    {
      case -1:
        for(UInt j = 0; j < size_row_; j++) {
          for(UInt i = 0; i < size_col_; i++)
            os << data_[j][i] << (j == size_row_ - 1 && i == size_col_ - 1 ? "" : " "); // space not for last element

          if(newline) os << std::endl;
          else os << " ";
        }
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
          {
            if(boost::is_complex<TYPE>::value)
            {
              Complex cval = (Complex) data_[j][i];
              os << cval.real() << "+" << cval.imag() << "i";
            }
            else
              os << data_[j][i];
            os << (i < size_col_-1 ? " " : "");
          }
          if(j < size_row_-1)
            os << "; ";
        }
        os << "]";
        break;
        //      std::cout << "dMat: \n " << dMat << std::endl;
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

  template<>
  void Matrix<Integer>::PerformRotation( const Matrix<Double>& R,  Matrix<Integer>& retMat ) const {
    EXCEPTION("Rotation only defined for double- and complex valued matrixes");
  }

  template<>
  void Matrix<UInt>::PerformRotation( const Matrix<Double>& R,  Matrix<UInt>& retMat ) const {
    EXCEPTION("Rotation only defined for double- and complex valued matrixes");
  }
  
  
  template<class TYPE>
  void Matrix<TYPE>::PerformRotation( const Matrix<Double>& R,  Matrix<TYPE>& retMat ) const {
    
    // Note; Currently the rotation only works for 3x3, 3x6, 6x3 and 6x6 matrices.
    // However, we should generalize the rotation also for 2x2, 2x4 and 4x4 matrices for the
    // 2D and axi case for mechanics.

    // get memory for transposed rotation matrix
    Matrix<Double> RT;
    RT.Resize(3,3);
    R.Transpose(RT);

    //get dimension of matrix
    UInt rowSize = size_row_;
    UInt colSize = size_col_;

    Matrix<TYPE> helpMat;

    if ( rowSize == 3 && colSize == 3) {
      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
      helpMat   = (*this) * RT;
      retMat = R * helpMat;
    }
    else if( (rowSize == 3 && colSize == 6) ||
             (rowSize == 6 && rowSize == 6 ) ) {
      // we also need Q;
      Matrix<Double> Q;

      // Composed Rotation Matrix
      // Ref.: M.Richter, "Entwicklung mechanischer Modelle zur analytischen
      // Beschreibung der Materialeigenschaften von textilbewehrtem Feinbeton",
      // Diss., Dresden, 2005, p. 27

      Q.Resize(6,6);  

      Q[0][0] = R[0][0]*R[0][0];
      Q[0][1] = R[0][1]*R[0][1];
      Q[0][2] = R[0][2]*R[0][2];
      Q[0][3] = 2.0*R[0][1]*R[0][2];
      Q[0][4] = 2.0*R[0][0]*R[0][2];
      Q[0][5] = 2.0*R[0][0]*R[0][1];

      Q[1][0] = R[1][0]*R[1][0];
      Q[1][1] = R[1][1]*R[1][1];
      Q[1][2] = R[1][2]*R[1][2];
      Q[1][3] = 2.0*R[1][1]*R[1][2];
      Q[1][4] = 2.0*R[1][0]*R[1][2];
      Q[1][5] = 2.0*R[1][0]*R[1][1];

      Q[2][0] = R[2][0]*R[2][0];
      Q[2][1] = R[2][1]*R[2][1];
      Q[2][2] = R[2][2]*R[2][2];
      Q[2][3] = 2.0*R[2][1]*R[2][2];
      Q[2][4] = 2.0*R[2][0]*R[2][2];
      Q[2][5] = 2.0*R[2][0]*R[2][1];

      Q[3][0] = R[1][0]*R[2][0];
      Q[3][1] = R[1][1]*R[2][1];
      Q[3][2] = R[1][2]*R[2][2];
      Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
      Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
      Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];

      Q[4][0] = R[0][0]*R[2][0];
      Q[4][1] = R[0][1]*R[2][1];
      Q[4][2] = R[0][2]*R[2][2];
      Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
      Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
      Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];

      Q[5][0] = R[0][0]*R[1][0];
      Q[5][1] = R[0][1]*R[1][1];
      Q[5][2] = R[0][2]*R[1][2];
      Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
      Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
      Q[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];


      //  std::cout << "R:\n" << R << std::endl;
      //  std::cout << "Q:\n" << Q << std::endl;
      //  std::cout << "Tensor orig:\n" << matTensor << std::endl;

      Matrix<Double> QT;
      QT.Resize(6,6);
      Q.Transpose(QT);

      if ( rowSize == 3 && colSize == 6 ) {
        helpMat   = (*this) * QT;
        retMat = R * helpMat;
      }
      else if (rowSize == 6 && colSize == 6 ) {
        helpMat   = (*this) * QT;
        retMat = Q * helpMat;
      }
      //  else {
      //    EXCEPTION("Cannot rotate tensor due to dimensions!");
      //  }
    } else {
      EXCEPTION("Tensor rotation currently only works for 3D matrices!");
    }
    

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


  /** Assigns a multiple of another matrix */
  template<class TYPE>
  void Matrix<TYPE>::Assign(const Matrix<TYPE>& other_mat, TYPE factor, bool size_tolerant)
  {
    if(size_row_ == other_mat.size_row_ && size_col_ == other_mat.size_col_)
      for(UInt k = 0, s = size_row_ * size_col_; k < s; ++k)
        data_[0][k] = factor * other_mat.data_[0][k];
    else
    {
      if(!size_tolerant)
        EXCEPTION("matrices do not match");
      // assure the non assigned area is zero
      if(size_row_ > other_mat.size_row_ || size_col_ > other_mat.size_col_)
        this->Init();
      for(unsigned int r = 0, nr = std::min(size_row_, other_mat.size_row_); r < nr; r++)
        for(unsigned int c = 0, nc = std::min(size_col_, other_mat.size_col_); c < nc; c++)
          data_[r][c] = factor * other_mat.data_[r][c];
    }
  }
  
  
  // perform matrix-matrix multiplication using BLAS (general case)
  template<class TYPE>
  void Matrix<TYPE>::Mult_Blas( const Matrix<TYPE>& mMat1,
                                Matrix<TYPE>& rMat1,
                                bool trans_a, bool trans_b,
                                TYPE alpha, TYPE beta, bool conjugate ) const {

#ifdef USE_BLAS
#ifdef CHECK_INDEX
    if((trans_a == true) && (trans_b == true)){
      if (size_row_ != mMat1.GetNumCols())
        EXCEPTION("incompatible dimension");
    } else if((trans_a == false) && (trans_b == true)){
      if (size_col_ != mMat1.GetNumCols())
        EXCEPTION("incompatible dimension");
    } else if((trans_a == true) && (trans_b == false)){
      if (size_row_ != mMat1.GetNumRows())
        EXCEPTION("incompatible dimension");
    } else {
      if (size_col_ != mMat1.GetNumRows())
        EXCEPTION("incompatible dimension");
    }
#endif

#ifdef CHECK_INITIALIZED
    if(size_row_ == 0 || size_col_ == 0)
      EXCEPTION("undefined Matrix");
    if(mMat1.GetNumRows() == 0 || mMat1.GetNumCols() ==0)
      EXCEPTION("undefined Matrix");
    if(rMat1.GetNumRows() == 0||rMat1.GetNumCols()==0)
      EXCEPTION("undefined Matrix");
#endif

    // because of the column-wise access of fortran the routine would calc
    // C^T = A^T * B^T if you give it C, A and B
    // but because A^T * B^T = (B * A)^T you can just calculate:
    // C = B * A
    // --> swap A and B

    // actually many transposed in complex are meant to be conjugate complex, e.g. B'(DB). Unfortunately
    // this is often not considered in CFS, so take care.
    char transp = conjugate && boost::is_complex<TYPE>::value ? 'c' : 't';

    char transa = trans_a ? transp : 'n';
    char transb = trans_b ? transp : 'n';

    // m would normally be the number of rows of C but Fortran accesses it as C^T so m is the number of columns
    //  in the same way n is now the number of rows
    //  k is the number of rows of op(B) and in our case op(B) = A or A^T so k would be rows of A in the first case and cols of A in the second one
    // but here you also have to remember the column wise access so cols and rows are swapped like for c

    int n = rMat1.GetNumRows();
    int m = rMat1.GetNumCols();
    int k = trans_a ? size_row_ : size_col_;

    TYPE* A = data_[0];
    TYPE* B = mMat1.data_[0];
    TYPE* C = rMat1.data_[0];

    // here you use the same properties as in the documentation of dgemm with the difference,
    // that our lda is LDB and ldb is LDA because we swapped A and B

    int lda = trans_a ? n : k;
    int ldb = trans_b ? k : m;
    int ldc = m;
    CallGEMM(&transb,&transa,&m,&n,&k,&alpha,B,&ldb,A,&lda,&beta,C,&ldc);
#else
    EXCEPTION("Compile with USE_BLAS = yes ");
#endif
   }
  

  template<class TYPE >
  void Matrix<TYPE>::CallGEMM(char* transb, char* transa, int* m, int* n, int* k, TYPE* alpha, TYPE* B, int* ldb, TYPE* A, int* lda, TYPE* beta, TYPE* C, int* ldc) const {
    EXCEPTION("wrong type");
  }

  template< >
  void Matrix<Double>::CallGEMM(char* transb, char* transa, int* m, int* n, int* k, Double* alpha, Double* B, int* ldb, Double* A, int* lda, Double* beta, Double* C, int* ldc) const {
    dgemm(transb,transa,m,n,k,alpha,B,ldb,A,lda,beta,C,ldc);
  }

  template< >
  void Matrix<Complex>::CallGEMM(char* transb, char* transa, int* m, int* n, int* k, Complex* alpha, Complex* B, int* ldb, Complex* A, int* lda, Complex* beta, Complex* C, int* ldc) const {
    zgemm(transb,transa,m,n,k,alpha,B,ldb,A,lda,beta,C,ldc);
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
  TYPE Matrix<TYPE>::ScalarProduct(const Matrix<TYPE>& other_mat) const
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
  
  // **************
  //   operator<<
  // **************
  template<class S>
  std::ostream &operator<< ( std::ostream &out, const Matrix<S> &mat ) {


    // Set output format
    out.setf( std::ios::scientific );
    const UInt oldPrec = out.precision();
    const UInt newPrec = 7;
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
  bool Matrix<TYPE>::IsHermitian(double eps) const
  {
    for(UInt r = 0; r < size_row_; r++)
      for(UInt c = r+1; c < size_col_; c++)
        if(!close(data_[r][c], data_[c][r], eps))
          return false;

    return true;
  }


  template< >
  bool Matrix<Complex>::IsHermitian(double eps) const
  {
    for(UInt r = 0; r < size_row_; r++)
    {
      if(!IsNoise(data_[r][r].imag()))
        return false;

      for(UInt c = r+1; c < size_col_; c++)
      {
        Complex u = data_[r][c]; // upper
        Complex l = data_[c][r]; // lower

        if(!close(u.real(), l.real()))
          return false;

        if(!close(std::abs(u.imag()), std::abs(l.imag())))
          return false;
      }
    }
    return true;
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
    
    std::complex<double> *lp_rhsVecf77;
    std::complex<double> *lp_sysVecf77;
    std::complex<double> *lp_workf77;
    std::complex<double> auxVal2;
    
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
    
    lp_rhsVecf77 = new std::complex<double>[size_row_*lp_nrRHS];
    lp_interchanges = new int[size_row_*size_row_];
    lp_workf77 = new std::complex<double>[size_row_*size_row_];
    lp_sysVecf77 = new std::complex<double>[size_row_*size_row_];
   

    // Convert CFS++ Vector<Complex> to Vector<std::complex<double>>
    for ( UInt count = 0; count < size_row_*bcols; ++count ) {
      lp_rhsVecf77[count] = lp_rhsVec[count];
      }
    
    for (UInt count = 0; count < size_row_*size_row_; ++count ) {
      lp_sysVecf77[count] = lp_sysVec[count];
    }
    
    for (UInt count=0, ss = lp_work.GetSize(); count < ss; ++count ) {
      lp_workf77[count] = lp_work[count];
    }
    
    switch (LAPACK_MATRIX_TYPE){
      
    case ZGESV:
      // solves systems with general system matrix
      zgesv(&lp_dim , &lp_nrRHS, lp_sysVecf77, &lp_lda, 
            lp_interchanges, lp_rhsVecf77, &lp_ldb, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZGESV reports invalid input parameter" );
      }    
      break;
    case ZSYSV:
      
      lp_lwork=192;
      // solves systems with symmetric system matrix
      zsysv(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77, 
            &lp_lda, lp_interchanges, lp_rhsVecf77, &lp_ldb,
            lp_workf77, &lp_lwork, &lp_info);

      if ( lp_info != 0 ) {
        EXCEPTION( "ZSYSV reports invalid input parameter" );
      }
      break;
    case ZHESV:
      lp_lwork=192;
      // solves systems with hermitian system matrix
      zhesv(&lp_matType, &lp_dim , &lp_nrRHS, lp_sysVecf77,
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
      lp_rhsVec[count] = lp_rhsVecf77[count];
          
    for ( UInt count = 0, ss = size_row_*size_row_; count < ss; ++count ) 
      lp_sysVec[count] = lp_sysVecf77[count];
          
    for (UInt count=0, ss = lp_work.GetSize(); count < ss; ++count )
      lp_work[count] = lp_workf77[count];
      
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
    std::complex<double> auxValC;

    std::complex<double> * lp_af77 = new std::complex<double>[size_row_*size_row_];
    double * lp_wf77 = new double[size_row_];
    std::complex<double> * lp_workf77 = new std::complex<double>[99];
    double * lp_rworkf77 = new double[3*size_row_-2];
      
    // Convert CFS++ Vector<Complex> to Vector<std::complex<double>>
    for ( UInt count = 0; count < size_row_; count++ ) 
      for ( UInt countC = 0; countC <size_row_; countC++ ) {
        lp_af77[countC*size_row_+count] = data_[count][countC];
      }
    
    for ( UInt count = 0; count < size_row_; count++ ) {
      lp_wf77[count] = lp_w[count];
    }
    
    for (UInt count=0; count < lp_rwork.GetSize();count++){
      lp_rworkf77[count] = lp_rwork[count];
    }
    
    zheev( &lp_jobz, &lp_uplo, &lp_N, lp_af77, &lp_lda, lp_wf77, 
           lp_workf77, &lp_lworkf77, lp_rworkf77 ,&lp_infof77); 
    
    // reconvert f772C++
    for (UInt count=0; count < lp_work.GetSize();count++)
      lp_work[count] = lp_workf77[count];
    
    for ( UInt count = 0; count < size_row_; count++ ) 
      lp_w[count] = lp_wf77[count];
    
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
    void Matrix<TYPE>::Invert_Lapack() {
    EXCEPTION("General case not implemented");
  }
  
  template<> void Matrix<Double>::Invert_Lapack() {
#ifdef CHECK_INDEX
    if( size_row_ != size_col_) {
      EXCEPTION("Can only invert square matrices");
    }
#endif

#ifndef USE_LAPACK
    EXCEPTIO("Compile with LAPACK support for matrix inversion");
#else
    
    
    int *ipiv = new int[size_row_];
    int n = size_row_;
    int lwork = size_row_ * size_row_;
    double *work = new double[lwork];
    int info;

    // calculate LU-factorization of block
    dgetrf(&n,&n,data_[0],&n,ipiv,&info);
    if( info != 0 ) {
      EXCEPTION("Error during LU-factorization of matrix. "
                << "Error value is " << info );
    }
    // invert matrix using previous LU factorization
    dgetri(&n,data_[0],&n,ipiv,work,&lwork,&info);
    if( info != 0 ) {
      EXCEPTION("Error during inversion of matrix. "
                << "Error value is " << info );
    }

    delete[] ipiv;
    delete[] work;
#endif
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
  void Matrix<TYPE>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                              bool zeroOtherPart ) {
    EXCEPTION( "Matrix::SetPart: Only Implemented for Real and Complex matrices!" );
  }
  
  template<>
  void Matrix<Double>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                                bool zeroOtherPart)
  {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());
    assert((part == Global::REAL));

    *this = partMatrix;
  }

  template<>
  void Matrix<Complex>::SetPart( Global::ComplexPart part, const Matrix<Double> & partMatrix,
                                 bool zeroOtherPart)
                                 {
    assert(size_col_ == partMatrix.GetNumCols());
    assert(size_row_ == partMatrix.GetNumRows());

    if( zeroOtherPart) {
      // ------------------
      //  Zero other part
      // ------------------
      switch(part)
      {
        case Global::REAL:
          for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
            for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
              data_[iRow][iCol]  = Complex( partMatrix[iRow][iCol], 0.0);
            }
          }
          break;
        case Global::IMAG:
          for ( UInt iRow = 0; iRow < size_row_; iRow++ ) {
            for ( UInt iCol = 0; iCol < size_col_; iCol++ ) {
              data_[iRow][iCol]  = Complex( 0.0, partMatrix[iRow][iCol] );
            }
          }
          break;
        default:
          EXCEPTION( "Matrix<Complex>::SetPart: Only possible for REAL or IMAG part!" );
      }
    } else {


      // ------------------
      //  Keep other part
      // ------------------
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
  void Matrix<TYPE>::SetSubMatrixByInd( const Matrix<TYPE>& subMat,
                                        const StdVector<UInt>& rowInd,
                                        const StdVector<UInt>& colInd ) {

    UInt nrows = rowInd.GetSize();
    UInt ncols = colInd.GetSize();
    for( UInt iRow = 0; iRow < nrows; iRow++ ) {
      for( UInt iCol = 0; iCol < ncols; iCol++ ) {
#ifdef CHECK_INDEX
      if( rowInd[iRow] > size_row_ || 
          colInd[iCol] > size_col_ ) {
        EXCEPTION("Index pair (" << rowInd[iRow] << ", " << colInd[iCol]
                  << ") larger than matrix size being " << size_row_
                  << " x " << size_col_);
      }
#endif
      data_[rowInd[iRow]][colInd[iCol]]= subMat[iRow][iCol];
      }
    }
  }
  
  template<class TYPE>void Matrix<TYPE>::
  GetSubMatrixByInd( Matrix<TYPE>&ret, 
                     const StdVector<UInt>& rowInd,
                     const StdVector<UInt>& colInd ) const {

    UInt nrows = rowInd.GetSize();
    UInt ncols = colInd.GetSize();
    ret.Resize(nrows,ncols);
    for( UInt iRow = 0; iRow < nrows; iRow++ ) {
      for( UInt iCol = 0; iCol < ncols; iCol++ ) {
#ifdef CHECK_INDEX
      if( rowInd[iRow] > size_row_ || 
          colInd[iCol] > size_col_ ) {
        EXCEPTION("Index pair (" << rowInd[iRow] << ", " << colInd[iCol]
                  << ") larger than matrix size being " << size_row_
                  << " x " << size_col_);
      }
#endif
        ret[iRow][iCol] = data_[rowInd[iRow]][colInd[iCol]];
      }
    }
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
  bool Matrix<TYPE>::IsSymmetric(double eps) const
  {
    for(UInt i = 1; i < size_row_; ++i)
      for(UInt j = i+1; j < size_col_; ++j)
        if(std::abs(data_[i][j] != data_[j][i]) > eps) return false;

    return true;
  }

  template<>
  bool Matrix<Complex>::IsSymmetric(double eps) const
  {
    for(UInt i = 1; i < size_row_; ++i)
      for(UInt j = i+1; j < size_col_; ++j)
        if(!close(data_[i][j], data_[j][i], eps)) return false;

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

  template<class TYPE>
  Matrix<TYPE> Conj( const Matrix<TYPE>& m) {
    return m;
  }

  template<>
  Matrix<Complex> Conj( const Matrix<Complex>& m) {
    UInt numRows = m.GetNumRows();
    UInt numCols = m.GetNumCols();
    Matrix<Complex> ret(numRows, numCols);
    for( UInt i = 0; i < numRows; i++ ) {
      for (UInt j = 0; j < numCols; j++ ) {
        ret[i][j] = std::conj(m[i][j]);
      }
    }
    return ret;
  }

  template<class TYPE>
  Matrix<TYPE> Herm( const Matrix<TYPE>& m ) {
    return m;
  }

  template<>
  Matrix<Complex> Herm( const Matrix<Complex>& m ) {
    UInt numRows = m.GetNumRows();
    UInt numCols = m.GetNumCols();
    Matrix<Complex> herm(numCols, numRows);
    for( UInt i = 0; i < numCols; i++ ) {
      for (UInt j = 0; j < numRows; j++ ) {
        herm[i][j] = std::conj(m[j][i]);
      }
    }
    return herm;
  }

  template Matrix<Double> Conj( const Matrix<Double>& m);
  template Matrix<Double> Herm( const Matrix<Double>& m);

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

