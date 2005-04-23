#ifndef FILE_CFSMATRIX_2004
#define FILE_CFSMATRIX_2004
 
#include <iostream>

#include "Utils/cfsvector.hh"
#include "Utils/tools.hh"

namespace CoupledField{

// Forward declarations
class CFSVector;


//! Abstract interface class for a general dense matrix
/*! Abstract interface class for a general dense matrix
  \note Altough this class provides a general interface
   to the matrix class , one should always prefer a cast into
   the current type of matrix, e.g. 
   dynamic_cast<Matrix<Double>*>(ptCFSMatrix)
*/
class CFSMatrix{
public:

  //! Constructor 
  //! creates an empty matrix of size 0x0
  CFSMatrix(){
    ENTER_IFCN("CFSMatrix::CFSMatrix");
  };

  //! Constructor
  //! creates a matrix of size nRows x nCols, initalized
  //! with zeroes
  /*! 
    \param nRows (input) Number of rows
    \param nCols (input) Number of columns)
  */
  CFSMatrix(const Integer nRows, const Integer nCols){};
  
  //! Default Copy Construcctor
  CFSMatrix(const CFSMatrix &m){};

  //! Destructor
  virtual ~CFSMatrix(){ENTER_IFCN("CFSMatrix::CFSMatrix");};

  //! Hard coded query if values are complex
  virtual Boolean IsComplex() = 0;

  //! Initialize matrix with a given entry.
  //! If no entry given, it gets initalized with zeroes
  /*!
    \param val (input,opt.) Entry the matrix gets initalized with
  */
  //! \note This method does not change the size of the matrix
  virtual void Init(const Double input = 0.0)
  {Error("CFSMatrix::Init(): Not implemented here",__FILE__,__LINE__);}
  
  //! Change size of quadratic matrix
  /*!
    \param size (input) Number of rows / columns
  */
  //! \note the matrix contains afterwards only zeroes
  virtual void Resize(const Integer size) = 0;
  
  //! Change size of general matrix 
  /*!
    \param nRows (input) Number of rows
    \param nCols (input) Number of columns)
  */
  //! \note the matrix contains afterwards only zeroes
   virtual void Resize(const Integer nRows, const Integer nCols) = 0;
 
  //! Get the number of rows
  virtual Integer GetSizeRow() const = 0;
  
  //! Get the number of columns
  virtual Integer GetSizeCol() const = 0;
  
  //! Set the entry 'val' at position (row,col) in the matrix
  /*!
    \param row (input) Row of entry
    \param col (input) Column of entry
    \param val (input) Value to be set
  */
  void SetEntry(const Integer row, const Integer col, const Double val)
  {Error("CFSMatrix::SetEntry(): Not implemented here",__FILE__,__LINE__);} 
  
  //! Add'val' to the matrix entry at position (row,col) in the matrix
  /*!
    \param row (input) Row of entry
    \param col (input) Column of entry
    \param val (input) Value to be added
  */
  virtual void AddToEntry(const Integer row, const Integer col, const Double val)
  {Error("CFSMatrix::AddToEntry(): Not implemented here",__FILE__,__LINE__);} 
  
  //! Get the entry 'val' at position (row,col) in the matrix
  /*!
    \param row (input) Row of entry
    \param col (input) Column of entry
    \param val (output) Variable the value is written into
  */ 
  virtual void GetEntry(const Integer row, const Integer col, Double & val) const
 {Error("CFSMatrix::GetEntry(): Not implemented here",__FILE__,__LINE__);} 
  
   //! Calculates the determinant (up to size 3)
  /*!
    \param val (output) Return value of the method
  */
  virtual void Determinant(Double & val) const
  {Error("CFSMatrix::Determinant(): Not implemented here",__FILE__,__LINE__);}   
 
  virtual void Determinant(Complex & val) const
  {Error("CFSMatrix::Determinant(): Not implemented here",__FILE__,__LINE__);} 
  
  //! Invert the matrix and store it in 'inv' (up to size 3)
  /*!
    \param inv (output) Inverse of the matrix
  */
  virtual void Invert(CFSMatrix & inv) const = 0;

  //! Transpose the matrix and store it in 'trans'
  /*!
    \param trans (output) Transposed matrix
  */
  //! \note The matrix itself gets not changed.
  //! \note If the transposed of a matrix is needed for a operation
  //! with a vector, the according function like 'MultT' should be used
  virtual void Transpose(CFSMatrix & trans) const = 0;

 
  //! Solves directly a small system of equations of the form Ax=b
  //! using LU-Decomposition (no pivoting - so dear numerical people, do not exhaust this method!)
  //! The Matrix A=LU contains afterwards the the values of L in the lower triangular,
  //! and the values of U in the upper part.
 virtual void DirectSolve(CFSVector & x, CFSVector & b) const
  {Error("CFSMatrix:DirectSolve(): Not implemented here",__FILE__,__LINE__);}

  //!
  virtual  void DyadicMult(CFSVector & vec1) const
  {Error("CFSMatrix:: DyadicMult(): Not implemented here",__FILE__,__LINE__);} 
  

  //!
  virtual void DyadicMult(CFSVector & vec1, CFSVector & vec2) const
  {Error("CFSMatrix:: DyadicMult: Not implemented here",__FILE__,__LINE__);} 
  
  
  //! copies a submatrix at the position (row, col) into subMat, 
  //! the amount of copied elements depends on the size of subMat
  virtual void GetSubMatrix(CFSMatrix & subMat, const Integer nRows, const Integer nCols) const = 0;

  //! overwrites the matrix elements at the position (row, col) with subMat
  //! in a rectangular (submatrix) way
  virtual void SetSubMatrix(CFSMatrix & subMat, const Integer nRows, const Integer nCols) = 0;

  //! converts a matrix into a vector, by appending successively all rows
  virtual void ConvertToVec_AppendRows(CFSVector & vec) const
  {Error("CFSMatrix::ConvertToVec_AppendRows(): Not implemented here",__FILE__,__LINE__);} 
   
  //! converts a matrix into a vector, by appending successively all columns
  virtual void ConvertToVec_AppendCols(CFSVector & vec) const
  {Error("CFSMatrix::ConvertToVec_AppendCols(): Not implemented here",__FILE__,__LINE__);} 
  
  //! scales the diagonal elements of a  matrix by a factor
  virtual void ScaleDiagElems(const Double factor)
  {Error("CFSMatrix::ScaleDiagElems(): Not implemented ere",__FILE__,__LINE__);} 
  
  virtual void ScaleDiagElems(const Complex factor)
  {Error("CFSMatrix::ScaleDiagElems(): Not implemented here",__FILE__,__LINE__);} 

  //!
  virtual void Add(const Double fac, const CFSMatrix & mat)
  {Error("CFSMatrix::Add(): Not implemented here",__FILE__,__LINE__);} 
  
  virtual void Add(const Complex fac, const CFSMatrix & mat)
  {Error("CFSMatrix::Add(): Not implemented here",__FILE__,__LINE__);} 
  
  //! Perform a matrix-vector multiplication rvec = this*mvec
  virtual void Mult(CFSVector & mvec, CFSVector & rvec) const
  {Error("CFSMatrix:Mult(): Not implemented here",__FILE__,__LINE__);}

  //! Perform a matrix-matrix multiplication rMat = this*mMat
  virtual void Mult(CFSMatrix & mMat, CFSMatrix & rMat) const
  {Error("CFSMatrix:Mult(): Not implemented here",__FILE__,__LINE__);}

  //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
  virtual void MultT(const CFSVector & mvec, CFSVector & rvec) const = 0;

  //! Perform a matrix-vector multiplication rvec += this*mvec
  virtual void MultAdd(const CFSVector & mvec, CFSVector & rvec) const = 0;
  
  //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
  virtual void MultTAdd(const CFSVector & mvec, CFSVector& rvec) const = 0;
  
  //! Perform a matrix-vector multiplication rvec -= this*mvec
  virtual void MultSub(const CFSVector & mvec, CFSVector & rvec) const = 0;


#define DECL_BASEMATRIX_FCN(TYPE)						\
  virtual void Init(const TYPE input = TYPE())					\
  {Error("CFSMatrix::Init(): Not implemented here",__FILE__,__LINE__);}		\
  virtual void SetEntry(const Integer row, const Integer col, const TYPE val)	\
  {Error("CFSMatrix::SetEntry(): Not implemented here",__FILE__,__LINE__);}	\
  virtual void AddToEntry(const Integer row, const Integer col, const TYPE val)	\
  {Error("CFSMatrix::AddToEntry(): Not implemented here",__FILE__,__LINE__);}	\
  virtual void GetEntry(const Integer row, const Integer col, TYPE & val) const	\
  {Error("CFSMatrix::GetEntry(): Not implemented here",__FILE__,__LINE__);}
  
DECL_BASEMATRIX_FCN(Integer)
DECL_BASEMATRIX_FCN(Complex)
   
 


};


} // namespace
#endif
