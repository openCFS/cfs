// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DENSEMATRIX_HH
#define FILE_DENSEMATRIX_HH
 
#include <iostream>
#include "MatVec/SingleVector.hh"
#include "Utils/tools.hh"

namespace CoupledField{

  // Forward declarations
  class SingleVector;


  //! Abstract interface class for a general dense matrix
  /*! Abstract interface class for a general dense matrix
    \note Altough this class provides a general interface
    to the matrix class , one should always prefer a cast into
    the current type of matrix, e.g. 
    dynamic_cast<Matrix<Double>*>(ptDenseMatrix)
  */
  class DenseMatrix{
  public:

    //! Constructor 
    //! creates an empty matrix of size 0x0
    DenseMatrix(){
    };

    //! Destructor
    virtual ~DenseMatrix() {};

    //! Constructor
    //! creates a matrix of size nRows x nCols, initalized
    //! with zeroes
    /*! 
      \param nRows (input) Number of rows
      \param nCols (input) Number of columns)
    */
    explicit DenseMatrix(const UInt nRows, const UInt nCols){};
  
    //! Default Copy Construcctor
    DenseMatrix(const DenseMatrix &m){};

    //! Destructor

    //! Get entry type of matrix
    virtual BaseMatrix::EntryType GetEntryType() = 0;

    //! Initialize matrix with zeroes
    //! \note This method does not change the size of the matrix
    virtual void Init() = 0;
  
    //! Change size of quadratic matrix
    /*!
      \param size (input) Number of rows / columns
    */
    virtual void Resize(const UInt size ) = 0;
  
    //! Change size of general matrix 
    /*!
      \param nRows (input) Number of rows
      \param nCols (input) Number of columns)
    */
    virtual void Resize(const UInt nRows, const UInt nCols ) = 0;
 
    //! Get the number of rows
    virtual UInt GetNumRows() const = 0;
  
    //! Get the number of columns
    virtual UInt GetNumCols() const = 0;
  
    //! Invert the matrix and store it in 'inv' (up to size 3)
    /*!
      \param inv (output) Inverse of the matrix
    */
    virtual void Invert(DenseMatrix & inv) const = 0;

    //! Transpose the matrix and store it in 'trans'
    /*!
      \param trans (output) Transposed matrix
    */
    //! \note The matrix itself gets not changed.
    //! \note If the transposed of a matrix is needed for a operation
    //! with a vector, the according function like 'MultT' should be used
    virtual void Transpose(DenseMatrix & trans) const = 0;
 
    //! Solves directly a small system of equations of the form Ax=b
    //! using LU-Decomposition (no pivoting - so dear numerical people, do not exhaust this method!)
    //! The Matrix A=LU contains afterwards the the values of L in the lower triangular,
    //! and the values of U in the upper part.
    virtual void DirectSolve(SingleVector & x, const SingleVector & b) const = 0;

    //!
    virtual  void DyadicMult(const SingleVector & vec1) = 0;

    //!
    virtual void DyadicMult(const SingleVector & vec1, const SingleVector & vec2) = 0;
  
  
    //! copies a submatrix at the position (row, col) into subMat, 
    //! the amount of copied elements depends on the size of subMat
    virtual void GetSubMatrix(DenseMatrix & subMat, const UInt nRows, const UInt nCols) const = 0;

    //! overwrites the matrix elements at the position (row, col) with subMat
    //! in a rectangular (submatrix) way
    virtual void SetSubMatrix(const DenseMatrix & subMat, const UInt nRows, const UInt nCols) = 0;

    //! converts a matrix into a vector, by appending successively all rows
    virtual void ConvertToVec_AppendRows(SingleVector & vec) const = 0;
   
    //! converts a matrix into a vector, by appending successively all columns
    virtual void ConvertToVec_AppendCols(SingleVector & vec) const = 0;
  
    //! Perform a matrix-vector multiplication rvec = this*mvec
    virtual void Mult(const SingleVector & mvec, SingleVector & rvec) const = 0;

    //! Perform a matrix-matrix multiplication rMat = this*mMat
    virtual void Mult(const DenseMatrix & mMat, DenseMatrix & rMat) const = 0;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    virtual void MultT(const SingleVector & mvec, SingleVector & rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec += this*mvec
    virtual void MultAdd(const SingleVector & mvec, SingleVector & rvec) const = 0;
  
    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    virtual void MultTAdd(const SingleVector & mvec, SingleVector& rvec) const = 0;
  
    //! Perform a matrix-vector multiplication rvec -= this*mvec
    virtual void MultSub(const SingleVector & mvec, SingleVector & rvec) const = 0;

    virtual std::string ToString(int level=0) = 0;

  private:    


  };


} // namespace
#endif
