// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//! \file generatematvec.hh
//! This is the header file for the several factories that can be used to
//! generate most matrix/vector objects. The exception being parallel
//! matrices. Note that <b>GenerateVectorObject()</b> will fall back on
//! <b>GenerateParVectorObject()</b>, if the input matrix is parallel.
//! \remark Note that the documentation of the functions in this file
//! will not appear on this webpage, but in the description of the namespace
//! OLAS to which the functions belong.


#ifndef OLAS_GENERATEMATVEC_HH
#define OLAS_GENERATEMATVEC_HH


#include "MatVec/basematrix.hh"
#include "MatVec/basevector.hh"
#include "MatVec/stdmatrix.hh"
#include "MatVec/sbmmatrix.hh"


namespace CoupledField {


  // Forward declarations
  class SBM_Matrix;

  //! Creates the associated vector to a given matrix

  //! Creates the associated vector to a given matrix, e.g.
  //! an RVector if an RMatrix is given. The result
  //! is returned as a BaseVector* and might have to be
  //! casted in to SingleVector or SBMVector afterwards.
  //! The vector is initialised to contain zeros on the scalar level.
  BaseVector* GenerateVectorObject(const BaseMatrix &m);

  //! Generation of a vector of specified type

  //! This function can be used to dynamically generate a vector object.
  //! The type and size of the generated vector object can be specified
  //! via the input parameters.
  //! \note
  //! - Since there is no VectorType enumeration data type in %OLAS
  //!   we use a combination of MatrixStorageType, MatrixEntryType and
  //!   an integer specifying the blockSize of the vector entries to
  //!   fix the type of the vector.
  //! - As the name implies SBM_Vectors, but also parallel vector objects,
  //!   cannot be generated via this method.
  //! - The vector is initialised so that the entries on the scalar level are
  //!   all identical to zero.
  BaseVector* GenerateSingleVectorObject( const BaseMatrix::StorageType sType,
                                          const BaseMatrix::EntryType eType,
                                          const UInt length );

  //! Craeate a vector by copy

  //! This functions will generate a copy of a sparse olas vector and
  //! return the pointer to the copy.
  SingleVector* CopySingleVectorObject( const SingleVector& origVec );


  //! Creates a StdMatrix object

  //! This routine will generate an StdMatrix object and return a pointer to
  //! it. Note that the matrix object is not ready to use before the sparsity
  //! pattern has not been set using SetSparsityPattern.
  //! \param etype the type of scalar matrix entries (Double or Complex)
  //! \param stype the storage layout of the matrix
  //! \param dof the number of degrees of freedom; this determines the
  //!            (block)type of the entries.
  //! \param nrows number of matrix rows
  //! \param ncols number of matrix columns
  //! \param fill number of non-zero matrix entries
  StdMatrix* GenerateStdMatrixObject( const BaseMatrix::EntryType etype,
                                      const BaseMatrix::StorageType stype,
                                      const Integer nrows,
                                      const Integer ncols,
                                      const Integer fill );

  //! Creates a StdMatrix object

  //! This routine will generate an StdMatrix object by copying an existing
  //! StdMatrix and return a pointer to the new matrix.
  StdMatrix* CopyStdMatrixObject( const StdMatrix &origMat );

  //! Creates a SBM_Matrix object

  //! This routine will generate an %SBM_Matrix object and return a pointer
  //! to it. The method will populate the %SBM_Matrix with sub-matrices and
  //! also set their sparsity pattern. This information and also the
  //! dimension of the sub-matrices is taken from the
  //! graphs administered by the specified graph manager object. Thus, the
  //! %SBM_Matrix returned is ready to use.
  //! \param eType        the type of scalar matrix entries
  //!                     (Double or Complex)
  //! \param rowDim       number of rows of the %SBM_Matrix
  //! \param colDim       number of columns of the %SBM_Matrix
  //! \param graphManager pointer to a GraphManager object for SBM_Matrices
  //! \param symmetric signals whether the %SBM_Matrix is symmetric or not
  //! \note
  //! - The block-size of the sub-matrices will always be one, since the
  //!   %SBM_Matrix class currently only allows sub-matrices with scalar
  //!   entries.
  //! - If symmetric is set to 'false' we select SPARSE_NONSYM as storage
  //!   type for all sub-matrices. If it is set to 'true' we select
  //!   SPARSE_SYM for the sub-matrices on the diagonal.
  SBM_Matrix* GenerateSBMMatrixObject( BaseMatrix::EntryType eType,
                                       UInt rowDim, UInt colDim,
                                       GraphManager *graphManager,
                                       bool symmetric );
}

#endif
