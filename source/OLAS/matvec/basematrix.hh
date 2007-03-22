// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEMATRIX_HH
#define OLAS_BASEMATRIX_HH


#include <iostream>

#include "matvec/typedefs.hh"
#include "utils/utils.hh"
#include "matvec/basevector.hh"

namespace OLAS {

  //! Generic matrix class

  //! All matrices in %OLAS are children of this
  //! class.
  class BaseMatrix {

  public:


    //! Default Constructor
    BaseMatrix(){
      ENTER_FCN("BaseMatrix::BaseMatrix");
    }

    //! Default Destructor

    //! Note that the default destructor must be virtual, so that if an
    //! inherited object is destroyed via a reference to BaseMatrix, the
    //! corresponding desctructor is called.
    virtual ~BaseMatrix(){
      ENTER_FCN("BaseMatrix::~BaseMatrix");
    }

    //! Return the block size of the matrix

    //! The method returns the block size of the matrix. By block size in
    //! this context we understand the number of degrees of freedom to which
    //! a matrix entry corresponds. A matrix using scalar entries e.g. has a
    //! block size of 1, while a matrix with 4x4 matrices as entries has a
    //! block size of 4.
    virtual Integer GetBlockSize() const = 0;

    //! Get the number of matrix rows

    //! The method returns the number of rows of the matrix. Note that this
    //! count is based upon the stored entry type. So e.g. for a 2x3 block
    //! matrix the returned value will be 2 independent of the size of the
    //! blocks.
    virtual Integer GetNrows() const = 0;

    //! Get the number of matrix columns

    //! The method returns the number of columns of the matrix. Note that this
    //! count is based upon the stored entry type. So e.g. for a 2x3 block
    //! matrix the returned value will be 3 independent of the size of the
    //! blocks.
    virtual Integer GetNcols() const = 0;

    //! Get the number of matrix rows on scalar level

    //! The method returns the number of rows of the matrix on the scalar
    //! level, i.e. the count is not based upon the stored entry type but
    //! accumulated over the scalar real/complex entries.
    //! This is different from the GetNrows() method.
    virtual Integer GetNrowsScalar() const = 0;

    //! Get the number of matrix columns on scalar level

    //! The method returns the number of columns of the matrix on the scalar
    //! level, i.e. the count is not based upon the stored entry type but
    //! accumulated over the scalar real/complex entries.
    //! This is different from the GetNcols() method.
    virtual Integer GetNcolsScalar() const = 0;

    //! Initialize matrix to zero
    
    //! This matrix sets all matrix entries to zero.
    virtual void Init() = 0;

    //! Add the multiple of a matrix to this matrix

    //! The method adds the multiple of a matrix to this matrix. In doing so
    //! the sparsity structure of the matrix mat is assumed to be identical to
    //! this matrix' structure.
    virtual void Add( const Double a, const BaseMatrix& mat ) = 0;
  
    //! Perform a matrix-vector multiplication rvec = this*mvec
    virtual void Mult(const BaseVector& mvec, BaseVector& rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec = transpose(this)*mvec
    virtual void MultT(const BaseVector& mvec, BaseVector& rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec += this*mvec
    virtual void MultAdd(const BaseVector& mvec, BaseVector& rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec += transpose(this)*mvec
    virtual void MultTAdd(const BaseVector& mvec, BaseVector& rvec) const = 0;

    //! Perform a matrix-vector multiplication rvec -= this*mvec
    virtual void MultSub(const BaseVector& mvec, BaseVector& rvec) const = 0;

    //! Scale all matrix entries by a constant factor

    //! This method allows to scale all non-zero matrix entries by a constant
    //! real-valued factor. It is currently pseudo-implemented in this class,
    //! since not all derived classes implement it, yet. The pseudo
    //! implementation will issue an error, if not over-written.
    virtual void Scale( Double factor ) {
      ENTER_FCN( "BaseMatrix::Scale" );
      (*error) << "BaseMatrix::Scale: Method must be implemented by derived "
               << "class, but is not!";
      Error( __FILE__, __LINE__ );
    };

    //! Compute the residual for a linear system with this matrix

    //! The method computes the residual \f$r=b-Ax\f$ where \f$A\f$ is the
    //! matrix represented by this matrix object.
    virtual void CompRes( BaseVector &r, const BaseVector &x,
			  const BaseVector& b ) const = 0;

    //! Determine maximum absolute value of diagonal entries

    //! This method determines the the maximal absolute value (on the scalar)
    //! level of the entries on the main diagonal of the matrix.
    virtual Double GetMaxDiag() const = 0;

    //!prints the matrix to *cla
    void Print(){Print(*cla);}
  
    //!prints the matrix to ostream os
    virtual void Print(std::ostream &os) const = 0;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file according to the
    //! MatrixMarket specifications. For details of the specification see
    //! http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    //! \param comment string to be inserted into file header (optional)
    virtual void Export( const Char *fname, const Char *comment = NULL )
      const = 0;


    // **** Matrix Data and Structure information ****
   
     
    //! Return the structural type of the matrix

    //! The method returns the structural type of the matrix. This is encoded
    //! as a value of the enumeration data type StructureType.
    virtual MatrixStructureType GetStructureType() const = 0;

    //! Return the Entry type of the matrix

    //! The method returns the entry type of the matrix on the scalar level.
    //! This is encoded as a value of the enumeration data type
    //! MatrixEntryType.
    virtual MatrixEntryType GetEntryType() const = 0;

    //! Is this a parallel matrix object?

    //! This method can be used to query, whether a matrix object is a
    //! sequential or a parallel object. By default the method will return
    //! false, since by default we work with sequential objects.
    //! Note that this method will exist, even if %OLAS was compiled
    //! without parallel support.
    virtual bool IsParallel() const {
      return false;
    }

  };

} // namespace

#endif // OLAS_BASEMATRIX_HH
