// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_BASEMATRIX_HH
#define OLAS_BASEMATRIX_HH


#include <iostream>
<<<<<<< HEAD

#include "General/Enum.hh"
=======
>>>>>>> origin/master
#include "General/defs.hh"
#include "General/Exception.hh"


namespace CoupledField {

  //! Forward class declaration
  class BaseVector;

  //! Generic matrix class

  //! All matrices in %OLAS are children of this
  //! class.
  class BaseMatrix {

  public:

    // ---- Enumerations for dynamic Matrix information ----

    //! Structural type of the matrix.

    //! This enumeration data type describes the structural type of the matrix.
    //! This is used to distinguish e.g. SBM_Matrices from StdMatrices and so
    //! on. It can take one of the following values:
    //! - NOSTRUCTURETYPE
    //! - SBM_MATRIX
    //! - SPARSE_MATRIX
    typedef enum { NOSTRUCTURETYPE, SBM_MATRIX, SPARSE_MATRIX, DENSE_MATRIX } StructureType;
    // see Enum<BaseMatrix::StructureType> MatrixStructureTypeEnum in Environment.hh

    //! Entry type of the matrix.

    //! This enumeration data type describes the entry type of the matrix.
    //! The entry type describes the scalar type of the entries e.g. Integer,
    //! Double or Complex. It can take one of the following values:
    //! - NOENTRYTYPE
    //! - INTEGER
    //! - DOUBLE
    //! - COMPLEX
    //! - UINT
    //! - F77REAL8
    //! - F77COMPLEX16
    typedef enum { NOENTRYTYPE, INTEGER, DOUBLE, COMPLEX,	UINT, F77REAL8, F77COMPLEX16 } EntryType;
    // see MatrixEntryTypeEnum

    //! Storage type of the matrix.

    //! This enumeration data type describes the storage format of a sparse
    //! matrix. It can take one of the following values
    //! - NOSTORAGETYPE
    //! - SPARSE_SYM
    //! - SPARSE_NONSYM
    //! - LAPACK_GBMATRIX
    typedef enum { NOSTORAGETYPE, SPARSE_SYM, SPARSE_NONSYM,
       LAPACK_GBMATRIX, DIAG, VAR_BLOCK_ROW} StorageType;
    // see MatrixStorageTypeEnum   

    //! Available output formats for matrices.

    //! This enumeration data type lists the supported  output formats of a
    //! sparse matrix. It can take one of the following values
    //! - MATRIX_MARKET
    //! - HARWELL_BOEING
    typedef enum { MATRIX_MARKET = 1, HARWELL_BOEING = 2, PLAIN} OutputFormat;
    // see MatrixOutputFormatEnum
    
    //! Default Constructor
    BaseMatrix(){
    }

    //! Default Destructor

    //! Note that the default destructor must be virtual, so that if an
    //! inherited object is destroyed via a reference to BaseMatrix, the
    //! corresponding desctructor is called.
    virtual ~BaseMatrix(){
    }

    //! Get the number of matrix rows

    //! The method returns the number of rows of the matrix. Note that this
    //! count is based upon the stored entry type. So e.g. for a 2x3 block
    //! matrix the returned value will be 2 independent of the size of the
    //! blocks.
    virtual UInt GetNumRows() const = 0;

    //! Get the number of matrix columns

    //! The method returns the number of columns of the matrix. Note that this
    //! count is based upon the stored entry type. So e.g. for a 2x3 block
    //! matrix the returned value will be 3 independent of the size of the
    //! blocks.
    virtual UInt GetNumCols() const = 0;

    //! Initialize matrix to zero

    //! This matrix sets all matrix entries to zero.
    virtual void Init() = 0;
    
    //! Return matrix as separated string
    virtual std::string ToString( char colSeparator = ' ',
                                  char rowSeparator = '\n' ) const = 0;

    //! Export the matrix to a file in MatrixMarket format

    //! The method will export the matrix to an ascii file.
    //! \param fname name of output file
    //! \param format matrix output format 1 = MatrixMarket, 2 = Harwell-Boeing (optional)
    //! \param comment string to be inserted into file header (optional)
    virtual void Export( const std::string& fname,
                         OutputFormat format = MATRIX_MARKET,
                         const std::string& comment = "") const = 0;


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
      EXCEPTION("BaseMatrix::Scale: Method must be implemented by derived " \
                "class, but is not!");
    };

    virtual void Scale( Complex factor ) {
      EXCEPTION("BaseMatrix::Scale: Method must be implemented by derived " \
                "class, but is not!");
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

    // **** Matrix Data and Structure information ****


    //! Return the structural type of the matrix

    //! The method returns the structural type of the matrix. This is encoded
    //! as a value of the enumeration data type StructureType.
    virtual BaseMatrix::StructureType GetStructureType() const = 0;

    //! The method returns the storage type of the matrix. This is encoded
    //! as a value of the enumeration data type StorageType.
    virtual BaseMatrix::StorageType GetStorageType() const = 0;

    //! Return the Entry type of the matrix

    //! The method returns the entry type of the matrix on the scalar level.
    //! This is encoded as a value of the enumeration data type
    //! MatrixEntryType.
    virtual BaseMatrix::EntryType GetEntryType() const = 0;

    //! Return memory requirement in Bytes
    
    //! This method returns the memory size of the matrix in Bytes
    virtual Double GetMemoryUsage() const = 0;

    /** gives some debug information */
    std::string ToInfoString() const;

  };

  // Function for determining matrix/vector entry type (i.e. Integer, Double,
  // Complex) for enum-type refer to Environment.hh

  //! Class for determining the type of a matrix/vector entry on scalar level

  //! This class is used to associate with a template that specifies a tiny
  //! matrix or vector the type of the matrix' or vector's entry on the
  //! scalar level. The class contains a single public attribute called
  //! <em>M_EntryType</em> of type MatrixEntryType for storing this
  //! information.
  template<class T>
  class EntryType{
  public:
    static const BaseMatrix::EntryType M_EntryType = BaseMatrix::NOENTRYTYPE;
  };

  template<> class EntryType<double> { public: static const BaseMatrix::EntryType M_EntryType = BaseMatrix::DOUBLE; };
  template<> class EntryType<Complex> { public: static const BaseMatrix::EntryType M_EntryType = BaseMatrix::COMPLEX; };
  template<> class EntryType<int> { public: static const BaseMatrix::EntryType M_EntryType = BaseMatrix::INTEGER; };
  template<> class EntryType<unsigned int> { public: static const BaseMatrix::EntryType M_EntryType = BaseMatrix::UINT; };

} // namespace

#endif // OLAS_BASEMATRIX_HH
