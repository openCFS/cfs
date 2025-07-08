#ifndef OLAS_SPARSITYPATTERNS
#define OLAS_SPARSITYPATTERNS

namespace CoupledField {

  //! Base class for administration of sparsity patterns

  //! This class forms the base class for the hierarchy of classes that
  //! are used in the administration of matrix sparsity patterns. Its sole
  //! purpose is to provide a common data type to be used by the PatternPool
  //! class.
  class BaseSparsityPattern {
  public:
    virtual ~BaseSparsityPattern() {};
  };


  //! Auxilliary class for administrating patterns of SCRS_Matrix objects
  //! in combination with the PatternPool class concept.
  class SCRS_Pattern : public BaseSparsityPattern {

  public:

    //! Default constructor
    SCRS_Pattern() {
      cidx_ = NULL;
      rptr_ = NULL;
      diagPtr_ = NULL;
    }

    //! Destructor
    virtual ~SCRS_Pattern() {
      delete [] ( cidx_ );
      delete [] ( rptr_ );
      delete [] ( diagPtr_);
    }

    //! Array with column indices of non-zero matrix entries
    UInt *cidx_;

    //! Array with starting indices of rows in cidx_ array
    UInt *rptr_;

    //! Array containing the indices of the diagonal matrix entries
    UInt *diagPtr_;
  };
  
  //! Auxilliary class for administrating patterns of CRS_Matrix objects
  //! in combination with the PatternPool class concept.
  class CRS_Pattern : public BaseSparsityPattern {

  public:

    //! Default constructor
    CRS_Pattern() {
      cidx_   = NULL;
      rptr_   = NULL;
      diagPtr_ = NULL;
    }

    //! Destructor
    virtual ~CRS_Pattern() {
     delete [] ( cidx_ );
     delete [] ( rptr_ );
     delete [] ( diagPtr_ );
    }

    //! Array with column indices of non-zero matrix entries
    UInt *cidx_;

    //! Array with starting indices of rows in cidx_ array
    UInt *rptr_;
    
    //! Array containing the indices of the diagonal matrix entries
    UInt *diagPtr_;
  };

}

#endif
