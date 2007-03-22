// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_SPARSITYPATTERNS
#define OLAS_SPARSITYPATTERNS

namespace OLAS {

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
      ENTER_FCN( "SCRS_Pattern::SCRS_Pattern" );
      cidx_ = NULL;
      rptr_ = NULL;
    }

    //! Destructor
    virtual ~SCRS_Pattern() {
      ENTER_FCN( "SCRS_Pattern::~SCRS_Pattern" );
      DeleteArray( cidx_ );
      DeleteArray( rptr_ );
    }

    //! Array with column indices of non-zero matrix entries
    Integer *cidx_;

    //! Array with starting indices of rows in cidx_ array
    Integer *rptr_;
  };

}

#endif
