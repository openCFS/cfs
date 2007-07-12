// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include "matvec/sbmvector.hh"


namespace OLAS {


  // *****************************************************
  //   Constructor for empty vectors of specified length
  // *****************************************************
  SBM_Vector::SBM_Vector( UInt size ) {
    ENTER_FCN( "SBM_Vector::SBM_Vector" );
    size_ = size;
    NewArray( subVec_, SparseVector*, size );
    for ( UInt i = 1; i <= size_; i++ ) {
      subVec_[i] = NULL;
    }
  }


  // *******************
  //   Deep destructor
  // *******************
  SBM_Vector::~SBM_Vector() {
    ENTER_FCN( "SBM_Vector::~SBM_Vector" );
    if ( subVec_ != NULL ) {
      for ( UInt i = 1; i <= size_; i++ ) {
	delete subVec_[i];
      }
      DeleteArray( subVec_ );
    }
  }


  // ***********************************************************
  //   Set number of vector entries and re-size internal array
  // ***********************************************************
  void SBM_Vector::SetSize( Integer size ) {
    ENTER_FCN( "SBM_Vector::SetSize" );
    if ( subVec_ != NULL ) {
      for ( UInt i = 1; i <= size_; i++ ) {
	delete subVec_[i];
      }
      DeleteArray( subVec_ );
    }
    size_ = (UInt)size;
    NewArray( subVec_, SparseVector*, size );
    for ( UInt i = 1; i <= size_; i++ ) {
      subVec_[i] = NULL;
    }
  }


  // *********
  //   Inner
  // *********

  // Double version
  void SBM_Vector::Inner( const BaseVector& vec, Double &retval ) const {

    ENTER_FCN( "SBM_Vector::Inner" );

    // Downcast BaseVector to SBM_Vector
    try {
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(vec);
      Double auxval;
      retval = 0.0;
      for ( UInt i = 1; i <= size_; i++ ) {
        if ( subVec_[i] != NULL && sbmvec.subVec_[i] != NULL ) {
          (*this)(i).Inner( sbmvec(i), auxval );
          retval += auxval;
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }

  }


  // Complex version
  void SBM_Vector::Inner( const BaseVector& vec, Complex &retval ) const {

    ENTER_FCN( "SBM_Vector::Inner" );

    // Downcast BaseVector to SBM_Vector
    try {
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(vec);
      Complex auxval;
      retval = 0.0;
      for ( UInt i = 1; i <= size_; i++ ) {
        if ( subVec_[i] != NULL && sbmvec.subVec_[i] != NULL ) {
          (*this)(i).Inner( sbmvec(i), auxval );
          retval += auxval;
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // ************************************
  //   Add vector to this vector object
  // ************************************
  void SBM_Vector::Add( const BaseVector &vec ) {
    ENTER_FCN( "SBM_Vector::Add" );

    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(vec);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL ) {
          if ( sbmvec.subVec_[i] != NULL ) {
            subVec_[i]->Add( sbmvec(i) );
          }
        }
        else {
          // A copy constructor is needed here
          if ( sbmvec.subVec_[i] != NULL ) {
            (*error) << "Is there already a copy constructor for "
                     << "StdMatrix objects?";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // ******************************
  //   Add SparseVector to ith entry
  // ******************************
  void SBM_Vector::AddToSubVector( const SparseVector &vec, const Integer i ) {

    ENTER_FCN( "SBM_Vector::AddToSubVector" );

    if ( subVec_[i] != NULL ) {
      subVec_[i]->Add( vec );
    }
    else {
      // A copy constructor is needed here
      (*error) << "Is there already a copy constructor for "
               << "StdMatrix objects?";
      Error( __FILE__, __LINE__ );
    }
  }


  // *************************************************
  //   Scale this vector object and add other vector
  // *************************************************

  // Double Version
  void SBM_Vector::Axpy( const Double alpha, const BaseVector &y ) {
    ENTER_FCN( "SBM_Vector::Axpy" );

    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(y);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL ) {
          if ( sbmvec.subVec_[i] != NULL ) {
            subVec_[i]->Axpy( alpha, sbmvec(i) );
          }
        }
        else {
          // A copy constructor is needed here
          if ( sbmvec.subVec_[i] != NULL ) {
            (*error) << "Is there already a copy constructor for "
                     << "StdMatrix objects?";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // Complex Version
  void SBM_Vector::Axpy( const Complex alpha, const BaseVector &y ) {
    ENTER_FCN( "SBM_Vector::Axpy" );

    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(y);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL ) {
          if ( sbmvec.subVec_[i] != NULL ) {
            subVec_[i]->Axpy( alpha, sbmvec(i) );
          }
        }
        else {
          // A copy constructor is needed here
          if ( sbmvec.subVec_[i] != NULL ) {
            (*error) << "Is there already a copy constructor for "
                     << "StdMatrix objects?";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // **********************************************************
  //   Add a scaled version of a vector to this vector object
  // **********************************************************

  // Double Version
  void SBM_Vector::Add( Double alpha, const BaseVector& v ) {

    ENTER_FCN( "SBM_Vector::Add" );
    
    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(v);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL ) {
          if ( sbmvec.subVec_[i] != NULL ) {
            subVec_[i]->Add( alpha, sbmvec(i) );
          }
        }
        else {
          // A copy constructor is needed here
          if ( sbmvec.subVec_[i] != NULL ) {
            (*error) << "Is there already a copy constructor for "
                     << "StdMatrix objects?";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // Complex Version
  void SBM_Vector::Add( Complex alpha, const BaseVector& v ) {

    ENTER_FCN( "SBM_Vector::Add" );
    
    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbmvec = dynamic_cast<const SBM_Vector&>(v);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL ) {
          if ( sbmvec.subVec_[i] != NULL ) {
            subVec_[i]->Add( alpha, sbmvec(i) );
          }
        }
        else {
          // A copy constructor is needed here
          if ( sbmvec.subVec_[i] != NULL ) {
            (*error) << "Is there already a copy constructor for "
                     << "StdMatrix objects?";
            Error( __FILE__, __LINE__ );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // ***************************************************************
  //   Replace this vector object by the sum of two scaled vectors
  // ***************************************************************
  void SBM_Vector::Add( Double alpha, const BaseVector& y,
			Double beta, const BaseVector& z ) {
    ENTER_FCN( "SBM_Vector::Add" );
    
    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbm_y = dynamic_cast<const SBM_Vector&>(y);
      const SBM_Vector& sbm_z = dynamic_cast<const SBM_Vector&>(z);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL && sbm_y.subVec_[i] != NULL &&
             sbm_z.subVec_[i] != NULL ) {
          subVec_[i]->Add( alpha, sbm_y(i), beta, sbm_z(i) );
        }
        else {
          (*error) << "Gosh, this will be tricky to implement, or not?";
          Error( __FILE__, __LINE__ );
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }

  void SBM_Vector::Add( Complex alpha, const BaseVector& y,
			Complex beta, const BaseVector& z ) {
    ENTER_FCN( "SBM_Vector::Add" );
    
    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbm_y = dynamic_cast<const SBM_Vector&>(y);
      const SBM_Vector& sbm_z = dynamic_cast<const SBM_Vector&>(z);

      // Let's hope both vectors have the same size
      for ( UInt i = 1; i <= size_; i++ ){
        if ( subVec_[i] != NULL && sbm_y.subVec_[i] != NULL &&
             sbm_z.subVec_[i] != NULL ) {
          subVec_[i]->Add( alpha, sbm_y(i), beta, sbm_z(i) );
        }
        else {
          (*error) << "Gosh, this will be tricky to implement, or not?";
          Error( __FILE__, __LINE__ );
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      Error( WRONG_CAST_MSG, __FILE__, __LINE__ );
    }
  }


  // ************************************************
  //   Compute Euclidean norm of this vector object
  // ************************************************
  Double SBM_Vector::NormEuclid() const {
    ENTER_FCN( "SBM_Vector::NormEuclid" );
    Double norm = 0.0;
    Double auxval;
    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        auxval = subVec_[i]->NormEuclid();
        norm += auxval*auxval;
      }
    }
    return sqrt( norm );
  }


  // *****************
  //   Export vector
  // *****************
  void SBM_Vector::Export( const Char *fname ) const {

    ENTER_FCN( "SBM_Vector::Export" );

    std::stringstream fileName;
    std::string outFile;

    for ( Integer j = 1; j <= size_; j++ ) {

      // construct file name
      fileName.str( "" );
      fileName << fname << '_' << j << ".vec";
      outFile = fileName.str();

      // export sub-vector
      if ( subVec_[j] != NULL ) {
        subVec_[j]->Export( outFile.c_str() );
      }
    }
  }


  // ********************
  //   ScalarDiv (real)
  // ********************
  void SBM_Vector::ScalarDiv( const Double factor ) {

    ENTER_IFCN( "SBM_Vector::ScalarDiv" );

    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarDiv( factor );
      }
    }
  }


  // *********************
  //   ScalarMult (real)
  // *********************
  void SBM_Vector::ScalarMult( const Double factor ) {

    ENTER_IFCN( "SBM_Vector::ScalarMult" );

    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarMult( factor );
      }
    }
  }


  // ***********************
  //   ScalarDiv (complex)
  // ***********************
  void SBM_Vector::ScalarDiv( const Complex factor ) {

    ENTER_IFCN( "SBM_Vector::ScalarDiv" );

    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarDiv( factor );
      }
    }
  }


  // ************************
  //   ScalarMult (complex)
  // ************************
  void SBM_Vector::ScalarMult( const Complex factor ) {

    ENTER_IFCN( "SBM_Vector::ScalarMult" );

    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarMult( factor );
      }
    }
  }

}
