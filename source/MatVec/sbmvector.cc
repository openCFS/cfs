// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include "MatVec/sbmvector.hh"
#include "generatematvec.hh"


namespace CoupledField {


  // *****************************************************
  //   Constructor for empty vectors of specified length
  // *****************************************************
  SBM_Vector::SBM_Vector( UInt size ) {
    size_ = size;
    myEntryType_ = BaseMatrix::NOENTRYTYPE;
    NEWARRAY( subVec_, SingleVector*, size );
    for ( UInt i = 1; i <= size_; i++ ) {
      subVec_[i] = NULL;
    }
  }


  // *******************
  //   Deep destructor
  // *******************
  SBM_Vector::~SBM_Vector() {
    if ( subVec_ != NULL ) {
      for ( UInt i = 1; i <= size_; i++ ) {
	delete subVec_[i];
      }
      DELETEARRAY( subVec_ );
    }
  }


  // ***********************************************************
  //   Set number of vector entries and re-size internal array
  // ***********************************************************
  void SBM_Vector::SetSize( UInt size ) {
    if ( subVec_ != NULL ) {
      for ( UInt i = 1; i <= size_; i++ ) {
	delete subVec_[i];
      }
      DELETEARRAY( subVec_ );
    }
    size_ = (UInt)size;
    NEWARRAY( subVec_, SingleVector*, size );
    for ( UInt i = 1; i <= size_; i++ ) {
      subVec_[i] = NULL;
    }
  }

  // **********************
  //   Assignment operator
  // **********************
  SBM_Vector& SBM_Vector::operator= ( const SBM_Vector &bVec ) {

    // first, delete sub-vectors
    if ( subVec_ != NULL ) {
      for ( UInt i = 1; i <= size_; i++ ) {
        delete subVec_[i];
      }
      DELETEARRAY( subVec_ );
    }
    size_ = 0;

    // then, copy sub-vectors from original vector class
    size_ = bVec.size_;
    NEWARRAY( subVec_, SingleVector*, size_ );
    for ( UInt i =  1; i <= size_; i++ ) {
      // determine entrytype
      subVec_[i] = CopySingleVectorObject( *bVec.subVec_[i]);
    }

    return *this;
  }



  // *********
  //   Inner
  // *********

  // Double version
  void SBM_Vector::Inner( const BaseVector& vec, Double &retval ) const {


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
      EXCEPTION( WRONG_CAST_MSG );
    }

  }


  // Complex version
  void SBM_Vector::Inner( const BaseVector& vec, Complex &retval ) const {


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
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // ***************
  //   Print vector
  // ***************
  std::string SBM_Vector::ToString( char separator ) const {
    std::stringstream os;
    for( UInt i = 1; i <= size_; i++ ) {
      if( subVec_[i] != NULL ) {
        os <<   "sub-Vector #" << i
           << "\n--------------\n";
        os <<  subVec_[i]->ToString( separator );
      }
    }
    return os.str();
  }

  // ************************************
  //   Add vector to this vector object
  // ************************************
  void SBM_Vector::Add( const BaseVector &vec ) {

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
            EXCEPTION( "Is there already a copy constructor for "
                     << "StdMatrix objects?" );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // ******************************
  //   Add SingleVector to ith entry
  // ******************************
  void SBM_Vector::AddToSubVector( const SingleVector &vec, const UInt i ) {


    if ( subVec_[i] != NULL ) {
      subVec_[i]->Add( vec );
    }
    else {
      // A copy constructor is needed here
      EXCEPTION( "Is there already a copy constructor for "
               << "StdMatrix objects?" );
    }
  }


  // *************************************************
  //   Scale this vector object and add other vector
  // *************************************************

  // Double Version
  void SBM_Vector::Axpy( const Double alpha, const BaseVector &y ) {

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
            EXCEPTION( "Is there already a copy constructor for "
                     << "StdMatrix objects?" );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // Complex Version
  void SBM_Vector::Axpy( const Complex alpha, const BaseVector &y ) {

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
            EXCEPTION( "Is there already a copy constructor for "
                     << "StdMatrix objects?" );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // **********************************************************
  //   Add a scaled version of a vector to this vector object
  // **********************************************************

  // Double Version
  void SBM_Vector::Add( Double alpha, const BaseVector& v ) {


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
            EXCEPTION( "Is there already a copy constructor for "
                     << "StdMatrix objects?" );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // Complex Version
  void SBM_Vector::Add( Complex alpha, const BaseVector& v ) {


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
            EXCEPTION( "Is there already a copy constructor for "
                     << "StdMatrix objects?" );
          }
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // ***************************************************************
  //   Replace this vector object by the sum of two scaled vectors
  // ***************************************************************
  void SBM_Vector::Add( Double alpha, const BaseVector& y,
			Double beta, const BaseVector& z ) {

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
          EXCEPTION( "Gosh, this will be tricky to implement, or not?" );
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }

  void SBM_Vector::Add( Complex alpha, const BaseVector& y,
			Complex beta, const BaseVector& z ) {

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
          EXCEPTION( "Gosh, this will be tricky to implement, or not?" );
        }
      }
    }

    // Treat downcast failure
    catch(std::bad_alloc &e) {
      EXCEPTION( WRONG_CAST_MSG );
    }
  }


  // ************************************************
  //   Compute Euclidean norm of this vector object
  // ************************************************
  Double SBM_Vector::NormL2() const {
    Double norm = 0.0;
    Double auxval;
    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        auxval = subVec_[i]->NormL2();
        norm += auxval*auxval;
      }
    }
    return sqrt( norm );
  }


  // *****************
  //   Export vector
  // *****************
  void SBM_Vector::Export( const char *fname ) const {


    std::stringstream fileName;
    std::string outFile;

    for ( UInt j = 1; j <= size_; j++ ) {

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


    for ( UInt i = 1; i <= size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarMult( factor );
      }
    }
  }

}
