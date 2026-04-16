#include <cmath>
#include "SBM_Vector.hh"
#include "generatematvec.hh"
#include "Utils/ToolsFull.hh"

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  SBM_Vector::SBM_Vector(  BaseMatrix::EntryType entryType ) {
    size_ = 0;
    myEntryType_ = entryType;
    ownSubVectors_ = true;
  }

  // *****************************************************
  //   Constructor for empty vectors of specified length
  // *****************************************************
  SBM_Vector::SBM_Vector( UInt size,
                          BaseMatrix::EntryType entryType ) {
    size_ = size;
    ownSubVectors_ = true;
    myEntryType_ = entryType;
    subVec_.Resize( size );
    
    // Create vectors, if entryType is set. Otherwise initialize with
    // NULL-pointer
    if( myEntryType_ != BaseMatrix::NOENTRYTYPE ) {
      for( UInt i = 0; i < size_; ++i ) {
        subVec_[i] = 
            GenerateSingleVectorObject(myEntryType_, 0 );
      }
    } else {
      subVec_.Init( NULL );
    }
  }


  // *************************
  //   Weak Copy Constructor
  // *************************
  SBM_Vector::SBM_Vector( const SBM_Vector& origVector, UInt numRows ) {
    
    size_ = numRows;
    myEntryType_ = origVector.myEntryType_;
    // This is a weak copy, so no ownership of the pointers
    ownSubVectors_ = false;
    subVec_.Resize( size_ );
    subVec_.Init( NULL );
    for ( UInt k = 0; k < subVec_.GetSize(); k++ ) {
      subVec_[k] = origVector.subVec_[k];
    }
  }
  
  
  // *******************
  //   Deep destructor
  // *******************
  SBM_Vector::~SBM_Vector() {
   
    // only delete if we have ownership
    if( ownSubVectors_ ) {
      for ( UInt i = 0; i < size_; i++ ) {
        delete subVec_[i];
      }
    }
  }

  // ******************************
  //   Reset the initial EntryType
  // ******************************
  void SBM_Vector::ResetEntryType(BaseMatrix::EntryType entryType){
    if(size_ == 0 && subVec_ == NULL && ownSubVectors_){
      myEntryType_ = entryType;
    }else{
      EXCEPTION("SBM_Vector::ResetEntryType Cannot reset the entry-type\n"
                "because the SBM_Vector was already modified");
    }
  }


  // ***********************************************************
  //   Set number of vector entries and re-size internal array
  // ***********************************************************
  void SBM_Vector::SetSize( UInt size ) {
    
    if( !ownSubVectors_ ) {
      EXCEPTION( "As this is a weak copy of a SBM-vector, I refuse to "
                 "change the size" );
    }
    
    for ( UInt i = 0; i < size_; i++ ) {
      delete subVec_[i];
    }
    
    size_ = (UInt)size;
    subVec_.Resize( size_ );
    subVec_.Init( NULL );
  }
  
  
  
  // **********************
  //   Perform, Resize
  // **********************
  void SBM_Vector::Resize( UInt newSize) {
    
    // If new and old size are the same, just leave
    if( newSize == size_)
      return;
    
    // check for weak copy of SBM vector
    if( !ownSubVectors_ ) {
      EXCEPTION( "As this is a weak copy of a SBM-vector, I refuse to "
                 "change the size" );
    }
    
    UInt keepSize = std::min(newSize, size_);
    
    // delete all vectors not needed anymore
    for( UInt i = keepSize; i < size_; ++i ) {
      delete subVec_[i];
    }
    
    subVec_.Resize( newSize );
    size_ = newSize;
    
    // Create new, empty sub-vectors, if the entryType is set
    for( UInt i = keepSize; i < size_; ++i ) {
      if( myEntryType_ != BaseMatrix::NOENTRYTYPE ) {
        subVec_[i] = 
            GenerateSingleVectorObject(myEntryType_, 0 );
      } else {
        subVec_[i] = NULL;
      }
    }
  }
  
  // **********************
  //   Assign Sub-vector
  // **********************
  void SBM_Vector::SetSubVector( SingleVector *subvec, UInt i ) {

    //Check if we are even allowed to set a new sub-vector
    if ( !ownSubVectors_ ) {
      EXCEPTION( "As this is a weak copy of a SBM-vector, I refuse to "
          "get overwritten" );
    }
    
    // Check index of the new vector
#ifdef CHECK_INDEX
    if( i >= size_ ) {
      EXCEPTION( "SetSubVector: Can not set sub-vector " << i
                 << " as size of the SBM-vector is " << size_ );
    }
#endif
    // Check if entryType of new vector is the same as of the SBM-vector 
    if( myEntryType_ != BaseMatrix::NOENTRYTYPE &&
        subvec->GetEntryType() != myEntryType_ ) {
      EXCEPTION( "Can not set sub-vector of type '"
          << MatrixEntryTypeEnum.ToString( subvec->GetEntryType() )
          << "', as the SBM-vector is of type '"
          << MatrixEntryTypeEnum.ToString( myEntryType_ ) << "'." )
    }
    
    delete subVec_[i];
    subVec_[i] = subvec;
    myEntryType_ = subvec->GetEntryType();
  }


  // **********************
  //   Assignment operator
  // **********************
  SBM_Vector& SBM_Vector::operator= ( const SBM_Vector &bVec ) {
    // check, if vectors have the same length
    // 1) same length: loop over sub-vectors and copy entries
    if( size_ == bVec.size_ ){
      // loop over all sub-vectors
      for( UInt i = 0; i < size_; ++i ) {
			
        // switch depending on entry type
        if( myEntryType_ == BaseMatrix::DOUBLE ) {
          Vector<Double> & myVec = dynamic_cast<Vector<Double>& >(*subVec_[i]);
          Vector<Double> & rVec = dynamic_cast<Vector<Double>& >(*bVec.subVec_[i]);
          myVec = rVec;
        } else if( myEntryType_ == BaseMatrix::COMPLEX ) {
          Vector<Complex> & myVec = dynamic_cast<Vector<Complex>& >(*subVec_[i]);
          Vector<Complex> & rVec = dynamic_cast<Vector<Complex>& >(*bVec.subVec_[i]);
          myVec = rVec;
        } else {
          EXCEPTION( "Can only copy double- and complex "
              << "valued sub-vectors");
        }
      }
    }
    
    // 2) different length: delete sub-vectors and copy all
    //    vectors from bVec
    else {
      // Only set subvectors, if this is not a cheap copy
      if( !ownSubVectors_ ) {
        EXCEPTION( "As this is a weak copy of a SBM-vector, I refuse to "
            "get overwritten" );
      }


      // first, delete sub-vectors
      for ( UInt i = 0; i < size_; i++ ) {
        delete subVec_[i];
      }

      size_ = 0;
      // then, copy sub-vectors from original vector class
      size_ = bVec.size_;
      subVec_.Resize( size_ );
      myEntryType_ = bVec.myEntryType_;
      for ( UInt i =  0; i < size_; i++ ) {
        subVec_[i] = CopySingleVectorObject( *bVec.subVec_[i]);
      }
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
      for ( UInt i = 0; i < size_; i++ ) {
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
      for ( UInt i = 0; i < size_; i++ ) {
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
  std::string SBM_Vector::ToString(ToStringFormat format, const std::string& sep, int digits) const
  {
    std::stringstream os;
    for( UInt i = 0; i < size_; i++ ) {
      if( subVec_[i] != NULL ) {
        os <<  "sub-Vector #" << i << "\n--------------\n";
        os <<  subVec_[i]->ToString(format, sep, digits);
        os << "\n";
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
      for ( UInt i = 0; i < size_; i++ ){
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
      for ( UInt i = 0; i < size_; i++ ){
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
      for ( UInt i = 0; i < size_; i++ ){
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
      for ( UInt i = 0; i < size_; i++ ){
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
      for ( UInt i = 0; i < size_; i++ ){
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
  void SBM_Vector::Add( Double alpha, const BaseVector& y, Double beta, const BaseVector& z ) {

    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbm_y = dynamic_cast<const SBM_Vector&>(y);
      const SBM_Vector& sbm_z = dynamic_cast<const SBM_Vector&>(z);

      // Let's hope both vectors have the same size
      for ( UInt i = 0; i < size_; i++ ){
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
      for ( UInt i = 0; i < size_; i++ ){
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
            Complex beta, const BaseVector& z , UInt h) {

    try {

      // The dynamic cast will only work if vec really is an SBM_Vector
      const SBM_Vector& sbm_y = dynamic_cast<const SBM_Vector&>(y);
      const SBM_Vector& sbm_z = dynamic_cast<const SBM_Vector&>(z);

      // Let's hope both vectors have the same size
      for ( UInt i = 0; i < size_; i++ ){
        if ( subVec_[i] != NULL && sbm_y.subVec_[i] != NULL &&
             sbm_z.subVec_[i] != NULL ) {
          if ( i == h){
            subVec_[i]->Add( alpha, sbm_y(i), beta, sbm_z(i) );
          }
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
    for ( UInt i = 0; i < size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        auxval = subVec_[i]->NormL2();
        norm += auxval*auxval;
      }
    }
    return sqrt( norm );
  }

  double SBM_Vector::NormL2(const SBM_Vector& other) const
  {
    assert(size_ == other.size_);
    double norm = 0.0;
    for(unsigned int i = 0; i < size_; i++)
      if(subVec_[i] != nullptr)
        norm += CoupledField::NormL2(subVec_[i], other.subVec_[i]); // need to give namespace to find from ToolsFull.hh?!

    return norm;
  }


  // *****************
  //   Export vector
  // *****************
  void SBM_Vector::Export(const std::string& fname, BaseMatrix::OutputFormat format) const {


    std::stringstream ss;
    std::string outFile;

    for ( UInt j = 0; j < size_; j++ ) {

      // construct file name
      if(size_ > 1)
        ss << '_' << j;
      outFile = fname + ss.str();
      ss.str(std::string());

      // export sub-vector
      if ( subVec_[j] != NULL ) {
        subVec_[j]->Export( outFile.c_str(), format );
      }
    }
  }

  // ********************
  //   ScalarDiv (real)
  // ********************
  void SBM_Vector::ScalarDiv( const Double factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarDiv( factor );
      }
    }
  }


  // *********************
  //   ScalarMult (real)
  // *********************
  void SBM_Vector::ScalarMult( const Double factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarMult( factor );
      }
    }
  }


  // ***********************
  //   ScalarDiv (complex)
  // ***********************
  void SBM_Vector::ScalarDiv( const Complex factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarDiv( factor );
      }
    }
  }


  // ************************
  //   ScalarMult (complex)
  // ************************
  void SBM_Vector::ScalarMult( const Complex factor ) {


    for ( UInt i = 0; i < size_; i++ ) {
      if ( subVec_[i] != NULL ) {
        subVec_[i]->ScalarMult( factor );
      }
    }
  }

}
