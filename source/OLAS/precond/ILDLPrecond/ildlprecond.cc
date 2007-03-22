// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <vector>

#include "precond/ILDLPrecond/ildlprecond.hh"

// By default, we assume that, if this class is to be debugged, then also the
// classes which implement the actual incomplete factorisation variants are to
// be debugged.
#ifdef DEBUG_ILDLPRECOND
#define DEBUG_ILDLKFACTORISER
#define DEBUG_ILDLTPFACTORISER
#define DEBUG_ILDLCNFACTORISER
#endif

// include source code of the factorisers
#include "precond/ILDLPrecond/ildl0factoriser.cc"
#include "precond/ILDLPrecond/ildlkfactoriser.cc"
#include "precond/ILDLPrecond/ildltpfactoriser.cc"
#include "precond/ILDLPrecond/ildlcnfactoriser.cc"

namespace OLAS {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  ILDLPrecond<T>::ILDLPrecond( const StdMatrix &stdMat, OLAS_Params *myParams,
                               OLAS_Report *myReport ) {

    ENTER_FCN( "ILDLPrecond::ILDLPrecond" );

    // Set pointers to communication objects
    this->myParams_ = myParams;
    this->myReport_ = myReport;

    // No factorisation was computed yet
    amFactorised_ = false;

    // No problem size known yet
    this->sysMatDim_ = 0;

    // Obtain and check variant information
    this->myParams_->GetEnumValue( "ILDLPRECOND_subType", myVariant_ );
    if ( myVariant_ != ILDL0 && myVariant_ != ILDLTP 
	 && myVariant_ != ILDLK && myVariant_ != ILDLCN ) {
      (*error) << "ILDLPrecond: The constructor detected an error in the "
               << "'ILDLPRECOND_subType' parameter! The value '"
               << Enum2String( myVariant_ ) << "' is not a valid variant.";
      Error( __FILE__, __LINE__ );
    }

    // Generate factorisation object
    GenerateFactoriser();

  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  ILDLPrecond<T>::~ILDLPrecond() {

    ENTER_FCN( "ILDLPrecond::~ILDLPrecond" );

    delete factoriser_;
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void ILDLPrecond<T>::Setup( StdMatrix &stdMat ) {

    ENTER_FCN( "ILDLPrecond::Setup" );

    // Test the storage layout
    MatrixStorageType sType = stdMat.GetStorageType();
    if ( sType != SPARSE_SYM ) {
      (*error) << "ILDLPrecond::Setup: The ILDLPrecond requires the system "
               << "matrix to be an SCRS_Matrix, i.e. sparseSym. The system "
               << "matrix you supplied is a matrix in '"
               << Enum2String( sType )
               << "' format.";
      Error( __FILE__, __LINE__ );
    }

    TRY_CAST {

      // Down-cast to SCRS_Matrix
      RefCast( stdMat, SCRS_Matrix<T>, scrsMat );

      // Get new problem size and perform consistency check
      if ( amFactorised_ == false ) {
        this->sysMatDim_ = scrsMat.GetNcols();
      }
      else {
        if ( this->myParams_->GetBoolValue( "newMatrixPattern" ) == false &&
             this->sysMatDim_ != scrsMat.GetNcols() ) {
          (*error) << "ILDLPrecond::Setup: newMatrixPattern = false, but "
                   << "matrix dimension changed from " << this->sysMatDim_ << " to "
                   << scrsMat.GetNcols();
          Error( __FILE__, __LINE__ );
        }
        else {
          this->sysMatDim_ = scrsMat.GetNcols();
          amFactorised_ = false;
        }
      }

      // Logging
      bool logging = this->myParams_->GetIntValue( "ILDLPRECOND_logging" ) > 0;
      if ( logging ) {
        (*cla) << " -------------------------------------------------------"
               << "-----------------------\n"
               << " ILDLPRECOND::SETUP\n + Factorisation of a "
               << scrsMat.GetNrows() << " x " << scrsMat.GetNcols()
               << " matrix (nnz = " << scrsMat.GetNnz() << ")"
               << std::endl;
      }

      // Let factoriser compute the factorisation
      factoriser_->Factorise( scrsMat, dataD_, rptrU_, cidxU_, dataU_, true );
      amFactorised_ = true;

      // Convert D to D^{-1}
      for ( UInt i = 1; i <= this->sysMatDim_; i++ ) {
        dataD_[i] = opType<T>::invert( dataD_[i] );
      }

      // If the user desires it, we will export the matrix factor
      // to a file in Matrix Market format
      if( this->myParams_->GetBoolValue( "ILDLPRECOND_saveFacToFile" ) == true ) {
        bool pattern = this->myParams_->GetBoolValue( "ILDLPRECOND_savePatternOnly");
        std::string filename;
        filename = this->myParams_->GetStringValue( "ILDLPRECOND_facFileName" );
        ExportFactorisation( filename.c_str(), pattern );
      }

      // finish log report
      if ( logging ) {

        (*cla) << "\n Change of fill-pattern:"
               << "\n + nnz in upper triangle of A: "
               << (scrsMat.GetNnz()+this->sysMatDim_) / 2
               << "\n + nnz in F = D + L^T:         "
               << rptrU_[this->sysMatDim_-1] + this->sysMatDim_ << " = "
               << this->sysMatDim_ << " + " << rptrU_[this->sysMatDim_-1]
               << "\n Memory requirements:"
               << "\n + array size (dataU):         " << cidxU_.size()
               << "\n + array capacity (dataU):     " << cidxU_.capacity()
               << '\n';
        (*cla) << " -------------------------------------------------------"
               << "-----------------------\n" << std::endl;
      }

    } CATCH_CAST;

  }


  // *********
  //   Apply
  // *********
  template<typename T>
  void ILDLPrecond<T>::Apply( const StdMatrix &stdMat, const StdVector &r,
                              StdVector &z ) const {

    ENTER_FCN( "ILDLPrecond::Solve" );

    bool logging = this->myParams_->GetIntValue( "ILDLPRECOND_logging" ) > 1;

    // Test that a factorisation is available, if not issue an error
    if ( amFactorised_ == false ) {
      (*error) << "ILDLPrecond::Apply: No factorisation available. "
               << "Call Setup() first!";
      Error( __FILE__, __LINE__ );
    }

    // Solve the problem
    TRY_CAST {
      ConstRefCast( r, Vector<T>, myR );
      RefCast( z, Vector<T>, myZ );

      // Logging
      if ( logging ) {
        (*cla) << " -------------------------------------------------------"
               << "-----------------------\n"
               << " ILDLPRECOND::APPLY: Solving linear system with "
               << stdMat.GetNcols() << " unknowns" << std::endl;
      }

      SolveLDLSystem( &(cidxU_[0]), &(rptrU_[0]), &(dataU_[0]),
                      &(dataD_[0]), myZ, myR, this->sysMatDim_ );

      // Logging
      if ( logging ) {
        (*cla) << " -------------------------------------------------------"
               << "-----------------------\n"
               << std::endl;
      }

    } CATCH_CAST;
  }


  // ***********************
  //   ExportFactorisation
  // ***********************
  template <typename T>
  void ILDLPrecond<T>::ExportFactorisation( const char *fname,
                                            bool patternOnly ) {

    ENTER_FCN( "ILDLPrecond::ExportFactorisation" );

    UInt i, j;
    T aux;

    // ====================
    //   Open output file
    // ====================
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      (*error) << "ILDLPrecond::ExportFactorisation: Cannot open file "
               << fname << " for writing!";
      Error( __FILE__, __LINE__ );
    }

    // =====================
    //   Write file header
    // =====================

    // Matrix Market Format Specification
    std::string myFormat;
    if ( patternOnly == true ) {
      myFormat = "pattern";
    }
    else {
      if ( EntryType<T>::M_EntryType == DOUBLE ) {
        myFormat = "real";
      }
      else if ( EntryType<T>::M_EntryType == COMPLEX ) {
        myFormat = "complex";
      }
      else {
        (*error) << "ILDLPrecond::ExportFactorisation: Cannot identify "
                 << "template parameter";
        Error( __FILE__, __LINE__ );
      }
    }
    fprintf( fp, "%%%%MatrixMarket matrix coordinate %s general\n",
             myFormat.c_str() );

    // Comment
    if ( patternOnly == true ) {
      fprintf( fp, "%%\n%% Sparsity pattern of LDL^T factorisation matrix " );
      fprintf( fp, "F = D + L^T\n" );
      fprintf( fp, "%% computed by ILDLPrecond\n%%\n" );
    }
    else {
      fprintf( fp, "%%\n%% LDL^T factorisation matrix F = D + L^T " );
      fprintf( fp, " computed by ILDLPrecond\n%%\n" );
    }

    // Information on number of rows, columns and entries
    fprintf( fp, "%d\t%d\t%d\n", this->sysMatDim_, this->sysMatDim_,
             rptrU_[this->sysMatDim_+1] - 1 + this->sysMatDim_ );

    // ======================
    //   Write entries of D
    // ======================
    for ( i = 1; i <= this->sysMatDim_; i++ ) {
      fprintf( fp, "%6d\t%6d\t", i, i );
      if ( patternOnly == false ) {
        aux = opType<T>::invert( dataD_[i] );
        opType<T>::ExportEntry( aux, 0, 0, fp );
      }
      fprintf( fp, "\n" );
    }

    // ======================
    //   Write entries of U
    // ======================
    for ( i = 1; i <= this->sysMatDim_; i++ ) {
      for ( j = rptrU_[i]; j < rptrU_[i+1]; j++ ) {
        fprintf( fp, "%6d\t%6d\t", i, cidxU_[j] );
        if ( patternOnly == false ) {
          opType<T>::ExportEntry( dataU_[j], 0, 0, fp );
        }
        fprintf( fp, "\n" );
      }
    }

    // =====================
    //   Close output file
    // =====================
    if ( fclose( fp ) == EOF ) {
      (*warning) << "ILDLPrecond::ExportFactorisation: Could not close file "
                 << fname << " after writing!";
      Warning( __FILE__, __LINE__ );
    }
  }


  // **********************
  //   GenerateFactoriser
  // **********************
  template<typename T>
  void ILDLPrecond<T>::GenerateFactoriser() {

    ENTER_FCN( "ILDLPrecond::GenerateFactoriser" );

    switch ( myVariant_ ) {

      // ILDLK factoriser
    case ILDL0:
      factoriser_ = new ILDL0Factoriser<T>( this->myParams_, this->myReport_ );
      break;

      // ILDLK factoriser
    case ILDLK:
      factoriser_ = new ILDLKFactoriser<T>( this->myParams_, this->myReport_ );
      break;

      // ILDLTP factoriser
    case ILDLTP:
      factoriser_ = new ILDLTPFactoriser<T>( this->myParams_, this->myReport_ );
      break;

      // ILDLK factoriser
    case ILDLCN:
      factoriser_ = new ILDLCNFactoriser<T>( this->myParams_, this->myReport_ );
      break;

      // Wrong preconditioner type
    default:
      (*error) << "GenerateFactoriserObject: preconditioner type = '"
               << Enum2String( myVariant_ ) << "' is no valid variant of an "
               << "ILDL preconditioner (internal error: the constructor "
               << "should already have detected this!)";
      Error( __FILE__, __LINE__ );
    }
  }

}
