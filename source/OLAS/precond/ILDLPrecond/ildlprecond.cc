#include <vector>

#include "MatVec/opdefs.hh"
#include "MatVec/StdMatrix.hh"
#include "MatVec/SCRS_Matrix.hh"
#include "DataInOut/ProgramOptions.hh"
#include "ildlprecond.hh"

// By default, we assume that, if this class is to be debugged, then also the
// classes which implement the actual incomplete factorisation variants are to
// be debugged.
#ifdef DEBUG_ILDLPRECOND
#define DEBUG_ILDLKFACTORISER
#define DEBUG_ILDLTPFACTORISER
#define DEBUG_ILDLCNFACTORISER
#endif

// include source code of the factorisers
#include "ildl0factoriser.hh"
#include "ildlkfactoriser.hh"
#include "ildltpfactoriser.hh"
#include "ildlcnfactoriser.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  template<typename T>
  ILDLPrecond<T>::ILDLPrecond( const StdMatrix &stdMat, PtrParamNode precondNode,
                               PtrParamNode olasInfo ) {

    // The ILDL preconditioner is currently not working
    WARN("The ILDL-Preconditioner is not yet ported from 1 to 0 based"
        "numbering, i.e. it will NOT work correctly!" );
        
    // Set pointers to communication objects
    this->xml_ = precondNode;
    this->infoNode_ = olasInfo->Get("ildl", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);

    // No factorisation was computed yet
    amFactorised_ = false;

    // No problem size known yet
    this->sysMatDim_ = 0;

    // Obtain and check variant information
    std::string subTypeStr = "NOPRECOND";
    
    myVariant_ = BasePrecond::precondType.Parse(subTypeStr);
    if ( myVariant_ != BasePrecond::ILDL0 &&
         myVariant_ != BasePrecond::ILDLTP &&
         myVariant_ != BasePrecond::ILDLK &&
         myVariant_ != BasePrecond::ILDLCN ) {
      
      std::string tmp;
      tmp = BasePrecond::precondType.ToString(myVariant_);
      
      EXCEPTION( "ILDLPrecond: The constructor detected an error in the "
               << "'ILDLPRECOND_subType' parameter! The value '"
               << tmp << "' is not a valid variant." );
    }

    this->xml_->GetValue("logging", logging_, ParamNode::INSERT ) ;
    
    // Generate factorisation object
    GenerateFactoriser();

  }


  // **************
  //   Destructor
  // **************
  template<typename T>
  ILDLPrecond<T>::~ILDLPrecond() {


    delete factoriser_;
  }


  // *********
  //   Setup
  // *********
  template<typename T>
  void ILDLPrecond<T>::Setup( BaseMatrix &mat ) {


    // Test the storage layout
    StdMatrix & stdMat = dynamic_cast<StdMatrix&>(mat);
    BaseMatrix::StorageType sType = stdMat.GetStorageType();
    if ( sType != BaseMatrix::SPARSE_SYM ) {
      EXCEPTION( "ILDLPrecond::Setup: The ILDLPrecond requires the system "
               << "matrix to be an SCRS_Matrix, i.e. sparseSym. The system "
               << "matrix you supplied is a matrix in '"
               << MatrixStorageTypeEnum.ToString( sType )
               << "' format." );
    }

    SCRS_Matrix<T>& scrsMat = dynamic_cast<SCRS_Matrix<T>&>(stdMat);


    // Get new problem size and perform consistency check
    if ( amFactorised_ == false ) {
      this->sysMatDim_ = scrsMat.GetNumCols();
    }
    else {
      // TODO: THIS CHECK DOES NOT MAKE SENSE IN MY OPINION SINCE
      //       'newMatrixPattern' is set to false in olasparams.cc
      //       and gets never changed elsewhere. A more intelligent
      //       test would ask the matrix if its pattern did change.

      bool newMatrixPattern = false;
      
      if ( newMatrixPattern == false &&
          this->sysMatDim_ != scrsMat.GetNumCols() ) {
        EXCEPTION( "ILDLPrecond::Setup: newMatrixPattern = false, but "
            << "matrix dimension changed from " << this->sysMatDim_ << " to "
            << scrsMat.GetNumCols() );
      }
      else {
        this->sysMatDim_ = scrsMat.GetNumCols();
        amFactorised_ = false;
      }
    }

    // Let factoriser compute the factorisation
    factoriser_->Factorise( scrsMat, dataD_, rptrU_, cidxU_, dataU_, true );
    amFactorised_ = true;

    // Convert D to D^{-1}
    for ( UInt i = 1; i <= this->sysMatDim_; i++ ) {
      dataD_[i] = OpType<T>::invert( dataD_[i] );
    }

    // If the user desires it, we will export the matrix factor
    // to a file in Matrix Market format
    std::string saveFacFile = "ildl_fac.out";
    std::string subTypeStr;
    
    subTypeStr = BasePrecond::precondType.ToString(myVariant_);

    if(this->xml_->Has("saveFacFile")) {
      bool pattern = false;
      this->xml_->GetValue("savePatternOnly", pattern, ParamNode::INSERT);
      this->xml_->GetValue("saveFacFile", saveFacFile, ParamNode::INSERT);
      ExportFactorisation( saveFacFile.c_str(), pattern );
    }
  }


  // *********
  //   Apply
  // *********
  template<typename T>
  void ILDLPrecond<T>::Apply( const BaseMatrix &stdMat, const BaseVector &r,
                              BaseVector &z ) const {



    // Test that a factorisation is available, if not issue an error
    if ( amFactorised_ == false ) {
      EXCEPTION( "ILDLPrecond::Apply: No factorisation available. "
               << "Call Setup() first!" );
    }

    // Solve the problem
    const Vector<T>& myR = dynamic_cast<const Vector<T>&>(r);
    Vector<T>& myZ = dynamic_cast<Vector<T>&>(z);

    this->SolveLDLSystem( &(cidxU_[0]), &(rptrU_[0]), &(dataU_[0]),
        &(dataD_[0]), myZ, myR, this->sysMatDim_ );
  }


  // ***********************
  //   ExportFactorisation
  // ***********************
  template <typename T>
  void ILDLPrecond<T>::ExportFactorisation( const char *fname,
                                            bool patternOnly ) {


    UInt i, j;
    T aux;

    // ====================
    //   Open output file
    // ====================
    FILE *fp = fopen( fname, "w" );
    if ( fp == NULL ) {
      EXCEPTION( "ILDLPrecond::ExportFactorisation: Cannot open file "
               << fname << " for writing!" );
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
      if ( EntryType<T>::M_EntryType == BaseMatrix::DOUBLE ) {
        myFormat = "real";
      }
      else if ( EntryType<T>::M_EntryType == BaseMatrix::COMPLEX ) {
        myFormat = "complex";
      }
      else {
        EXCEPTION( "ILDLPrecond::ExportFactorisation: Cannot identify "
            << "template parameter" );
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
        aux = OpType<T>::invert( dataD_[i] );
        OpType<T>::ExportEntry( aux, 0, 0, fp );
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
          OpType<T>::ExportEntry( dataU_[j], 0, 0, fp );
        }
        fprintf( fp, "\n" );
      }
    }

    // =====================
    //   Close output file
    // =====================
    if ( fclose( fp ) == EOF ) {
      WARN("ILDLPrecond::ExportFactorisation: Could not close file "
           << fname << " after writing!");
    }
  }


  // **********************
  //   GenerateFactoriser
  // **********************
  template<typename T>
  void ILDLPrecond<T>::GenerateFactoriser() {


    switch ( myVariant_ ) {

      // ILDLK factoriser
    case BasePrecond::ILDL0:
      factoriser_ = new ILDL0Factoriser<T>( this->xml_, this->infoNode_ );
      break;

      // ILDLK factoriser
    case BasePrecond::ILDLK:
      factoriser_ = new ILDLKFactoriser<T>( this->xml_, this->infoNode_ );
      break;

      // ILDLTP factoriser
    case BasePrecond::ILDLTP:
      factoriser_ = new ILDLTPFactoriser<T>( this->xml_, this->infoNode_ );
      break;

      // ILDLK factoriser
    case BasePrecond::ILDLCN:
      factoriser_ = new ILDLCNFactoriser<T>( this->xml_, this->infoNode_ );
      break;

      // Wrong preconditioner type
    default:
      std::string tmp;
      tmp = BasePrecond::precondType.ToString(myVariant_);      
      
      EXCEPTION( "GenerateFactoriserObject: preconditioner type = '"
          << tmp << "' is no valid variant of an "
          << "ILDL preconditioner (internal error: the constructor "
          << "should already have detected this!)" );
    }
  }

  // Explicit template instantiation
  template class ILDLPrecond<Double>;
  template class ILDLPrecond<Complex>;
  
}
