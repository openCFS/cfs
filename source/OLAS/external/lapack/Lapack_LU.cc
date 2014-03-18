#include <cstring>

#include "MatVec/BLASLAPACKInterface.hh"

#include "OLAS/external/lapack/Lapack_LU.hh"

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  Lapack_LU::Lapack_LU() {
    EXCEPTION( "Default constructor of LAPACK_LU is forbidden!" );
  }


  // *****************************
  //   Constructor (for factory)
  // *****************************
  Lapack_LU::Lapack_LU( PtrParamNode solverNode, PtrParamNode olasInfo )
    : pivots_(NULL), facmat_(NULL), amFactorised_(false) {

    xml_ = solverNode;
    infoNode_ = olasInfo->Get("lapackLU");
    
    // Initialise pointers for LAPACK workspaces
    workspaceF77REAL8_     = NULL;
    workspaceF77COMPLEX16_ = NULL;
    workspaceInt_          = NULL;

    // Initialise attributes related to scaling
    row_scalings_ = NULL;
    col_scalings_ = NULL;
    scalingType_ = 'N';

  }


  // *************************
  //   Alternate Constructor
  // *************************
  Lapack_LU::Lapack_LU( BaseMatrix &mat,
                        PtrParamNode solverNode, 
			                  PtrParamNode olasInfo )
    : pivots_(NULL), facmat_(NULL), amFactorised_(false) {

    xml_ = solverNode;
    infoNode_ = olasInfo->Get("lapackLU");

    // Initialise pointers for LAPACK workspaces
    workspaceF77REAL8_     = NULL;
    workspaceF77COMPLEX16_ = NULL;
    workspaceInt_          = NULL;

    // Perform factorisation
    Setup( mat );
  }


  // ***************************
  //   Deep Default Destructor
  // ***************************
  Lapack_LU::~Lapack_LU() {
    delete facmat_;
    delete[] pivots_;
    delete[] workspaceF77REAL8_;
    delete[] workspaceF77COMPLEX16_;
    delete[] workspaceInt_;
    delete[] row_scalings_;
    delete[] col_scalings_;
  }


  // ********************************
  //   Setup: Perform Factorisation
  // ********************************
  void Lapack_LU::Setup( BaseMatrix &sysmat, PtrParamNode analysis_step ) {
    PrivateSetup( sysmat );
  }


  // ********************************
  //   Setup: Perform Factorisation
  // ********************************
  void Lapack_LU::PrivateSetup( const BaseMatrix &sysmat ) {


    // Report to logfile
    (*cla) << " --------------------------------------\n"
	   << " LAPACK_LU: Starting factorisation" << std::endl;

    const StdMatrix& stdmat = dynamic_cast<const StdMatrix&>(sysmat);
      // Check that we have the correct matrix type
      BaseMatrix::StorageType mtype = stdmat.GetStorageType();
      if ( mtype != BaseMatrix::LAPACK_GBMATRIX ) {
        EXCEPTION( "Expected a LAPACK_GBMATRIX" );
      }

      // Allocate memory for pivot indices
      if ( amFactorised_ == true ) {
	delete[] pivots_;
      }
      pivots_ = new int[stdmat.GetNumCols()];
      ASSERTMEM(pivots_,stdmat.GetNumCols());

      // Get the entry type to figure out which Factorisation method to call
      BaseMatrix::EntryType etype = stdmat.GetEntryType();

      // Call appropriate factorisation routine
      switch( etype ) {

      case BaseMatrix::F77REAL8:
	FactoriseF77REAL8( stdmat );
	break;

      case BaseMatrix::F77COMPLEX16:
	FactoriseF77COMPLEX16( stdmat );
	break;

      default:
        EXCEPTION( "Matrix entry type not valid for a LAPACK matrix" );
      }

    // now we have a (new) factorisation
    amFactorised_ = true;

    // Report to logfile
    (*cla) << " LAPACK_LU: Finished factorisation\n"
	   << " --------------------------------------" << std::endl;

  }


  // *******************************************
  //   Factorisation (real & double precision)
  // *******************************************
  void Lapack_LU::FactoriseF77REAL8( const StdMatrix &stdmat ) {


    // Down-cast matrix
    const LapackGBMatrix<double,double> &mat =
      dynamic_cast<const LapackGBMatrix<double,double>&>(stdmat);

    // Obtain some matrix information
    int lp_ncols = (int)mat.GetNumCols();
    int lp_nrows = (int)mat.GetNumRows();
    int lp_wlower = (int)mat.GetLowerBandwidth();
    int lp_wupper = (int)mat.GetUpperBandwidth();
    int lp_ldab = lp_wlower + lp_wupper + 1;
    double *matdata = mat.GetDataPointer0();

    // status flag
    int lp_info = 0;

    // =======================================================
    //   Scale matrix to improve condition number (optional)
    // =======================================================
    bool tryScaling = true;
      
    if(xml_) {
      xml_->GetValue("tryScaling", tryScaling, ParamNode::INSERT );
    }
    
    if ( tryScaling ) {

      // Allocate memory for scaling factors
      row_scalings_ = new double[lp_nrows];
      ASSERTMEM( row_scalings_, lp_nrows );
      col_scalings_ = new double[lp_ncols];
      ASSERTMEM( col_scalings_, lp_ncols );

      // output parameters (condition numbers)
      double lp_rowcond;
      double lp_colcond;
      double lp_amax;

      // First compute a possible scaling
      dgbequ( &lp_nrows, &lp_ncols, &lp_wlower, &lp_wupper,
              matdata, &lp_ldab, row_scalings_, col_scalings_,
              &lp_rowcond, &lp_colcond, &lp_amax, &lp_info );

      // Process status flag
      if ( lp_info < 0 ) {
        EXCEPTION( "DGBEQU reports invalid input parameter" );
      }
      else if ( lp_info > lp_nrows ) {
        EXCEPTION( "DBGEQU reports that column "
            << lp_info - lp_nrows << " is exactly zero" );
      }
      else if ( lp_info > 0 ) {
        EXCEPTION( "DBGEQU reports that row " << lp_info
            << " is exactly zero" );
      }
      else {
        (*cla) << " DGBEQU reports: rowcnd = "   << lp_rowcond
	         << "\n                 colcnd = " << lp_colcond
	         << "\n                 amax   = " << lp_amax
	         << std::endl;
      }

      // Scale matrix
      char lp_equed;
      dlaqgb( &lp_nrows, &lp_ncols, &lp_wlower, &lp_wupper,
              matdata, &lp_ldab, row_scalings_, col_scalings_,
              &lp_rowcond, &lp_colcond, &lp_amax, &lp_equed );

      // Process result
      if ( lp_equed == 'N' ) {
	(*cla) << " DLAQGB performed no scaling" << std::endl;
      }
      else if ( lp_equed == 'R' ) {
	(*cla) << " DLAQGB performed row scaling" << std::endl;
      }
      else if ( lp_equed == 'C' ) {
	(*cla) << " DLAQGB performed column scaling" << std::endl;
      }
      else if ( lp_equed == 'B' ) {
	(*cla) << " DLAQGB performed row & column scaling" << std::endl;
      }
      else {
        EXCEPTION( "DLAQGB returned with unexpected parameter value" );
      }
      (*cla) << " DLAQGB reports: rowcnd = "   << lp_rowcond
	     << "\n                 colcnd = " << lp_colcond
	     << "\n                 amax   = " << lp_amax
	     << std::endl;

      // Save type of scaling performed for re-scaling of solution
      scalingType_ = lp_equed;

    }
    
    // ========================
    //   Generate matrix copy
    // ========================

    // We have to generate a copy of the matrix containing additional
    // wlower_ super-diagonals in the beginning as work-space and for
    // the result of the factorisation
    LapackGBMatrix<double,double> *facmat = new
    LapackGBMatrix<double,double>( lp_nrows, lp_ncols, lp_wlower,
                                     lp_wupper, BaseMatrix::F77REAL8 );
    if ( facmat == NULL ) {
      EXCEPTION( "Memory allocation for new LapackGBMatrix failed" );
    }

    // Upcast matrix to store as class attribute
    facmat_ = dynamic_cast<LapackBaseMatrix*>(facmat);

    // Copy original matrix to rows wlower+1 - 2*wlower+wupper+1
    // Since the data_ array is in FORTRAN storage layout we do this by
    // copying the individual columns of the original matrix into the lower
    // parts of the corresponding columns in the new matrix
    double *facmatdata = facmat->GetDataPointer0();
    Integer off1 = facmat->nrowsact_;
    Integer off2 = mat.nrowsact_;

    for ( Integer j = 0; j < lp_ncols; j++ ) {

      // copy j-th column
      std::memcpy( facmatdata + j * off1 + lp_wlower, matdata + j * off2,
		   off2 * sizeof(double) );
    }

    // ============================
    //   Compute LU factorisation
    // ============================

    // Call factorisation routine
    lp_ldab = 2 * lp_wlower + lp_wupper + 1;
    dgbtrf( &lp_nrows, &lp_ncols, &lp_wlower, &lp_wupper, facmatdata,
            &lp_ldab, pivots_, &lp_info );

    // Process status flag
    if ( lp_info < 0 ) {
      EXCEPTION( "DGBTRF reports invalid input parameter" );
    }
    else if ( lp_info > 0 ) {

#ifdef DEBUG_LAPACK_LU
      (*debug) << " Problems: DGBTRF returned with info = " << lp_info
	       << std::endl;
      (*debug) << " Matrix entry is: "
	       << facmatdata[facmat->Index(lp_info,lp_info)-1]
	       << std::endl;
#endif

      EXCEPTION( "DBGTRF reports zero pivot element" );
    }

  }


  // **********************************************
  //   Factorisation (complex & double precision)
  // **********************************************
  void Lapack_LU::FactoriseF77COMPLEX16( const StdMatrix &stdmat ) {


    // Down-cast matrix
    const LapackGBMatrix<std::complex<double>,Complex> &mat =
      dynamic_cast<const LapackGBMatrix<std::complex<double>,Complex>&>(stdmat);

    // Obtain some matrix information
    int lp_ncols = (int)mat.GetNumCols();
    int lp_nrows = (int)mat.GetNumRows();
    int lp_wlower = (int)mat.GetLowerBandwidth();
    int lp_wupper = (int)mat.GetUpperBandwidth();
    int lp_ldab = lp_wlower + lp_wupper + 1;
    std::complex<double> *matdata = mat.GetDataPointer0();

    // status flag
    int lp_info = 0;

    // =======================================================
    //   Scale matrix to improve condition number (optional)
    // =======================================================

    bool tryScaling = true;
      
    if(xml_) {
      xml_->GetValue("tryScaling", tryScaling, ParamNode::INSERT );
    }
    
    if ( tryScaling ) {

      // Allocate memory for scaling factors
      row_scalings_ = new double[lp_nrows];
      ASSERTMEM( row_scalings_, lp_nrows );
      col_scalings_ = new double[lp_ncols];
      ASSERTMEM( col_scalings_, lp_ncols );

      // output parameters (condition numbers)
      double lp_rowcond;
      double lp_colcond;
      double lp_amax;

      // First compute a possible scaling
      zgbequ( &lp_nrows, &lp_ncols, &lp_wlower, &lp_wupper,
              matdata, &lp_ldab, row_scalings_, col_scalings_,
              &lp_rowcond, &lp_colcond, &lp_amax, &lp_info );

      // Process status flag
      if ( lp_info < 0 ) {
        EXCEPTION( "ZGBEQU reports invalid input parameter" );
      }
      else if ( lp_info > lp_nrows ) {
        EXCEPTION( "ZBGEQU reports zero column" );
      }
      else if ( lp_info > 0 ) {
        EXCEPTION( "ZBGEQU reports zero row" );
      }
      else {
	(*cla) << " ZGBEQU reports: rowcnd = "   << lp_rowcond
	       << "\n                 colcnd = " << lp_colcond
	       << "\n                 amax   = " << lp_amax
	       << std::endl;
      }

      // Scale matrix
      char lp_equed;
      zlaqgb( &lp_nrows, &lp_ncols, &lp_wlower, &lp_wupper,
              matdata, &lp_ldab, row_scalings_, col_scalings_,
              &lp_rowcond, &lp_colcond, &lp_amax, &lp_equed );

      // Process result
      if ( lp_equed == 'N' ) {
	(*cla) << " ZLAQGB performed no scaling" << std::endl;
      }
      else if ( lp_equed == 'R' ) {
	(*cla) << " ZLAQGB performed row scaling" << std::endl;
      }
      else if ( lp_equed == 'C' ) {
	(*cla) << " ZLAQGB performed column scaling" << std::endl;
      }
      else if ( lp_equed == 'B' ) {
	(*cla) << " ZLAQGB performed row & column scaling" << std::endl;
      }
      else {
        EXCEPTION( "ZLAQGB returned with unexpected parameter value" );
      }
      (*cla) << " ZLAQGB reports: rowcnd = "   << lp_rowcond
             << "\n                 colcnd = " << lp_colcond
             << "\n                 amax   = " << lp_amax
	     << std::endl;

      // Save type of scaling performed for re-scaling of solution
      scalingType_ = lp_equed;

    }

    // ========================
    //   Generate matrix copy
    // ========================

    // We have to generate a copy of the matrix containing additional
    // wlower_ super-diagonals in the beginning as work-space and for
    // the result of the factorisation
    LapackGBMatrix<std::complex<double>,Complex> *facmat = new
      LapackGBMatrix<std::complex<double>,Complex>( lp_nrows, lp_ncols, lp_wlower,
					    lp_wupper, BaseMatrix::F77COMPLEX16 );
    if ( facmat == NULL ) {
      EXCEPTION( "Memory allocation for new LapackGBMatrix failed" );
    }

    // Upcast matrix to store as class attribute
    facmat_ = dynamic_cast<LapackBaseMatrix*>(facmat);

    // Copy original matrix to rows wlower+1 - 2*wlower+wupper+1
    // Since the data_ array is in FORTRAN storage layout we do this by
    // copying the individual columns of the original matrix into the lower
    // parts of the corresponding columns in the new matrix
    std::complex<double> *facmatdata = facmat->GetDataPointer0();
    Integer off1 = facmat->nrowsact_;
    Integer off2 = mat.nrowsact_;

    for ( Integer j = 0; j < lp_ncols; j++ ) {

      // copy j-th column
      std::memcpy( facmatdata + j * off1 + lp_wlower, matdata + j * off2,
		   off2 * sizeof(std::complex<double>) );
    }

    // ============================
    //   Compute LU factorisation
    // ============================

    // Call factorisation routine
    lp_ldab = 2 * lp_wlower + lp_wupper + 1;
    zgbtrf( &lp_nrows, &lp_ncols, &lp_wlower, &lp_wupper, facmatdata,
            &lp_ldab, pivots_, &lp_info );

    // Process status flag
    if ( lp_info < 0 ) {
      EXCEPTION( "ZGBTRF reports invalid input parameter" );
    }
    else if ( lp_info > 0 ) {

#ifdef DEBUG_LAPACK_LU
      (*debug) << " Problems: ZGBTRF returned with info = " << lp_info
	       << std::endl;
      (*debug) << " Matrix entry is: "
	       << facmatdata[facmat->Index(lp_info,lp_info)-1]
	       << std::endl;
#endif

      EXCEPTION( "ZBGTRF reports zero pivot element" );
    }

  }


  // ***********************
  //   Solve linear system
  // ***********************
  void Lapack_LU::Solve( const BaseMatrix &sysmat, 
      const BaseVector &rhs, BaseVector &sol, PtrParamNode analysis_step ) {


    // Are we expected to be verbose?
    bool logging = false;

    // Report to logfile
    if ( logging == true ) {
      (*cla) << " --------------------------------------\n"
      << " LAPACK_LU: Computing solution" << std::endl;
    }

    if ( facmat_ == NULL || amFactorised_ == false ) {

      // If the two indicators are consistent call Setup()
      if ( facmat_ == NULL && amFactorised_ == false ) {
        PrivateSetup( sysmat );
      }

      // The two indicators disagree, so complain
      else if ( amFactorised_ == false ) {
        EXCEPTION( "LAPACKLU: Internal error. facmat_ <> NULL but amFactorised_ = "
            "false!" );
      }
      else {
        EXCEPTION( "LAPACKLU: Internal error. facmat_ = NULL but amFactorised_ = "
            "true!" );
      }
    }

    const StdMatrix& stdmat = dynamic_cast<const StdMatrix&>(sysmat);
    // Check that we have the correct matrix type
    BaseMatrix::StorageType mtype = stdmat.GetStorageType();
    if ( mtype != BaseMatrix::LAPACK_GBMATRIX ) {
      EXCEPTION( "Expected a LAPACK_GBMATRIX" );
    }

    // Get the entry type to figure out which Factorization method to call
    BaseMatrix::EntryType etype = stdmat.GetEntryType();

    // Call appropriate solution routine
    switch( etype ) {

    case BaseMatrix::F77REAL8:
      SolveF77REAL8( rhs, sol, &sysmat );
      break;

    case BaseMatrix::F77COMPLEX16:
      SolveF77COMPLEX16( rhs, sol, &sysmat );
      break;

    default:
      EXCEPTION( "Matrix entry type not valid for a LAPACK matrix" );
    }

    // Report to logfile
    if ( logging == true ) {
      (*cla) << " --------------------------------------" << std::endl;
    }

  }


  // *********************************************
  //   Solution method (real & double precision)
  // *********************************************
  void Lapack_LU::SolveF77REAL8( const BaseVector &rhs, BaseVector &sol,
				 const BaseMatrix *mat ) {


    // Are we expected to be verbose?
    bool logging = false;

    // Some variables for LAPACK
    char lp_trans = 'N';
    int lp_info = 0;
    int lp_one = 1;

    // Downcast matrices
    LapackGBMatrix<double,double> *facmat =
      dynamic_cast<LapackGBMatrix<double,double> *>(facmat_);
    const LapackGBMatrix<double,double> *mymat =
      dynamic_cast<const LapackGBMatrix<double,double> *>(mat);

    // This should not have failed, but better test
    if ( facmat == NULL ) {
      EXCEPTION( "Downcast of facmat_ failed!" );
    }
    if ( facmat == NULL ) {
      EXCEPTION( "Downcast of system matrix failed!" );
    }

    // Obtain data pointers
    double *facmatdata = facmat->GetDataPointer0();
    const double *matdata = mymat->GetDataPointer0();

    // Obtain some matrix information
    int lp_ncols  = (int)facmat->GetNumCols();
    int lp_wlower = (int)facmat->GetLowerBandwidth();
    int lp_wupper = (int)facmat->GetUpperBandwidth();
    int lp_ldab   = lp_wlower + lp_wupper + 1;
    int lp_ldabf  = 2 * lp_wlower + lp_wupper + 1;

    // Downcast vectors and get data pointers
    const double *lp_rhs;
    double *lp_sol;
    const Vector<Double>& myrhs = dynamic_cast<const Vector<Double>&>(rhs);
    Vector<Double>& mysol = dynamic_cast<Vector<Double>&>(sol);

    lp_rhs = myrhs.GetPointer() + 1;
    lp_sol = mysol.GetPointer() + 1;

    // ====================
    //   Compute solution
    // ====================
   
    // Copy rhs entries into solution vector

    // If row scaling occured, the system matrix was replaced by 
    // diag(R) * A, where R is the vector of row scaling factors.
    // We will get the correct solution, if we transform also the
    // right hand side and solve the system diag(R) * A * x = diag(R) * y.
    if ( scalingType_ == 'R' || scalingType_ == 'B' ) {
      for ( Integer i = 0; i < lp_ncols; i++ ) {
//	std::cerr << "R(" << i << ") = " << row_scalings_[i] << std::endl;
	lp_sol[i] = row_scalings_[i] * lp_rhs[i];
      }
    }

    // If now scaling was performed we simply copy the vector
    else {
      for ( Integer i = 0; i < lp_ncols; i++ ) {
	lp_sol[i] = lp_rhs[i];
      }
    }

    // Perform backward/forward substitution
    dgbtrs( &lp_trans, &lp_ncols, &lp_wlower, &lp_wupper, &lp_one,
            facmatdata, &lp_ldabf, pivots_, lp_sol, &lp_ncols, &lp_info );

    // Process status flag
    if ( lp_info != 0 ) {
      EXCEPTION( "DGBTRS reports invalid input parameter" );
    }
    else if ( logging == true ) {
      (*cla) << " LAPACK_LU: Solution went fine" << std::endl;
    }

    // ==============================
    //   Refine solution (optional)
    // ==============================
    bool refineSol = true;
      
    if(xml_) {
      xml_->GetValue("refineSol", refineSol, ParamNode::INSERT );
    }
    
    if ( refineSol ) {

      // Prepare some parameters
      double lp_ferr = 0;
      double lp_berr = 0;

      // If not yet done allocate workspaces
      if ( workspaceF77REAL8_ == NULL ) {
	workspaceF77REAL8_ = new double[3 * lp_ncols];
      }
      if ( workspaceInt_ == NULL ) {
	workspaceInt_ = new int[lp_ncols];
      } 
      if ( workspaceF77REAL8_ == NULL || workspaceInt_ == NULL ) {
        EXCEPTION( "Memory allocation for iterative refinement step failed" );
      }

      // Perform iterative refinement
      dgbrfs( &lp_trans, &lp_ncols, &lp_wlower, &lp_wupper, &lp_one,
              matdata, &lp_ldab, facmatdata, &lp_ldabf, pivots_, lp_rhs,
              &lp_ncols, lp_sol, &lp_ncols, &lp_ferr, &lp_berr,
              workspaceF77REAL8_, workspaceInt_, &lp_info );

      // Process results
      if ( lp_info != 0 ) {
        EXCEPTION( "DBGRFS reports invalid input parameter" );
      }
      else if ( logging == true ) {
	(*cla) << " LAPACK_LU: Iterative refinement suceeded\n"
	       << " DGBRFS reports:"
	       << "\n Estimated forward  error bound = " << lp_ferr
	       << "\n Estimated backward error bound = " << lp_berr
	       << std::endl;
      }
    }

    // ============================================
    //   Check whether solution must be re-scaled
    // ============================================

    // If column scaling occured, the solution we now have is the one
    // for A * diag(C) * z = y. Here C is the vector of column scaling
    // factors. We, however, need A * x = y, which can be computed
    // by the identity x = diag(C) * z.
    if ( scalingType_ == 'C' || scalingType_ == 'B' ) {
      for ( Integer i = 0; i < lp_ncols; i++ ) {
	lp_sol[i] *= col_scalings_[i];
      }
    }

    // ===================
    //   Generate Report
    // ===================

    // Now this currently is of dubious value, since the two things queried
    // from olasReport are actually meaningless in the context of a direct
    // solver. Nevertheless we supply some values for consistency

    PtrParamNode out = infoNode_->Get(ParamNode::PN_PROCESS)->Get("solver", ParamNode::APPEND);
    out->Get("numIter")->SetValue(-1);
    out->Get("finalNorm")->SetValue(-1.0);

  }


  // ************************************************
  //   Solution method (complex & double precision)
  // ************************************************
  void Lapack_LU::SolveF77COMPLEX16( const BaseVector &rhs, BaseVector &sol,
				     const BaseMatrix *mat ) {


    // Are we expected to be verbose?
    bool logging = false;

    // Some variables for LAPACK
    char lp_trans = 'N';
    int lp_info = 0;
    int lp_one = 1;

    // Downcast matrices
    LapackGBMatrix<std::complex<double>,std::complex<double> > *facmat =
      dynamic_cast<LapackGBMatrix<std::complex<double>,std::complex<double> > *>
      (facmat_);
    const LapackGBMatrix<std::complex<double>,std::complex<double> > *mymat =
      dynamic_cast<const LapackGBMatrix<std::complex<double>,std::complex<double> > *>
      (mat);

    // This should not have failed, but better test
    if ( facmat == NULL ) {
      EXCEPTION( "Downcast of facmat_ failed!" );
    }
    if ( facmat == NULL ) {
      EXCEPTION( "Downcast of system matrix failed!" );
    }

    // Obtain data pointers
    std::complex<double> *facmatdata = facmat->GetDataPointer0();
    const std::complex<double> *matdata = mymat->GetDataPointer0();

    // Obtain some matrix information
    int lp_ncols  = (int)facmat->GetNumCols();
    int lp_wlower = (int)facmat->GetLowerBandwidth();
    int lp_wupper = (int)facmat->GetUpperBandwidth();
    int lp_ldab   = lp_wlower + lp_wupper + 1;
    int lp_ldabf  = 2 * lp_wlower + lp_wupper + 1;

    // Downcast vectors and get data pointers
    //
    // NOTE: For the moment we generate arrays of type F77COMPLEX16
    //       from the data arrays contained in the sol and rhs vectors
    //       We must copy sol back afterwards
    std::complex<double> *lp_rhs;
    std::complex<double> *lp_sol;
    NEWARRAY( lp_rhs, std::complex<double>, lp_ncols );
    NEWARRAY( lp_sol, std::complex<double>, lp_ncols );
    Complex auxVal1;
    std::complex<double> auxVal2;

    for ( int count = 1; count <= lp_ncols; count++ ) {
      sol.GetEntry( count, auxVal1 );
      lp_sol[count] = auxVal1;

      rhs.GetEntry( count, auxVal1 );
      lp_rhs[count] = auxVal1;
    }

    lp_rhs += 1;
    lp_sol += 1;

    // ====================
    //   Compute solution
    // ====================
   
    // Copy rhs entries into solution vector

    // If row scaling occured, the system matrix was replaced by 
    // diag(R) * A, where R is the vector of row scaling factors.
    // We will get the correct solution, if we transform also the
    // right hand side and solve the system diag(R) * A * x = diag(R) * y.
    if ( scalingType_ == 'R' || scalingType_ == 'B' ) {
      for ( Integer i = 0; i < lp_ncols; i++ ) {
        std::cerr << "R(" << i << ") = " << row_scalings_[i] << std::endl;
        lp_sol[i] = row_scalings_[i] * lp_rhs[i];
      }
    }

    // If now scaling was performed we simply copy the vector
    else {
      for ( Integer i = 0; i < lp_ncols; i++ ) {
	lp_sol[i] = lp_rhs[i];
      }
    }

    // Perform backward/forward substitution
    zgbtrs( &lp_trans, &lp_ncols, &lp_wlower, &lp_wupper, &lp_one,
            facmatdata, &lp_ldabf, pivots_, lp_sol, &lp_ncols, &lp_info );

    // Process status flag
    if ( lp_info != 0 ) {
      EXCEPTION( "DGBTRS reports invalid input parameter" );
    }
    else if ( logging == true ) {
      (*cla) << " LAPACK_LU: Solution went fine" << std::endl;
    }

    // ==============================
    //   Refine solution (optional)
    // ==============================
    bool refineSol = true;
      
    if(xml_) {
      xml_->GetValue("refineSol", refineSol, ParamNode::INSERT);
    }
    
    if ( refineSol ) {

      // Prepare some parameters
      double lp_ferr = 0;
      double lp_berr = 0;

      // If not yet done allocate workspaces
      if ( workspaceF77COMPLEX16_ == NULL ) {
	workspaceF77COMPLEX16_ = new std::complex<double>[3 * lp_ncols];
      }
      if ( workspaceF77REAL8_ == NULL ) {
	workspaceF77REAL8_ = new double[lp_ncols];
      } 
      if ( workspaceF77COMPLEX16_ == NULL || workspaceF77REAL8_ == NULL ) {
        EXCEPTION( "Memory allocation for iterative refinement step failed" );
      }

      // Perform iterative refinement
      zgbrfs( &lp_trans, &lp_ncols, &lp_wlower, &lp_wupper, &lp_one,
              matdata, &lp_ldab, facmatdata, &lp_ldabf, pivots_, lp_rhs,
              &lp_ncols, lp_sol, &lp_ncols, &lp_ferr, &lp_berr,
              workspaceF77COMPLEX16_, workspaceF77REAL8_, &lp_info );

      // Process results
      if ( lp_info != 0 ) {
        EXCEPTION( "ZBGRFS reports invalid input parameter" );
      }
      else if ( logging == true ) {
	(*cla) << " Iterative refinement suceeded!\n ZGBRFS reports:"
	       << "\n Estimated forward  error bound = " << lp_ferr
	       << "\n Estimated backward error bound = " << lp_berr
	       << std::endl;
      }
    }

    // =============================
    //   Copy solution into vector
    // =============================
    for ( int count = 1; count <= lp_ncols; count++ ) {
      auxVal2 = lp_sol[count-1];
      sol.SetEntry( count, auxVal2 );
    }

    // ============================================
    //   Check whether solution must be re-scaled
    // ============================================

    // If column scaling occured, the solution we now have is the one
    // for A * diag(C) * z = y. Here C is the vector of column scaling
    // factors. We, however, need A * x = y, which can be computed
    // by the identity x = diag(C) * z.
    if ( scalingType_ == 'C' || scalingType_ == 'B' ) {
      for ( Integer i = 0; i < lp_ncols; i++ ) {
#if 0
	lp_sol[i].real() *= col_scalings_[i];
	lp_sol[i].imag() *= col_scalings_[i];
#endif
        lp_sol[i] *= col_scalings_[i];
      }
    }

    // ===================
    //   Generate Report
    // ===================

    // Now this currently is of dubious value, since the two things queried
    // from olasReport are actually meaningless in the context of a direct
    // solver. Nevertheless we supply some values for consistency
    /* FIXME this pollutes the xml file and needs to be updated to provide more information
    PtrParamNode out = infoNode_->Get(ParamNode::PN_PROCESS)->Get("solver", ParamNode::APPEND);
    out->Get("numIter")->SetValue(-1);
    out->Get("finalNorm")->SetValue(-1.0);
    */

  }

}
