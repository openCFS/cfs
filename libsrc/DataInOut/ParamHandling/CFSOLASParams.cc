#ifdef USE_OLAS

#include <string.h>

#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "CFSOLASParams.hh"

namespace CoupledField {

  // *************
  //   SetParams
  // *************
  void CFSOLASParams::SetParams( std::string pdename,
				 BaseParamHandler *cfs,
				 OLAS_Params *olas, bool overrideExpert ) {

    ENTER_FCN( "CFSOLASParams::SetParams" );

    // First determine the type of solver for this PDE
    std::string sTypeString;
    SolverType sType;
    cfs->Get( "type", sTypeString, pdename, "solver" );
    if ( sTypeString == "CG" ) {
      sType = CG;
    }
    else if ( sTypeString == "hyprePCG" ) {
      sType = HYPRE_PCG;
    }
    else if ( sTypeString == "hypreGMRES" ) {
      sType = HYPRE_GMRES;
    }
    else if ( sTypeString == "hypreBICGSTAB" ) {
      sType = HYPRE_BICGSTAB;
    }
    else if ( sTypeString == "lapackLU" ) {
      sType = LAPACK_LU;
    }
    else {
      std::string errmsg = "Solver " + sTypeString;
      errmsg += " not supported yet.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Now determine the type of preconditioner for this PDE
    std::string pTypeString;
    PrecondType pType;
    cfs->Get( "precond", pTypeString, pdename, "solver" );
    if ( pTypeString == "noPrecond" ) {
      pType = NOPRECOND;
    }
    else if ( pTypeString == "Id" ) {
      pType = ID;
    }
    else if ( pTypeString == "Jacobi" ) {
      pType = ID;
    }
    else if ( pTypeString == "SSOR" ) {
      pType = SSOR;
    }
    else if ( pTypeString == "MG" ) {
      pType = MG;
    }
    else if ( pTypeString == "hypreJACOBI" ) {
      pType = HYPRE_JACOBI;
    }
    else if ( pTypeString == "hypreSPAI" ) {
      pType = HYPRE_SPAI;
    }
    else if ( pTypeString == "hypreILU" ) {
      pType = HYPRE_ILU;
    }
    else if ( pTypeString == "hypreAMG" ) {
      pType = HYPRE_AMG;
    }
    else {
      std::string errmsg = "Preconditioner " + pTypeString;
      errmsg += " not supported yet.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // Next determine the matrix type for this PDE
    std::string mMatString;
    MatrixStorageType mType;
    cfs->Get( "storage", mMatString, pdename, "matrix" );
    if ( mMatString == "sparseSym" ) {
      mType = SPARSE_SYM;
    }
    else if ( mMatString == "sparseNonSym" ) {
      mType = SPARSE_NONSYM;
    }
    else if ( mMatString == "skylineSym" ) {
      mType = SKYLINE_SYM;
    }
    else if ( mMatString == "skylineNonSym" ) {
      mType = SKYLINE_NONSYM;
    }
    else if ( mMatString == "hypreMatrix" ) {
      mType = HYPRE_MATRIX;
    }
    else if ( mMatString == "lapackGBMatrix" ) {
      mType = LAPACK_GBMATRIX;
    }
    else if ( mMatString == "lapackPBMatrix" ) {
      mType = LAPACK_PBMATRIX;
    }
    else if ( mMatString == "expertsChoice" ) {
      if ( overrideExpert ) {
	std::string errmsg = "You cannot specify expertsChoice as storage ";
	errmsg += "'type and set overrideExpert! This would leave the storage";
	errmsg += " type undefined!";
	Info->Error( errmsg, __FILE__, __LINE__ );
      }
      mType = NOSTORAGETYPE;
    }
    else {
      std::string errmsg = "Matrix storage type '" + mMatString;
      errmsg += "' not supported yet.";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }


    std::string eMatString;
    MatrixEntryType eType;
    cfs->Get( "entry", eMatString, pdename, "matrix" );
    if ( eMatString == "double" ) {
      eType = DOUBLE;
    }
    else if ( eMatString == "complex" ) {
      eType = COMPLEX;
    }

    // Let expert module modify the settings
    if ( !overrideExpert ) {
      CFSOLASParams::Expert( sType, pType, mType, eType );
    }

    // Insert information into OLAS_Params object
    olas->SetValue( "Solver", sType );
    olas->SetValue( "Precond", pType );
    olas->SetValue( "MatrixStructureType", STDMATRIX );
    olas->SetValue( "MatrixStorageType", mType );
    olas->SetValue( "MatrixEntryType"  , eType );

    // Set special parameters for solver and preconditioner
    CFSOLASParams::SetSolverParams( pdename, cfs, olas, sType );
    CFSOLASParams::SetPrecondParams( pdename, cfs, olas, pType );

    // For debugging
    // olas->ShowPool( OLAS_Params::INT_POOL     , std::cerr );
    // olas->ShowPool( OLAS_Params::DOUBLE_POOL  , std::cerr );
    // olas->ShowPool( OLAS_Params::BOOLEAN_POOL , std::cerr );
    // olas->ShowPool( OLAS_Params::STRING_POOL  , std::cerr );
    // olas->ShowPool( OLAS_Params::ENUM_POOL    , std::cerr );

    // Output Matrix
    StdVector<std::string> doExport;
    cfs->GetList( "file", doExport, pdename, "exportLinSys" );
    if ( doExport.GetSize() == 1 ) {
      olas->SetValue( "exportLinSys", true );
      olas->SetValue( "exportLinSysBaseName", doExport[0] );
    }

  }


  // *******************
  //   SetSolverParams
  // *******************
  void CFSOLASParams::SetSolverParams( std::string pdename,
				       BaseParamHandler *cfs,
				       OLAS_Params *olas,
				       SolverType sType ) {

    ENTER_FCN( "CFSOLASParams::SetSolverParams" );

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    StdVector<std::string> list;

    switch( sType ) {

    case CG:
      cfs->GetList( "tol", list, pdename, "PCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CG_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hyprePCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CG_maxIter", atoi(list[0].c_str()) );
      }
      break;

    case HYPRE_PCG:
      cfs->GetList( "tol", list, pdename, "hyprePCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREPCG_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hyprePCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREPCG_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "hyprePCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREPCG_logging", atoi(list[0].c_str()) );
      }
      cfs->GetList( "printLevel", list, pdename, "hyprePCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREPCG_printLevel", atoi(list[0].c_str()) );
      }
      break;

    case HYPRE_GMRES:
      cfs->GetList( "tol", list, pdename, "hypreGMRES" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREGMRES_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hypreGMRES" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREGMRES_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "maxKrylovDim", list, pdename, "hypreGMRES" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREGMRES_maxKrylovDim", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "hypreGMRES" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREGMRES_logging", atoi(list[0].c_str()) );
      }
      cfs->GetList( "printLevel", list, pdename, "hypreGMRES" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREGMRES_printLevel", atoi(list[0].c_str()) );
      }
      break;

    case HYPRE_BICGSTAB:
      cfs->GetList( "tol", list, pdename, "hypreBICGSTAB" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "hypreBICGSTAB" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "hypreBICGSTAB" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_logging", atoi(list[0].c_str()) );
      }
      cfs->GetList( "printLevel", list, pdename, "hypreBICGSTAB" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREBICGSTAB_printLevel", atoi(list[0].c_str()) );
      }
      break;

    }

  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::SetPrecondParams( std::string pdename,
					BaseParamHandler *cfs,
					OLAS_Params *olas,
					PrecondType pType ) {

    ENTER_FCN( "CFSOLASParams::SetPrecondParams" );

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    StdVector<std::string> list;

    switch( pType ) {

    case HYPRE_ILU:
      cfs->GetList( "level", list, pdename, "hypreILU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "EUCLID_level",  atoi(list[0].c_str()) );
      }
      cfs->GetList( "stats", list, pdename, "hypreILU" );
      if( list.GetSize() == 1 ) {
	bool stats = (list[0] == "yes");
	olas->SetValue( "EUCLID_stats", stats );
      }
      cfs->GetList( "memory", list, pdename, "hypreILU" );
      if( list.GetSize() == 1 ) {
	Integer memory = (list[0] == "yes");
	olas->SetValue( "EUCLID_memory", memory );
      }
      break;

    case HYPRE_SPAI:
      cfs->GetList( "thresh", list, pdename, "hypreSPAI" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARASAILS_thresh", atof(list[0].c_str()) );
      }
      cfs->GetList( "levels", list, pdename, "hypreSPAI" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARASAILS_levels", atoi(list[0].c_str()) );
      }
      cfs->GetList( "filter", list, pdename, "hypreSPAI" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARASAILS_filter", atof(list[0].c_str()) );
      }
      cfs->GetList( "symmetry", list, pdename, "hypreSPAI" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARASAILS_symmetry", atoi(list[0].c_str()) );
      }
      cfs->GetList( "loadBalance", list, pdename, "hypreSPAI" );
      if( list.GetSize() == 1 ) {
	Integer balance = (list[0] == "yes");
	olas->SetValue( "PARASAILS_loadBalance", balance );
      }
      cfs->GetList( "reuse", list, pdename, "hypreSPAI" );
      if( list.GetSize() == 1 ) {
	Integer balance = (list[0] == "yes");
	olas->SetValue( "PARASAILS_reuse", balance );
      }
      break;

    case HYPRE_AMG:
      cfs->GetList( "maxLevels", list, pdename, "hypreAMG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "BOOMERAMG_maxLevels", atoi(list[0].c_str()) );
      }
      cfs->GetList( "alpha", list, pdename, "hypreAMG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "BOOMERAMG_alpha", atof(list[0].c_str()) );
      }
      cfs->GetList( "numSweeps", list, pdename, "hypreAMG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "BOOMERAMG_numSweeps", atoi(list[0].c_str()) );
      }
      cfs->GetList( "cycleType", list, pdename, "hypreAMG" );
      if( list.GetSize() == 1 ) {
	Integer cycleType = list[0] == "W-cycle" ? 2 : 1;
	olas->SetValue( "BOOMERAMG_cycleType", cycleType );
      }
      break;
    }
  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::Expert( SolverType &sType, PrecondType &pType,
			      MatrixStorageType &mType,
			      MatrixEntryType &eType ) {

    ENTER_FCN( "CFSOLASParams::Expert" );

    std::string warn;

    // ======================
    //  Precondtioner stuff
    // ======================

    // For direct solver we need no preconditioner
    if ( sType == DIRECT || sType == LAPACK_LU && pType != NOPRECOND ) {
      warn = "Expert: Re-setting preconditioner type to 'NOPRECOND'";
      Info->Warning( warn );
      pType = NOPRECOND;
    }

    // Hypre solvers only work together with hypre preconditioners
    if ( sType == HYPRE_PCG || sType == HYPRE_GMRES || sType == HYPRE_BICGSTAB
	 && !( pType == NOPRECOND  || pType == HYPRE_AMG ||
	       pType == HYPRE_SPAI || pType == HYPRE_ILU ) ) {
      warn = "Expert: Re-setting preconditioner type to 'NOPRECOND'";
      Info->Warning( warn );
      pType = NOPRECOND;
    }

    // ===============
    //  Matrix stuff
    // ===============

    // Lapack solvers want their own matrix format
    if ( sType == LAPACK_LU ) {
      if ( mType != LAPACK_GBMATRIX ) {
	if ( mType != NOSTORAGETYPE ) {
	  warn = "Expert: Re-setting matrix storage type to 'LAPACK_GBMATRIX'";
	  Info->Warning( warn );
	}
	mType = LAPACK_GBMATRIX;
      }

      // The following block is useless, as long as OLAS does not really
      // use F77xxx entry type specifiers

      // if ( eType == DOUBLE ) {
      //   if ( sizeof(Double) == sizeof(float) ) {
      //     eType = F77FLOAT;
      //     warn = "Expert: Re-setting matrix entry type to F77FLOAT";
      //     Info->Warning( warn );
      //   }
      //   else {
      //     eType = F77DOUBLE;
      //     warn = "Expert: Re-setting matrix entry type to F77DOUBLE";
      //     Info->Warning( warn );
      //   }
      // }
      // else if ( eType == COMPLEX ) {
      //   if ( sizeof(Complex) == sizeof(std::complex<float>) ) {
      //     eType = F77COMPLEX8;
      //     warn = "Expert: Re-setting matrix entry type to F77COMPLEX8";
      //     Info->Warning( warn );
      //   }
      //   else {
      //     eType = F77COMPLEX16;
      //     warn = "Expert: Re-setting matrix entry type to F77COMPLEX16";
      //     Info->Warning( warn );
      //   }
      // }
    }

    // Hypre solvers want their own matrix format
    if ( sType == HYPRE_PCG || sType == HYPRE_GMRES ||
	 sType == HYPRE_BICGSTAB ) {
      if ( mType != HYPRE_MATRIX ) {
	if ( mType != NOSTORAGETYPE ) {
	  warn = "Expert: Re-setting matrix storage type to 'HYPRE_MATRIX'";
	  Info->Warning( warn );
	}
	mType = HYPRE_MATRIX;
      }
      if ( eType != DOUBLE ) {
	warn = "Expert: Re-setting matrix entry type to DOUBLE";
	Info->Warning( warn );
	eType = DOUBLE;
      }
    }

    // If no storage type was yet assigned use plain CRS
    if ( mType == NOSTORAGETYPE ) {
      mType = SPARSE_NONSYM;
    }

  }
}
#endif
