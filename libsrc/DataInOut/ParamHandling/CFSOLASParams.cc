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
    if ( sTypeString == "expertsChoice" ) {
      if ( overrideExpert ) {
	(*error) << "You cannot specify expertsChoice as solver type "
		 << "and at the same time specify overrideExpert! "
		 << "This would leave the solver type undefined!";
	Error( __FILE__, __LINE__ );
      }
      sTypeString = "no solver";
    }
    OLAS::String2Enum( sTypeString, sType );

    // Now determine the type of preconditioner for this PDE
    std::string pTypeString;
    PrecondType pType;
    cfs->Get( "precond", pTypeString, pdename, "solver" );
    OLAS::String2Enum( pTypeString, pType );

    // Next determine the matrix type for this PDE
    std::string mMatString;
    MatrixStorageType mType;
    cfs->Get( "storage", mMatString, pdename, "matrix" );
    if ( mMatString == "expertsChoice" ) {
      if ( overrideExpert ) {
	(*error) << "You cannot specify expertsChoice as storage type "
		 << "and at the same time specify overrideExpert! "
		 << "This would leave the storage type undefined!";
	Error( __FILE__, __LINE__ );
      }
      mMatString = "noStorageType";
    }
    OLAS::String2Enum( mMatString, mType );

    // Determine matrix entry type
    std::string eMatString;
    MatrixEntryType eType;
    cfs->Get( "entry", eMatString, pdename, "matrix" );
    if ( eMatString == "double" ) {
      eType = DOUBLE;
    }
    else if ( eMatString == "complex" ) {
      eType = COMPLEX;
    }

    // Type of re-odering
    std::string orderString;
    ReorderingType orderType;
    cfs->Get( "reordering", orderString, pdename, "matrix" );
    OLAS::String2Enum( orderString, orderType );

    // Let expert module modify the settings
    if ( !overrideExpert ) {
      CFSOLASParams::Expert( cfs, pdename, sType, pType, mType, eType,
			     orderType );
    }

    // Insert information into OLAS_Params object
    olas->SetValue( "Solver"             , sType );
    olas->SetValue( "Precond"            , pType );
    olas->SetValue( "MatrixStructureType", STDMATRIX );
    olas->SetValue( "MatrixStorageType"  , mType );
    olas->SetValue( "MatrixEntryType"    , eType );
    olas->SetValue( "GRAPH_reordering"   , orderType );

    // Set special parameters for solver and preconditioner
    CFSOLASParams::SetSolverParams( pdename, cfs, olas, sType );
    CFSOLASParams::SetPrecondParams( pdename, cfs, olas, pType );

    // For debugging (please do not delete)
    // olas->ShowPool( OLAS_Params::INT_POOL     , std::cerr );
    // olas->ShowPool( OLAS_Params::DOUBLE_POOL  , std::cerr );
    // olas->ShowPool( OLAS_Params::BOOLEAN_POOL , std::cerr );
    // olas->ShowPool( OLAS_Params::STRING_POOL  , std::cerr );
    // olas->ShowPool( OLAS_Params::ENUM_POOL    , std::cerr );

    // Output Matrix
    StdVector<std::string> doExport;
    cfs->GetList( "baseName", doExport, pdename, "exportLinSys" );
    if ( doExport.GetSize() == 1 ) {
      olas->SetValue( "exportLinSys", true );
      olas->SetValue( "exportLinSysBaseName", doExport[0] );
    }
    else {
      olas->SetValue( "exportLinSys", false );
    }

    // For hypre solvers
    olas->SetValue( "AccountForPenalty", false );

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
      cfs->GetList( "tol", list, pdename, "cg" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CG_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "cg" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CG_maxIter", atoi(list[0].c_str()) );
      }
      break;

    case GMRES:
      cfs->GetList( "tol", list, pdename, "gmres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "GMRES_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "gmres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "GMRES_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "maxKrylovDim", list, pdename, "gmres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "GMRES_maxKrylovDim", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "gmres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "GMRES_logging", (list[0] == "yes") );
      }
      break;

    case MINRES:
      cfs->GetList( "tol", list, pdename, "minres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "MINRES_epsilon", atof(list[0].c_str()) );
      }
      cfs->GetList( "maxIter", list, pdename, "minres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "MINRES_maxIter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "minres" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "MINRES_logging", (list[0] == "yes") );
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
      cfs->GetList( "twoNorm", list, pdename, "hyprePCG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "HYPREPCG_twoNorm", (list[0] == "yes") );
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

    case LAPACK_LU:
      cfs->GetList( "tryScaling", list, pdename, "lapackLU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "LAPACKLU_tryScaling", (list[0] == "yes") );
      }
      cfs->GetList( "refineSol", list, pdename, "lapackLU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "LAPACKLU_refineSol", (list[0] == "yes") );
      }
      cfs->GetList( "logging", list, pdename, "lapackLU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "LAPACKLU_logging", (list[0] == "yes") );
      }
      break;

    case LAPACK_LL:
      cfs->GetList( "logging", list, pdename, "lapackLL" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "LAPACKLL_logging", (list[0] == "yes") );
      }
      break;

    case LU_SOLVER:
      cfs->GetList( "logging", list, pdename, "directLU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "LUSOLVER_logging", (list[0] == "yes") );
      }
      cfs->GetList( "saveFacFile", list, pdename, "directLU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CROUT_saveFacToFile", true );
	olas->SetValue( "CROUT_facFileName", list[0] );
      }
      break;

    case LDL_SOLVER:
      cfs->GetList( "logging", list, pdename, "directLDL" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "LDLSOLVER_logging", (list[0] == "yes") );
      }
      cfs->GetList( "saveFacFile", list, pdename, "directLDL" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CROUT_saveFacToFile", true );
	olas->SetValue( "CROUT_facFileName", list[0] );
      }
      break;

    case PARDISO:
      cfs->GetList( "posDef", list, pdename, "pardiso" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARDISO_posDef", (list[0] == "yes") );
      }
      cfs->GetList( "hermitean", list, pdename, "pardiso" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARDISO_hermitean", (list[0] == "yes") );
      }
      cfs->GetList( "symStruct", list, pdename, "pardiso" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARDISO_symStructure", (list[0] == "yes") );
      }
      cfs->GetList( "ndOrdering", list, pdename, "pardiso" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARDISO_ndOdering", (list[0] == "yes") );
      }
      cfs->GetList( "logging", list, pdename, "pardiso" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "PARDISO_logging", (list[0] == "yes") );
      }
      break;

      // If this point is reached, it indicates that something is broken.
      // Probably not all solvers that SetParams allows are yet implemented
      // here?
    default:
      Error( "Internal error. Maybe a missing implementation of a case?",
	     __FILE__, __LINE__ );
      break;

    }

    // Now set the stopping rule
    std::string stopCrit;
    cfs->Get( "type", stopCrit, pdename, "stoppingRule" );
    StopCritType stopRule;
    OLAS::String2Enum( stopCrit, stopRule );
    olas->SetValue( "StoppingCriterion", stopRule );

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

      // The easiest cases. No parameters exist.
    case ID:
    case NOPRECOND:
    case JACOBI:
    case ILU0:
    case HYPRE_JACOBI:
      break;

    case HYPRE_ILU:
      cfs->GetList( "level", list, pdename, "hypreILU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "EUCLID_level",  atoi(list[0].c_str()) );
      }
      cfs->GetList( "stats", list, pdename, "hypreILU" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "EUCLID_stats", (list[0] == "yes") );
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

    case MG:
      cfs->GetList( "maxCoarseDepend", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_MaxCoarseDependency", atoi(list[0].c_str()) );
      }
      cfs->GetList( "minSystemSize", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_MinSystemSize", atoi(list[0].c_str()) );
      }
      cfs->GetList( "numPreSmooth", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_NumPreSmoothing", atoi(list[0].c_str()) );
      }
      cfs->GetList( "numPostSmooth", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_NumPostSmoothing", atoi(list[0].c_str()) );
      }
      cfs->GetList( "cycleParam", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_CycleParameter", atoi(list[0].c_str()) );
      }
      cfs->GetList( "alpha", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_CoarseningAlpha", atof(list[0].c_str()) );
      }
      cfs->GetList( "strongDiagRatio", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_StrongDiagRatio", atof(list[0].c_str()) );
      }
      cfs->GetList( "forceFineRatio", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_ForceFineRatio", atof(list[0].c_str()) );
      }
      cfs->GetList( "directSolver", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	if ( list[0] == "lapackLU" ) {
	  olas->SetValue( "AMG_DirectSolver", LAPACK_LU );
	}
	else {
	  olas->SetValue( "AMG_DirectSolver", NOSOLVER );
	}
      }
      cfs->GetList( "logging", list, pdename, "MG" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "AMG_logging", (list[0] == "yes") );
      }
      break;

    case ILUK:
      cfs->GetList( "level", list, pdename, "ILUK" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILUK_level", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "ILUK" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILUK_logging", (list[0] == "yes") );
      }
      cfs->GetList( "saveFacFile", list, pdename, "ILUK" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "CROUT_saveFacToFile", true );
	olas->SetValue( "CROUT_facFileName", list[0] );
      }
      break;

    case ILDLK:

      // Set sub-type
      olas->SetValue( "ILDLPRECOND_subType", ILDLK );

      // Now parameters
      cfs->GetList( "level", list, pdename, "ILDLK" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLPRECOND_level", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "ILDLK" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLKPRECOND_logging", (list[0] == "yes") );
      }
      cfs->GetList( "saveFacFile", list, pdename, "ILDLK" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
	olas->SetValue( "ILDLPRECOND_facFileName", list[0] );
        cfs->GetList( "savePatternOnly", list, pdename, "ILDLK" );
        if( list.GetSize() == 1 ) {
          olas->SetValue( "ILDLPRECOND_facPatternOnly", (list[0] == "yes") );
        }
      }
      break;

    case ILDLTP:

      // Set sub-type
      olas->SetValue( "ILDLPRECOND_subType", ILDLTP );

      // Now parameters
      cfs->GetList( "threshold", list, pdename, "ILDLTP" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLPRECOND_tau", atof(list[0].c_str()) );
      }
      cfs->GetList( "fillVal", list, pdename, "ILDLTP" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLPRECOND_fillVal", atoi(list[0].c_str()) );
      }
      cfs->GetList( "logging", list, pdename, "ILDLTP" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLTPPRECOND_logging", (list[0] == "yes") );
      }
      cfs->GetList( "saveFacFile", list, pdename, "ILDLTP" );
      if( list.GetSize() == 1 ) {
	olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
	olas->SetValue( "ILDLPRECOND_facFileName", list[0] );
        cfs->GetList( "savePatternOnly", list, pdename, "ILDLTP" );
        if( list.GetSize() == 1 ) {
          olas->SetValue( "ILDLPRECOND_facPatternOnly", (list[0] == "yes") );
        }
      }
      break;

      // If this point is reached, it indicates that something is broken.
      // Probably not all preconditioners that SetParams allows are yet
      // implemented here?
    default:
      Error( "Internal error. Maybe a missing implementation of a case?",
	     __FILE__, __LINE__ );
      break;
    }
  }


  // *******************
  //   Expert's Choice
  // *******************
  void CFSOLASParams::Expert( BaseParamHandler *cfs,
			      std::string pdename,
			      SolverType &sType,
			      PrecondType &pType,
			      MatrixStorageType &mType,
			      MatrixEntryType &eType,
			      ReorderingType &rType ) {

    ENTER_FCN( "CFSOLASParams::Expert" );

    std::string warn;

    // ==============
    //  Solver stuff
    // ==============

    // If no solver was specified use a direct one
    // Currently always use LU until LDL is available
    if ( sType == NOSOLVER ) {
      sType = LU_SOLVER;
    }


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
    if ( (sType == HYPRE_PCG || sType == HYPRE_GMRES ||
	  sType == HYPRE_BICGSTAB)

	 && !( pType == NOPRECOND    ||
	       pType == ID           ||
	       pType == HYPRE_JACOBI ||
	       pType == HYPRE_SPAI   ||
	       pType == HYPRE_ILU    ||
	       pType == HYPRE_AMG ) ) {

      warn = "Expert: Re-setting preconditioner type to 'ID'";
      Info->Warning( warn );
      pType = ID;
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
	else {
	  Info->PrintF( pdename,
			"Expert: Using LAPACK_GBMATRIX as storage type\n" );
	}
	mType = LAPACK_GBMATRIX;
      }
    }

    // The direct solver LU_SOLVER expects a CRS matrix
    if ( sType == LU_SOLVER ) {
      if ( mType != SPARSE_NONSYM ) {
	if ( mType != NOSTORAGETYPE ) {
	  (*warning) << "Expert: Changing matrix storage type from "
		     << Enum2String( mType ) << " to SPARSE_NONSYM for "
		     << "direct solver";
	  Warning( __FILE__, __LINE__ );
	}
	else {
	  Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
			"type for direct solver\n" );
	}
	mType = SPARSE_NONSYM;
      }
    }

    // Hypre solvers want their own matrix format
    if ( sType == HYPRE_PCG || sType == HYPRE_GMRES ||
	 sType == HYPRE_BICGSTAB ) {
      if ( mType != HYPRE_MATRIX ) {
	if ( mType != NOSTORAGETYPE ) {
	  warn = "Expert: Re-setting matrix storage type to 'HYPRE_MATRIX'\n";
	  Info->Warning( warn );
	}
	else {
	  Info->PrintF( pdename, "Expert: Using HYPRE_MATRIX as storage "
			"type\n" );
	}
	mType = HYPRE_MATRIX;
      }
      if ( eType != DOUBLE ) {
	warn = "Expert: Re-setting matrix entry type to DOUBLE\n";
	Info->Warning( warn );
	eType = DOUBLE;
      }
    }

    // If no storage type was yet assigned use plain CRS
    if ( mType == NOSTORAGETYPE ) {
      mType = SPARSE_NONSYM;
      Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage type\n" );
    }

    // In case of a harmonic analysis or parameter identification
    // we assume complex matrix entries by default
    std::string analysis;
    cfs->Get( "type", analysis, "analysis" );
    if ( analysis == "harmonic"
	 && eType != COMPLEX ) {
      eType = COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
		    "for HARMONIC analysis\n" );
    }

    if ( analysis == "paramIdent" && eType != COMPLEX ) {
      eType = COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
		    "for parameter identification\n" );
    }

    // ============
    //  Reordering
    // ============
    if ( sType == LAPACK_LU || sType == LU_SOLVER || sType == LDL_SOLVER ||
	 pType == ILU0 || pType == ILUK ) {
      if ( rType == NOREORDERING ) {
	Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
		      "'SLOAN'" );
	rType = SLOAN;
      }
      else {
	Info->PrintF( pdename, "Expert: Using SLOAN for re-ordering\n" );
	rType = SLOAN;
      }
    }
  }

}
