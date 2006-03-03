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
                                 OLAS_Params *olas, 
                                 AnalysisType analysisType,
                                 bool overrideExpert ) {

    ENTER_FCN( "CFSOLASParams::SetParams" );

    // We need three vectors for our parameter queries
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // First determine approach for handling inhomogeneous Dirichlet
    // boundary conditions
    bool usingPenalty = true;
    {
      std::string aux;
      keyVec  = "linearSystems", "system", "setup", "idbcHandling";
      attrVec = "", "name", "";
      valVec  = "", pdename, "";
      params->Get( keyVec, attrVec, valVec, aux );
      usingPenalty = aux == "penalty" ? true : false;
    }

    // Determine the type of solver for this PDE
    std::string sTypeString;
    OLAS::SolverType sType;

    keyVec  = "linearSystems", "system", "solver", "type";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";
    cfs->Get( keyVec, attrVec, valVec, sTypeString );

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
    OLAS::PrecondType pType;

    keyVec  = "linearSystems", "system", "solver", "precond";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";
    cfs->Get( keyVec, attrVec, valVec, pTypeString );

    OLAS::String2Enum( pTypeString, pType );

    // Now determine, if we are running an eigenfrequency analysis and if
    // yes, get the type of eigenvalue solver
    std::string esTypeString;
    OLAS::EigenSolverType esType = OLAS::NOEIGENSOLVER;
    
    if ( analysisType == EIGENFREQUENCY ) {
      keyVec  = "linearSystems", "system", "eigenSolver", "type";
      attrVec = "", "name", "";
      valVec  = "", pdename, "";
      cfs->Get( keyVec, attrVec, valVec, esTypeString );
      OLAS::String2Enum( esTypeString, esType );
    }

    // Before we try to determine matrix properties, we first should
    // find out, what type of linear system we are using. We do this
    // by trying to determine the symmetry attribute of the sbmMatrix
    // element
    bool sbmSymmetry = false;
    bool stdSystem = true;
    StdVector<std::string> sbmSymmetric;
    keyVec  = "linearSystems", "system", "sbmMatrix", "symmetric";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";
    cfs->GetList( keyVec, attrVec, valVec, sbmSymmetric );

    // If the list is not empty we got a sbmMatrix
    if ( sbmSymmetric.GetSize() == 1 ) {
      sbmSymmetry = sbmSymmetric[0] == "yes";
      stdSystem = false;
    }
    else if ( sbmSymmetric.GetSize() > 1 ) {
      (*error) << "Found multiple 'symmtric' attributes for (multiple?) "
               << "'sbmMatrix' elements! We are screwed!";
      Error( __FILE__, __LINE__ );
    }

    // Set matrix element type string for further queries
    std::string matElement = stdSystem == true ? "matrix" : "sbmMatrix";

    // Determine matrix entry type
    std::string eMatString;
   OLAS:: MatrixEntryType eType;

    keyVec  = "linearSystems", "system", matElement, "entry";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";
    cfs->Get( keyVec, attrVec, valVec, eMatString );

    if ( eMatString == "double" ) {
      eType = OLAS::DOUBLE;
    }
    else if ( eMatString == "complex" ) {
      eType = OLAS::COMPLEX;
    }

    // Following stuff only for required for standard systems
    OLAS::MatrixStorageType mType  = OLAS::NOSTORAGETYPE;
    OLAS::ReorderingType orderType = OLAS::NOREORDERING;
    bool allowChangeOfReordering = false;

    if ( stdSystem == true ) {

      // Next determine the matrix type for this PDE
      std::string mMatString;

      keyVec  = "linearSystems", "system", "matrix", "storage";
      attrVec = "", "name", "";
      valVec  = "", pdename, "";
      cfs->Get( keyVec, attrVec, valVec, mMatString );

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

      // Type of re-odering
      std::string orderString;

      keyVec  = "linearSystems", "system", "matrix", "reordering";
      attrVec = "", "name", "";
      valVec  = "", pdename, "";
      cfs->Get( keyVec, attrVec, valVec, orderString );

      if ( orderString == "expertsChoice" ) {
        if ( overrideExpert ) {
          (*error) << "You cannot specify expertsChoice as re-ordering "
                   << "strategy and at the same time specify overrideExpert! "
                   << "This would leave the re-ordering strategy undefined!";
          Error( __FILE__, __LINE__ );
        }
        orderString = "noReordering";
        allowChangeOfReordering = true;
      }
      OLAS::String2Enum( orderString, orderType );
    }

    // Let expert module modify the settings
    if ( overrideExpert == false && stdSystem == true ) {
      CFSOLASParams::Expert( cfs, pdename, sType, pType, mType, eType,
                             orderType, analysisType,
                             allowChangeOfReordering );
    }

    // Insert information into OLAS_Params object
    olas->SetValue( "Solver"                 , sType              );
    olas->SetValue( "Precond"                , pType              );
    olas->SetValue( "EigenSolver"            , esType             );
    olas->SetValue( "MatrixStructureType"    , OLAS::STDMATRIX    );
    olas->SetValue( "MatrixStorageType"      , mType              );
    olas->SetValue( "MatrixEntryType"        , eType              );
    olas->SetValue( "GRAPH_reordering"       , orderType          );
    olas->SetValue( "UsingPenaltyFormulation", usingPenalty       );
    olas->SetValue( "SystemName"             , pdename            );
    olas->SetValue( "SBM_Symmetry"           , sbmSymmetry        );

    // Set special parameters for solver and preconditioner
    CFSOLASParams::SetSolverParams( pdename, cfs, olas, sType );
    CFSOLASParams::SetPrecondParams( pdename, cfs, olas, pType );
    CFSOLASParams::SetEigenSolverParams( pdename, cfs, olas, esType );

    // For debugging (please do not delete)
    // olas->ShowPool( OLAS_Params::INT_POOL     , std::cerr );
    // olas->ShowPool( OLAS_Params::DOUBLE_POOL  , std::cerr );
    // olas->ShowPool( OLAS_Params::BOOLEAN_POOL , std::cerr );
    // olas->ShowPool( OLAS_Params::STRING_POOL  , std::cerr );
    // olas->ShowPool( OLAS_Params::ENUM_POOL    , std::cerr );

    // Output Matrix
    keyVec  = "linearSystems", "system", "exportLinSys", "baseName";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";
    StdVector<std::string> doExport;
    cfs->GetList( keyVec, attrVec, valVec, doExport );
    if ( doExport.GetSize() == 1 ) {
      olas->SetValue( "exportLinSys", true );
      olas->SetValue( "exportLinSysBaseName", doExport[0] );
    }
    else {
      olas->SetValue( "exportLinSys", false );
    }

    // For hypre solvers
    olas->SetValue( "AccountForPenalty", false );

    // Consistency test
    //
    // Currently only harmonic and paramIdent analysis can deal with 
    // complex material parameters.
    // For static and transient this will lead to a seg fault.
    // So we check this to avoid this unclear problem
    StdVector<std::string> matType;
    std::string analysis;
    cfs->GetList( "type", matType, pdename, "materialDataType" );
    if ( matType.GetSize() > 0 ) {
      if ( matType.GetSize() != 1 ) {
        (*error) << "materialDataType not unique, GetList returned "
                 << matType.GetSize() << " values";
        Error( __FILE__, __LINE__ );
      }
      cfs->Get( "type", analysis, "analysis" );
      if ( matType[0] == "imagMaterialParameter" &&
           analysis != "harmonic" && analysis != "paramIdent" &&
           analysis!="multiHarmonic") {
        (*error) << "XML-file specifies material parameters with imaginary "
                 << "part for an analysis\n of type '" << analysis << "'.\n\n"
                 << " Complex parameters are currently only implemented for "
                 << "a 'harmonic',\n 'multiHarmonic', or 'paramIdent' "
                 << "analysis, however.";
        Error( __FILE__, __LINE__ );
      }
    }
  }


  // *******************
  //   SetSolverParams
  // *******************
  void CFSOLASParams::SetSolverParams( std::string pdename,
                                       BaseParamHandler *cfs,
                                       OLAS::OLAS_Params *olas,
                                       OLAS::SolverType sType ) {

    ENTER_FCN( "CFSOLASParams::SetSolverParams" );

    // We need three vectors for our parameter queries
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // Basic part of query will always be the same
    keyVec  = "linearSystems", "system", "", "";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    StdVector<std::string> list;

    switch( sType ) {

    case OLAS::CG:
      keyVec[2] = "cg";

      keyVec[3] = "tol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "CG_epsilon", atof(list[0].c_str()) );
      }

      keyVec[3] = "maxIter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "CG_maxIter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "resDirectly";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "CG_resDirectly", atoi(list[0].c_str()) );
      }
      break;

    case OLAS::GMRES:
      keyVec[2] = "gmres";

      keyVec[3] = "tol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "GMRES_epsilon", atof(list[0].c_str()) );
      }

      keyVec[3] = "maxIter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "GMRES_maxIter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "maxKrylovDim";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "GMRES_maxKrylovDim", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "GMRES_logging", (list[0] == "yes") );
      }
      break;

    case OLAS::MINRES:
      keyVec[2] = "gmres";

      keyVec[3] = "tol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "MINRES_epsilon", atof(list[0].c_str()) );
      }

      keyVec[3] = "maxIter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "MINRES_maxIter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "MINRES_logging", (list[0] == "yes") );
      }
      break;

    case OLAS::HYPRE_PCG:
      keyVec[2] = "hyprePCG";

      keyVec[3] = "tol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREPCG_epsilon", atof(list[0].c_str()) );
      }

      keyVec[3] = "maxIter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREPCG_maxIter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREPCG_logging", atoi(list[0].c_str()) );
      }

      keyVec[3] = "printLevel";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREPCG_printLevel", atoi(list[0].c_str()) );
      }

      keyVec[3] = "twonorm";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREPCG_twoNorm", (list[0] == "yes") );
      }
      break;

    case OLAS::HYPRE_GMRES:
      keyVec[2] = "hypreGMRES";

      keyVec[3] = "tol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREGMRES_epsilon", atof(list[0].c_str()) );
      }

      keyVec[3] = "maxIter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREGMRES_maxIter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "maxKrylovDim";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREGMRES_maxKrylovDim", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREGMRES_logging", atoi(list[0].c_str()) );
      }

      keyVec[3] = "printLevel";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREGMRES_printLevel", atoi(list[0].c_str()) );
      }
      break;

    case OLAS::HYPRE_BICGSTAB:
      keyVec[2] = "hypreBICGSTAB";

      keyVec[3] = "tol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREBICGSTAB_epsilon", atof(list[0].c_str()) );
      }

      keyVec[3] = "maxIter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREBICGSTAB_maxIter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREBICGSTAB_logging", atoi(list[0].c_str()) );
      }

      keyVec[3] = "printLevel";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "HYPREBICGSTAB_printLevel", atoi(list[0].c_str()) );
      }
      break;

    case OLAS::LAPACK_LU:
      keyVec[2] = "lapackLU";

      keyVec[3] = "tryScaling";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LAPACKLU_tryScaling", (list[0] == "yes") );
      }

      keyVec[3] = "refineSol";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LAPACKLU_refineSol", (list[0] == "yes") );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LAPACKLU_logging", (list[0] == "yes") );
      }
      break;

    case OLAS::LAPACK_LL:
      keyVec[2] = "lapackLL";

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LAPACKLL_logging", (list[0] == "yes") );
      }
      break;

    case OLAS::LU_SOLVER:
      keyVec[2] = "directLU";

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LUSOLVER_logging", (list[0] == "yes") );
      }

      keyVec[3] = "saveFacFile";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "CROUT_saveFacToFile", true );
        olas->SetValue( "CROUT_facFileName", list[0] );
      }
      break;

    case OLAS::LDL_SOLVER:
      keyVec[2] = "directLDL";

      keyVec[3] = "itRefSteps";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LDLSOLVER_itRefSteps", atoi(list[0].c_str()) );
      }

      keyVec[3] = "itRefVerbosity";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LDLSOLVER_itRefVerbosity", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LDLSOLVER_logging", (list[0] == "yes") );
      }

      keyVec[3] = "saveFacFile";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "LDLSOLVER_saveFacToFile", true );
        olas->SetValue( "LDLSOLVER_facFileName", list[0] );
        keyVec[3] = "savePatternOnly";
        cfs->GetList( keyVec, attrVec, valVec, list );
        if( list.GetSize() == 1 ) {
          olas->SetValue( "LDLSOLVER_facPatternOnly", (list[0] == "yes") );
        }
      }
      break;

    case OLAS::PARDISO:
      keyVec[2] = "pardiso";

      keyVec[3] = "posDef";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARDISO_posDef", (list[0] == "yes") );
      }

      keyVec[3] = "hermitean";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARDISO_hermitean", (list[0] == "yes") );
      }

      keyVec[3] = "symStruct";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARDISO_symStructure", (list[0] == "yes") );
      }

      keyVec[3] = "ordering";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        OLAS::ReorderingType ordering;
        OLAS::String2Enum( list[0], ordering );
        olas->SetValue( "PARDISO_ordering", ordering );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARDISO_logging", (list[0] == "yes") );
      }

      keyVec[3] = "stats";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARDISO_stats", (list[0] == "yes") );
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
    keyVec  = "linearSystems", "system", "solver", "stoppingRule", "type";
    attrVec = "", "name", "", "";
    valVec  = "", pdename, "", "";
    std::string stopCrit;
    cfs->Get( keyVec, attrVec, valVec, stopCrit );
    OLAS::StopCritType stopRule;
    OLAS::String2Enum( stopCrit, stopRule );
    olas->SetValue( "StoppingCriterion", stopRule );
  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::SetPrecondParams( std::string pdename,
                                        BaseParamHandler *cfs,
                                        OLAS::OLAS_Params *olas,
                                        OLAS::PrecondType pType ) {

    ENTER_FCN( "CFSOLASParams::SetPrecondParams" );

    // We need three vectors for our parameter queries
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // Basic part of query will always be the same
    keyVec  = "linearSystems", "system", "", "";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    StdVector<std::string> list;

    switch( pType ) {

      // The easiest cases. No parameters exist.
    case OLAS::ID:
    case OLAS::NOPRECOND:
    case OLAS::JACOBI:
    case OLAS::ILU0:
    case OLAS::IC0:
    case OLAS::HYPRE_JACOBI:
      break;

    case OLAS::HYPRE_ILU:
      keyVec[2] = "hypreILU";

      keyVec[3] = "level";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "EUCLID_level",  atoi(list[0].c_str()) );
      }

      keyVec[3] = "stats";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "EUCLID_stats", (list[0] == "yes") );
      }

      keyVec[3] = "memory";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        Integer memory = (list[0] == "yes");
        olas->SetValue( "EUCLID_memory", memory );
      }
      break;

    case OLAS::HYPRE_SPAI:
      keyVec[2] = "hypreSPAI";

      keyVec[3] = "thresh";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARASAILS_thresh", atof(list[0].c_str()) );
      }

      keyVec[3] = "levels";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARASAILS_levels", atoi(list[0].c_str()) );
      }

      keyVec[3] = "filter";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARASAILS_filter", atof(list[0].c_str()) );
      }

      keyVec[3] = "symmetry";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "PARASAILS_symmetry", atoi(list[0].c_str()) );
      }

      keyVec[3] = "loadBalance";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        Integer balance = (list[0] == "yes");
        olas->SetValue( "PARASAILS_loadBalance", balance );
      }

      keyVec[3] = "reuse";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        Integer balance = (list[0] == "yes");
        olas->SetValue( "PARASAILS_reuse", balance );
      }
      break;

    case OLAS::HYPRE_AMG:
      keyVec[2] = "hypreAMG";

      keyVec[3] = "maxLevels";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "BOOMERAMG_maxLevels", atoi(list[0].c_str()) );
      }

      keyVec[3] = "alpha";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "BOOMERAMG_alpha", atof(list[0].c_str()) );
      }

      keyVec[3] = "numSweeps";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "BOOMERAMG_numSweeps", atoi(list[0].c_str()) );
      }

      keyVec[3] = "cycleType";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        Integer cycleType = list[0] == "W-cycle" ? 2 : 1;
        olas->SetValue( "BOOMERAMG_cycleType", cycleType );
      }
      break;

    case OLAS::MG:
      keyVec[2] = "MG";

      keyVec[3] = "maxCoarseDepend";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_MaxCoarseDependency", atoi(list[0].c_str()) );
      }

      keyVec[3] = "minSystemSize";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_MinSystemSize", atoi(list[0].c_str()) );
      }

      keyVec[3] = "numPreSmooth";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_NumPreSmoothing", atoi(list[0].c_str()) );
      }

      keyVec[3] = "numPostSmooth";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_NumPostSmoothing", atoi(list[0].c_str()) );
      }

      keyVec[3] = "cycleParam";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_CycleParameter", atoi(list[0].c_str()) );
      }

      keyVec[3] = "alpha";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_CoarseningAlpha", atof(list[0].c_str()) );
      }

      keyVec[3] = "strongDiagRatio";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_StrongDiagRatio", atof(list[0].c_str()) );
      }

      keyVec[3] = "forceFineRatio";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_ForceFineRatio", atof(list[0].c_str()) );
      }

      keyVec[3] = "directSolver";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        if ( list[0] == "lapackLU" ) {
          olas->SetValue( "AMG_DirectSolver", OLAS::LAPACK_LU );
        }
        else if ( list[0] == "directLDL" ) {
          olas->SetValue( "AMG_DirectSolver", OLAS::LDL_SOLVER );
        }
        else if ( list[0] == "pardiso" ) {
          olas->SetValue( "AMG_DirectSolver", OLAS::PARDISO );
        }
        else {
          olas->SetValue( "AMG_DirectSolver", OLAS::NOSOLVER );
        }
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "AMG_logging", (list[0] == "yes") );
      }
      break;

    case OLAS::ILUK:
      keyVec[2] = "ILUK";

      keyVec[3] = "level";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILUK_level", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILUK_logging", (list[0] == "yes") );
      }

      keyVec[3] = "saveFacFile";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "CROUT_saveFacToFile", true );
        olas->SetValue( "CROUT_facFileName", list[0] );
      }
      break;

    case OLAS::ILDLK:
      keyVec[2] = "ILDLK";

      keyVec[3] = "level";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_level", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLKPRECOND_logging", (list[0] == "yes") );
      }

      keyVec[3] = "saveFacFile";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
        olas->SetValue( "ILDLPRECOND_facFileName", list[0] );
        keyVec[3] = "savePatternOnly";
        cfs->GetList( keyVec, attrVec, valVec, list );
        if( list.GetSize() == 1 ) {
          olas->SetValue( "ILDLPRECOND_facPatternOnly", (list[0] == "yes") );
        }
      }
      break;

    case OLAS::ILDLTP:
      keyVec[2] = "ILDLTP";

      keyVec[3] = "threshold";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_tau", atof(list[0].c_str()) );
      }

      keyVec[3] = "fillVal";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_fillVal", atoi(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLTPPRECOND_logging", (list[0] == "yes") );
      }

      keyVec[3] = "saveFacFile";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
        olas->SetValue( "ILDLPRECOND_facFileName", list[0] );
        keyVec[3] = "savePatternOnly";
        cfs->GetList( keyVec, attrVec, valVec, list );
        if( list.GetSize() == 1 ) {
          olas->SetValue( "ILDLPRECOND_facPatternOnly", (list[0] == "yes") );
        }
      }
      break;

    case OLAS::ILDLCN:
      keyVec[2] = "ILDLCN";

      keyVec[3] = "threshold";
      cfs->GetList( keyVec, attrVec, valVec, list );
      cfs->GetList( "threshold", list, pdename, "ILDLCN" );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_tau", atof(list[0].c_str()) );
      }

      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLCNPRECOND_logging", (list[0] == "yes") );
      }

      keyVec[3] = "saveFacFile";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
        olas->SetValue( "ILDLPRECOND_facFileName", list[0] );
        keyVec[3] = "savePatternOnly";
        cfs->GetList( keyVec, attrVec, valVec, list );
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
  //   SetEigenSolverParams
  // *******************
  void CFSOLASParams::SetEigenSolverParams( std::string pdename, 
                                            BaseParamHandler *cfs,
                                            OLAS::OLAS_Params *olas, 
                                            OLAS::EigenSolverType sType ) {
    
    ENTER_FCN( "CFSOLASParams::SetEigenSolverParams" );

    // We need three vectors for our parameter queries
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    
    // Basic part of query will always be the same
    keyVec  = "linearSystems", "system", "", "";
    attrVec = "", "name", "";
    valVec  = "", pdename, "";

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    StdVector<std::string> list;

    switch (sType) {
    case OLAS::ARPACK:
      keyVec[2] = "arpack";
      
      // tolerance
      keyVec[3] = "tolerance";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ARPACK_tolerance", atof(list[0].c_str()) );
      }

      // maximum number of iterations
      keyVec[3] = "maxIt";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ARPACK_maxIt", atoi(list[0].c_str()) );
      }

      // number of arnoldi vectors
      keyVec[3] = "numVec";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ARPACK_numVec", atoi(list[0].c_str()) );
      }

      // which eigenvalues / vectorss
      keyVec[3] = "which";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ARPACK_which", list[0] );
      }

      // logging
      keyVec[3] = "logging";
      cfs->GetList( keyVec, attrVec, valVec, list );
      if( list.GetSize() == 1 ) {
        olas->SetValue( "ARPACK_logging", (list[0] == "yes") );
      }

      break;

    case OLAS::SUBSPACE :
      // Not yet implemented
      break;

    case OLAS::NOEIGENSOLVER :
      // Nothing to do here
      break;
      
      // If this point is reached, it indicates that something is broken.
      // Probably not all solvers that SetParams allows are yet implemented
      // here?

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
                              OLAS::SolverType &sType,
                              OLAS::PrecondType &pType,
                              OLAS::MatrixStorageType &mType,
                              OLAS::MatrixEntryType &eType,
                              OLAS::ReorderingType &rType,
                              AnalysisType analysisType,
                              bool allowChangeOfReordering ) {

    ENTER_FCN( "CFSOLASParams::Expert" );

    std::string warn;

    // ==============
    //  Solver stuff
    // ==============

    // If no solver was specified use a direct one
    // Currently always use LU until LDL is available
    if ( sType == OLAS::NOSOLVER ) {
      // sType = LU_SOLVER;
      sType = OLAS::LDL_SOLVER;
    }


    // =======================
    //  Preconditioner stuff
    // =======================

    // For direct solver we need no preconditioner
    if ( sType == OLAS::DIRECT || sType == OLAS::LAPACK_LU 
         && pType != OLAS::NOPRECOND ) {
      warn = "Expert: Re-setting preconditioner type to 'NOPRECOND'";
      Info->Warning( warn );
      pType = OLAS::NOPRECOND;
    }
 
    // Hypre solvers only work together with hypre preconditioners
    if ( (sType == OLAS::HYPRE_PCG || sType == OLAS::HYPRE_GMRES ||
          sType == OLAS::HYPRE_BICGSTAB)

         && !( pType == OLAS::NOPRECOND    ||
               pType == OLAS::ID           ||
               pType == OLAS::HYPRE_JACOBI ||
               pType == OLAS::HYPRE_SPAI   ||
               pType == OLAS::HYPRE_ILU    ||
               pType == OLAS::HYPRE_AMG ) ) {

      warn = "Expert: Re-setting preconditioner type to 'ID'";
      Info->Warning( warn );
      pType = OLAS::ID;
    }

    // ===============
    //  Matrix stuff
    // ===============

    // Lapack solvers want their own matrix format
    if ( sType == OLAS::LAPACK_LU ) {
      if ( mType !=OLAS:: LAPACK_GBMATRIX ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          warn = "Expert: Re-setting matrix storage type to ";
          warn += "'LAPACK_GBMATRIX'";
          Info->Warning( warn );
        }
        else {
          Info->PrintF( pdename,
                        "Expert: Using LAPACK_GBMATRIX as storage type\n" );
        }
        mType = OLAS::LAPACK_GBMATRIX;
      }
    }

    // The direct solver LU_SOLVER expects a CRS matrix
    else if ( sType == OLAS::LU_SOLVER ) {
      if ( mType != OLAS::SPARSE_NONSYM ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          (*warning) << "Expert: Changing matrix storage type from "
                     << Enum2String( mType ) << " to SPARSE_NONSYM for "
                     << Enum2String( sType ) << " solver";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for direct solver\n" );
        }
        mType = OLAS::SPARSE_NONSYM;
      }
    }

    // The direct solver LDL_SOLVER expects an SCRS matrix
    else if ( sType == OLAS::LDL_SOLVER ) {
      if ( mType != OLAS::SPARSE_SYM ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          (*warning) << "Expert: Changing matrix storage type from "
                     << Enum2String( mType ) << " to SPARSE_SYM for "
                     << Enum2String( sType ) << " solver";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_SYM as storage "
                        "type for directLDL solver\n" );
        }
        mType = OLAS::SPARSE_SYM;
      }
    }

    // Hypre solvers want their own matrix format
    else if ( sType == OLAS::HYPRE_PCG || sType == OLAS::HYPRE_GMRES ||
              sType == OLAS::HYPRE_BICGSTAB ) {
      if ( mType != OLAS::HYPRE_MATRIX ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          warn = "Expert: Re-setting matrix storage type to 'HYPRE_MATRIX'\n";
          Info->Warning( warn );
        }
        else {
          Info->PrintF( pdename, "Expert: Using HYPRE_MATRIX as storage "
                        "type\n" );
        }
        mType = OLAS::HYPRE_MATRIX;
      }
      if ( eType != OLAS::DOUBLE ) {
        warn = "Expert: Re-setting matrix entry type to DOUBLE\n";
        Info->Warning( warn );
        eType = OLAS::DOUBLE;
      }
    }

    // ILU-type preconditioners expect CRS matrices
    if ( pType == OLAS::ILU0 || pType == OLAS::ILUK ) {
      if ( mType != OLAS::SPARSE_NONSYM ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          (*warning) << "Expert: Changing matrix storage type from "
                     << Enum2String( mType ) << " to SPARSE_NONSYM for "
                     << Enum2String( pType ) << " preconditioner";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for preconditioner\n" );
        }
        mType = OLAS::SPARSE_NONSYM;
      }
    }

    // ILDL-type preconditioners expect SCRS matrices
    if ( pType == OLAS::ILDLK ) {
      if ( mType != OLAS::SPARSE_SYM ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          (*warning) << "Expert: Changing matrix storage type from "
                     << Enum2String( mType ) << " to SPARSE_SYM for "
                     << Enum2String( pType ) << " preconditioner";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_SYM as storage "
                        "type for preconditioner\n" );
        }
        mType = OLAS::SPARSE_SYM;
      }
    }

    // The MG preconditioner requires a CRS matrix
    if ( pType == OLAS::MG ) {
      if ( mType != OLAS::SPARSE_NONSYM ) {
        if ( mType != OLAS::NOSTORAGETYPE ) {
          (*warning) << "Expert: Changing matrix storage type from "
                     << Enum2String( mType ) << " to SPARSE_NONSYM for "
                     << Enum2String( pType ) << " preconditioner";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for MG preconditioner\n" );
        }
        mType = OLAS::SPARSE_NONSYM;
      }
    }

    // If no storage type was yet assigned use plain SCRS
    if ( mType == OLAS::NOSTORAGETYPE ) {
      mType = OLAS::SPARSE_SYM;
      Info->PrintF( pdename, "Expert: Using SPARSE_SYM as matrix "
                    "storage type\n" );
    }

    // In case of a harmonic analysis or parameter identification
    // we assume complex matrix entries by default
    if ( analysisType == HARMONIC
         && eType != OLAS::COMPLEX ) {
      eType = OLAS::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for HARMONIC analysis\n" );
    }

    std::string analysis;
    cfs->Get( "type", analysis, "analysis" );
    if ( analysis == "paramIdent" && eType != OLAS::COMPLEX ) {
      eType = OLAS::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for parameter identification\n" );
    }

    if ( analysisType == MULTIHARMONIC && eType != OLAS::COMPLEX ) {
      eType = OLAS::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for parameter identification\n" );
    }

    if ( analysisType == MULTIHARMONIC ) {
      if (( analysisType == HARMONIC
            ||analysisType == MULTIHARMONIC) &&
          eType != OLAS::COMPLEX ) {
        eType = OLAS::COMPLEX;
        Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type, "
                      "harmonic part of multi-sequence analysis\n" );
      }
    }


    // NOTE: We only try to determine a re-ordering, if the caller did
    //       explicitely allow us to do so
    if ( allowChangeOfReordering == true ) {


      // ========================
      //  Reordering for Solvers
      // ========================

      // For the sparse direct solvers implemented in OLAS, use METIS
      // re-ordering if available and SLOAN otherwise

#ifdef USE_METIS
      if ( sType == OLAS::LU_SOLVER || sType == OLAS::LDL_SOLVER ) {
        if ( rType == OLAS::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'METIS'\n" );
          rType = OLAS::METIS;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to METIS\n" );
          rType = OLAS::METIS;
        }
      }
#else
      if ( sType == OLAS::LU_SOLVER || sType == OLAS::LDL_SOLVER ) {
        if ( rType == OLAS::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = OLAS::SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = OLAS::SLOAN;
        }
      }
#endif

      // For the LAPACK solvers use SLOAN re-ordering
      if ( sType == OLAS::LAPACK_LU || sType ==OLAS:: LAPACK_LL ) {
        if ( rType == OLAS::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = OLAS::SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = OLAS::SLOAN;
        }
      }

      // For Pardiso do not re-order (since Pardiso does this better itself)
      if ( sType == OLAS::PARDISO && rType != OLAS::NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = OLAS::NOREORDERING;
      }


      // ================================
      //  Reordering for Preconditioners
      // ================================

      // For all advanced ILU type preconditioners use SLOAN re-ordering
      if ( pType == OLAS::ILUK || pType == OLAS::ILDLK 
           || pType == OLAS::HYPRE_ILU ) {
        if ( rType == OLAS::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = OLAS::SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = OLAS::SLOAN;
        }
      }

      // For ILU0 and JACOBI we do not use re-ordering.
      else if ( ( pType == OLAS::ILU0 || pType == OLAS::JACOBI 
                  || pType == OLAS::HYPRE_JACOBI )
                && rType != OLAS::NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = OLAS::NOREORDERING;
      }
    }
  }

}
