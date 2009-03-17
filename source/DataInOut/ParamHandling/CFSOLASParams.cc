// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string.h>

#include <def_use_metis.hh>

#include "OLAS/algsys/olasparams.hh"
#include "General/environment.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "Domain/domain.hh"
#include "Driver/singleDriver.hh"
#include "Driver/assemble.hh"



namespace CoupledField {

  // *************
  //   SetParams
  // *************
  void CFSOLASParams::SetParams( std::string pdename,
                                 ParamNode *cfs,
                                 OLAS_Params *olas,
                                 BasePDE::AnalysisType analysisType,
                                 Assemble * assemble,
                                 bool overrideExpert ) {

    // Fetch setup node
    ParamNode * setupNode = NULL;
    if( cfs )
      setupNode = cfs->Get("setup", false );


    // First determine approach for handling inhomogeneous Dirichlet
    // boundary conditions
    bool usingPenalty = true;
    {
      std::string aux = "penalty";
      if( setupNode )
        setupNode->Get("idbcHandling", aux, false );

      usingPenalty = aux == "penalty" ? true : false;
    }

    // Determine the type of solver for this PDE
    std::string sTypeString = "expertsChoice";
    SolverType sType;
    ParamNode * solverNode = NULL;
    if( cfs )
      solverNode = cfs->Get("solver", false);

    if( solverNode)
      solverNode->Get("type", sTypeString, false );
    if ( sTypeString == "expertsChoice" ) {
      if ( overrideExpert ) {
        EXCEPTION( "You cannot specify expertsChoice as solver type "
                   << "and at the same time specify overrideExpert! "
                   << "This would leave the solver type undefined!" );
      }
      sTypeString = "no solver";
    }
    String2Enum( sTypeString, sType );

    // Now determine the type of preconditioner for this PDE
    std::string pTypeString = "noPrecond";
    PrecondType pType;

    if( solverNode )
      solverNode->Get("precond", pTypeString, false );

    String2Enum( pTypeString, pType );

    // Now determine, if we are running an eigenfrequency analysis and if
    // yes, get the type of eigenvalue solver
    std::string esTypeString = "expertsChoice";
    EigenSolverType esType = NOEIGENSOLVER;

    if ( analysisType == BasePDE::EIGENFREQUENCY ) {
      if( cfs) {
        ParamNode * eigenSolverNode = cfs->Get("eigenSolver", false );
        if( eigenSolverNode )
          eigenSolverNode->Get("type", esTypeString, false );
      }
      if ( esTypeString == "expertsChoice" ) {
        if ( overrideExpert ) {
          EXCEPTION( "You cannot specify expertsChoice as EigenSolver type "
                     << "and at the same time specify overrideExpert! "
                     << "This would leave the solver type undefined!" );
        }
        esTypeString = "no eigensolver";
      }

      String2Enum( esTypeString, esType );
    }

    // Before we try to determine matrix properties, we first should
    // find out, what type of linear system we are using. We do this
    // by trying to determine the symmetry attribute of the sbmMatrix
    // element
    bool sbmSymmetry = false;
    bool stdSystem = true;
    std::string sbmSymmetric;
    ParamNode * matrixNode = NULL;

    if( cfs ) {
      matrixNode = cfs->Get("matrix", false);
      if( !matrixNode ) {
        matrixNode = cfs->Get("sbmMatrix", false);
        if( matrixNode ) {
          stdSystem = false;
          matrixNode->Get("symmetric", sbmSymmetric);
          sbmSymmetry = sbmSymmetric  == "yes";
        }
      }
    }

    // Determine matrix entry type
    BaseMatrix::EntryType eType = BaseMatrix::DOUBLE;
    if( matrixNode ) {
      std::string eMatString =  matrixNode->Get("entry")->AsString();
      if ( eMatString == "double" ) {
        eType = BaseMatrix::DOUBLE;
      }
      else if ( eMatString == "complex" ) {
        eType = BaseMatrix::COMPLEX;
      }
    }
    
    // Determine calculation of condition number
    bool calcCondition = false;
    if( matrixNode ) { 
      matrixNode->Get("calcConditionNumber", calcCondition, false );
    }

    // Following stuff only for required for standard systems
    BaseMatrix::StorageType mType  = BaseMatrix::NOSTORAGETYPE;
    ReorderingType orderType = NOREORDERING;
    bool allowChangeOfReordering = false;


    if ( stdSystem == true ) {
      std::string mMatString = "expertsChoice";
      if( matrixNode ) {
        // Next determine the matrix type for this PDE
        matrixNode->Get( "storage", mMatString );
      }
      if ( mMatString == "expertsChoice" ) {
        if ( overrideExpert ) {
          EXCEPTION( "You cannot specify expertsChoice as storage type "
                     << "and at the same time specify overrideExpert! "
                     << "This would leave the storage type undefined!" );
        }
        mMatString = "noStorageType";
      }
      mType = BaseMatrix::storageType.Parse( mMatString );

      // Type of re-odering
      std::string orderString = "expertsChoice";
      if( matrixNode ) {
        matrixNode->Get( "reordering", orderString );
      }
      if ( orderString == "expertsChoice" ) {
        if ( overrideExpert ) {
          EXCEPTION( "You cannot specify expertsChoice as re-ordering "
                     << "strategy and at the same time specify overrideExpert! "
                     << "This would leave the re-ordering strategy undefined!" );
        }
        orderString = "noReordering";
        allowChangeOfReordering = true;
      }
      String2Enum( orderString, orderType );
    }

    // Let expert module modify the settings
    if ( overrideExpert == false && stdSystem == true ) {
      CFSOLASParams::Expert( cfs, pdename, esType, sType, pType, mType, eType,
                             orderType, analysisType, assemble,
                             allowChangeOfReordering );
    }

    // Insert information into OLAS_Params object
    olas->SetValue( "Solver"                 , sType              );
    olas->SetValue( "Precond"                , pType              );
    olas->SetValue( "EigenSolver"            , esType             );
    olas->SetValue( "MatrixStructureType"    , BaseMatrix::SPARSE_MATRIX    );
    olas->SetValue( "MatrixStorageType"      , mType              );
    olas->SetValue( "MatrixEntryType"        , eType              );
    olas->SetValue( "GRAPH_reordering"       , orderType          );
    olas->SetValue( "UsingPenaltyFormulation", usingPenalty       );
    olas->SetValue( "SystemName"             , pdename            );
    olas->SetValue( "SBM_Symmetry"           , sbmSymmetry        );
    olas->SetValue( "CalcConditionNumber"    , calcCondition      );

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
  }


  // *******************
  //   SetSolverParams
  // *******************
  void CFSOLASParams::SetSolverParams( std::string pdename,
                                       ParamNode *cfs,
                                       OLAS_Params *olas,
                                       SolverType sType ) {


    // check if base element for solver is present
    if( !cfs) return;

    ParamNode * bsNode = cfs->Get("solver", false );
    if ( !bsNode ) return;

    // parameter node for solver
    ParamNode * sNode = NULL;
    std::string val;

    switch( sType ) {

    case DIAGSOLVER:
      break;

    case CG:
      break;

    case GMRES:
      sNode = bsNode->Get("gmres", false);
      if(!sNode) return;

      if( sNode->Has("tol") )
        olas->SetValue( "GMRES_epsilon",
                        sNode->Get("tol")->AsDouble() );

      if( sNode->Has("maxIter") )
        olas->SetValue( "GMRES_maxIter",
                        sNode->Get( "maxIter")->AsInt() );

      if( sNode->Has("maxKrylovDim") )
        olas->SetValue( "GMRES_maxKrylovDim",
                        sNode->Get("maxKrylovDim")->AsInt() );

      if( sNode->Has("logging") )
        olas->SetValue( "GMRES_logging",
                        sNode->Get("logging")->AsBool() );
      break;

    case MINRES:
      sNode = bsNode->Get("minres", false);
      if(!sNode) return;

      if( sNode->Has("tol") )
        olas->SetValue( "MINRES_epsilon",
                        sNode->Get("tol")->AsDouble() );

      if( sNode->Has("maxIter") )
        olas->SetValue( "MINRES_maxIter",
                        sNode->Get("maxIter")->AsInt() );

      if( sNode->Has("logging") )
        olas->SetValue( "MINRES_logging",
                        sNode->Get("logging")->AsBool() );
      break;

    case LAPACK_LU:
      sNode = bsNode->Get("lapackLU", false);
      if(!sNode) return;

      if( sNode->Has("tryScaling") )
        olas->SetValue( "LAPACKLU_tryScaling",
                        sNode->Get("tryScaling")->AsBool() );

      if( sNode->Has("refineSol") )
        olas->SetValue( "LAPACKLU_refineSol",
                        sNode->Get("refineSol")->AsBool() );

      if( sNode->Has( "logging") )
        olas->SetValue( "LAPACKLU_logging",
                        sNode->Get( "logging")->AsBool() );
      break;

    case LAPACK_LL:
      sNode = bsNode->Get("lapackLL", false);
      if(!sNode) return;

      if( sNode->Has( "logging") )
        olas->SetValue( "LAPACKLL_logging",
                        sNode->Get("logging")->AsBool() );
      break;

    case LU_SOLVER:
      sNode = bsNode->Get("directLU", false);
      if(!sNode) return;

      if( sNode->Has( "logging") )
        olas->SetValue( "LUSOLVER_logging",
                        sNode->Get( "logging")->AsBool() );

      if( sNode->Has("saveFacFile") ) {
        olas->SetValue("CROUT_saveFacToFile", true );
        olas->SetValue("CROUT_facFileName",
                       sNode->Get("saveFacFile")->AsString() );
      }
      break;

    case LDL_SOLVER:
      sNode = bsNode->Get("directLDL", false);
      if(!sNode) return;

      if( sNode->Has( "itRefSteps") )
        olas->SetValue("LDLSOLVER_itRefSteps",
                       sNode->Get( "itRefSteps")->AsInt() );

      if( sNode->Has("itRefVerbosity") )
        olas->SetValue( "LDLSOLVER_itRefVerbosity",
                       sNode->Get("itRefVerbosity")->AsInt() );

      if( sNode->Has("logging") )
        olas->SetValue("LDLSOLVER_logging",
                       sNode->Get("logging")->AsBool() );

      if( sNode->Has("saveFacFile") ) {
        olas->SetValue("LDLSOLVER_saveFacToFile", true );
        olas->SetValue("LDLSOLVER_facFileName",
                       sNode->Get("saveFacFile")->AsString() );
      }

      if( sNode->Has("savePatternOnly" ) ) {
        olas->SetValue( "LDLSOLVER_facPatternOnly",
                        sNode->Get("savePatternOnly")->AsBool());
      }
      break;

    case PARDISO:
      sNode = bsNode->Get("pardiso", false);
      if(!sNode) return;

      if( sNode->Has("posDef") )
        olas->SetValue("PARDISO_posDef",
                       sNode->Get("posDef")->AsBool() );

      if( sNode->Has("hermitean") )
        olas->SetValue("PARDISO_hermitean",
                       sNode->Get("hermitean")->AsBool() );

      if( sNode->Has("symStruct") )
        olas->SetValue("PARDISO_symStructure",
                       sNode->Get("symStruct")->AsBool() );

      if( sNode->Has("ordering") ) {
        ReorderingType ordering;
        String2Enum( sNode->Get("ordering")->AsString(), ordering );
        olas->SetValue( "PARDISO_ordering", ordering );
      }

      if( sNode->Has("logging") )
        olas->SetValue("PARDISO_logging",
                       sNode->Get("logging")->AsBool() );

      if( sNode->Has("stats") )
        olas->SetValue("PARDISO_stats",
                       sNode->Get("stats")->AsBool() );
      break;

    case ILUPACK_SOLVER:
      break;

       // If this point is reached, it indicates that something is broken.
       // Probably not all solvers that SetParams allows are yet implemented
       // here?
     default:
       EXCEPTION( "Internal error. Maybe a missing implementation of a case?" );
       break;

    }


    // Now set the stopping rule
    std::string stopCrit = "relNormRes0";
    ParamNode * stopRuleNode = bsNode->Get("stoppingRule", false );
    if( stopRuleNode ) {
      stopCrit = stopRuleNode->Get("type")->AsString();
    }
    StopCritType stopRule;
    String2Enum( stopCrit, stopRule );
    olas->SetValue( "StoppingCriterion", stopRule );
  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::SetPrecondParams( std::string pdename,
                                        ParamNode *cfs,
                                        OLAS_Params *olas,
                                        PrecondType pType ) {



    // check for base element of solver
    if( !cfs) return;
    ParamNode * bsNode = cfs->Get("solver", false );
    if( !bsNode ) return;

    // parameter node for preconditioner
    ParamNode * pNode = NULL;

    switch( pType ) {

      // The easiest cases. No parameters exist.
    case ID:
    case NOPRECOND:
    case JACOBI:
    case ILU0:
    case IC0:
      break;

    case MG:

      pNode = bsNode->Get("MG", false );
      if( !pNode) return;

      if( pNode->Has( "maxCoarseDepend" ) )
        olas->SetValue( "AMG_MaxCoarseDependency",
                        pNode->Get( "maxCoarseDepend")->AsInt() );

      if( pNode->Has( "minSystemSize") )
        olas->SetValue("AMG_MinSystemSize",
                       pNode->Get( "minSystemSize")->AsInt() );

      if( pNode->Has("numPreSmooth" ) )
        olas->SetValue("AMG_NumPreSmoothing",
                       pNode->Get("numPreSmooth")->AsInt() );

      if( pNode->Has("numPostSmooth") )
        olas->SetValue("AMG_NumPostSmoothing",
                       pNode->Get("numPostSmooth")->AsInt() );

      if( pNode->Has("cycleParam") )
        olas->SetValue("AMG_CycleParameter",
                       pNode->Get("cycleParam")->AsInt() );

      if( pNode->Has("alpha") )
        olas->SetValue("AMG_CoarseningAlpha",
                       pNode->Get("alpha")->AsDouble() );

      if( pNode->Has("strongDiagRatio") )
        olas->SetValue("AMG_StrongDiagRatio",
                       pNode->Get("strongDiagRatio")->AsDouble() );

      if( pNode->Has("forceFineRatio") )
        olas->SetValue("AMG_ForceFineRatio",
                       pNode->Get("forceFineRatio")->AsDouble() );

      if( pNode->Has("directSolver") ) {
        std::string solverString;
        pNode->Get("directSolver", solverString );
        if ( solverString == "lapackLU" )
          olas->SetValue( "AMG_DirectSolver", LAPACK_LU );
        else if( solverString == "directLDL" )
          olas->SetValue( "AMG_DirectSolver", LDL_SOLVER );
        else if( solverString == "pardiso" )
          olas->SetValue( "AMG_DirectSolver", PARDISO );
        else
          olas->SetValue( "AMG_DirectSolver", NOSOLVER );
      }

      if( pNode->Has("logging") )
        olas->SetValue("AMG_logging",
                       pNode->Get("logging")->AsBool() );
      break;


    case ILUK:

      pNode = bsNode->Get("ILUK", false );
      if( !pNode ) return;


      if( pNode->Has("level") )
        olas->SetValue("ILUK_level",
                       pNode->Get("level")->AsInt() );

      if( pNode->Has("logging") )
        olas->SetValue("ILUK_logging",
                       pNode->Get("logging")->AsBool() );

      if( pNode->Has("saveFacFile") ) {
        olas->SetValue( "CROUT_saveFacToFile", true );
        olas->SetValue( "CROUT_facFileName",
                        pNode->Get("saveFacFile")->AsString() );
      }
      break;

    case ILDLK:

      pNode = bsNode->Get("ILDLK" );

      if( pNode->Has("level") )
        olas->SetValue("ILDLPRECOND_level",
                       pNode->Get("level")->AsInt() );

      if( pNode->Has("logging") )
        olas->SetValue("ILDLKPRECOND_logging",
                       pNode->Get("logging")->AsBool() );

      if( pNode->Has("saveFacFile") ) {
        olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
          olas->SetValue("ILDLPRECOND_facFileName",
                         pNode->Get("saveFacFile")->AsString() );
        if( pNode->Has("savePatternOnly") )
          olas->SetValue("ILDLPRECOND_facPatternOnly",
                         pNode->Get("savePatternOnly")->AsBool() );
      }
      break;

    case ILDLTP:

      pNode = bsNode->Get("ILDLTP");
      if( !pNode ) return;


      if( pNode->Has("threshold") )
        olas->SetValue("ILDLPRECOND_tau",
                       pNode->Get("threshold")->AsDouble() );

      if( pNode->Has("fillVal") )
        olas->SetValue("ILDLPRECOND_fillVal",
                       pNode->Get("fillVal")->AsInt() );

      if( pNode->Has("logging") )
        olas->SetValue("ILDLTPPRECOND_logging",
                       pNode->Get("logging")->AsBool() );

      if( pNode->Has("saveFacFile") ) {
        olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
          olas->SetValue("ILDLPRECOND_facFileName",
                         pNode->Get("saveFacFile")->AsString() );
        if( pNode->Has("savePatternOnly") )
          olas->SetValue("ILDLPRECOND_facPatternOnly",
                         pNode->Get("savePatternOnly")->AsBool() );
      }
      break;

    case ILDLCN:

      pNode = bsNode->Get("ILDLCN");
      if( !pNode) return;


      if( pNode->Has("threshold") )
        olas->SetValue("ILDLPRECOND_tau",
                       pNode->Get("threshold")->AsDouble() );

      if( pNode->Has("loging") )
        olas->SetValue("ILDLCNPRECOND_logging",
                       pNode->Get("loging")->AsBool() );

      if( pNode->Has("saveFacFile") ) {
        olas->SetValue( "ILDLPRECOND_saveFacToFile", true );
        olas->SetValue("ILDLPRECOND_facFileName",
                       pNode->Get("saveFacFile")->AsString() );
        if( pNode->Has("savePatternOnly") )
          olas->SetValue("ILDLPRECOND_facPatternOnly",
                         pNode->Get("savePatternOnly")->AsBool() );
      }
      break;

      case ILDL0:
        break;


      // If this point is reached, it indicates that something is broken.
      // Probably not all preconditioners that SetParams allows are yet
      // implemented here?
    default:
      EXCEPTION( "Internal error. Maybe a missing implementation of a case?" );
      break;
    }
  }

  // *******************
  //   SetEigenSolverParams
  // *******************
  void CFSOLASParams::SetEigenSolverParams( std::string pdename,
                                            ParamNode *cfs,
                                            OLAS_Params *olas,
                                            EigenSolverType sType ) {


    // check, if eigenfrequency solver is present at all
    if( !cfs) return;
    ParamNode * besNode = cfs->Get("eigenSolver", false );
    if( !besNode ) return;

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    ParamNode * esNode = NULL;

    switch (sType) {
    case ARPACK:
      esNode = besNode->Get("arpack", false );
      if( !esNode) return;

      // tolerance
      if( esNode->Has("tolerance") )
        olas->SetValue( "ARPACK_tolerance",
                        esNode->Get("tolerance")->AsDouble() );

      // maximum number of iterations
      if( esNode->Has("maxIt") )
        olas->SetValue( "ARPACK_maxIt",
                        esNode->Get("maxIt")->AsInt() );

      // number of arnoldi vectors
      if( esNode->Has("numVec") )
        olas->SetValue( "ARPACK_numVec",
                        esNode->Get("numVec")->AsInt() );

      // which eigenvalues / vectorss
      if( esNode->Has("which") )
        olas->SetValue( "ARPACK_which",
                        esNode->Get("which")->AsString() );

      // logging
      if( esNode->Has("logging") )
        olas->SetValue( "ARPACK_logging",
                        esNode->Get("logging")->AsBool() );
      break;

    case SUBSPACE :
      // Not yet implemented
      break;

    case NOEIGENSOLVER :
      // Nothing to do here
      break;

      // If this point is reached, it indicates that something is broken.
      // Probably not all solvers that SetParams allows are yet implemented
      // here?

    default:
      EXCEPTION( "Internal error. Maybe a missing implementation of a case?");
      break;
    }
  }


  // *******************
  //   Expert's Choice
  // *******************
  void CFSOLASParams::Expert( ParamNode *cfs,
                              std::string pdename,
                              EigenSolverType &esType,
                              SolverType &sType,
                              PrecondType &pType,
                              BaseMatrix::StorageType &mType,
                              BaseMatrix::EntryType &eType,
                              ReorderingType &rType,
                              BasePDE::AnalysisType analysisType,
                              Assemble * assemble,
                              bool allowChangeOfReordering ) {


    std::string warn;

    // ===========================
    //  Determine Matrix symmetry
    // ===========================
    bool symmetricMat = assemble->IsFEMatSymmetric();

    // ==============
    //  EigenSolver stuff
    // ==============
    // If no eigenvalue solver was specified, use arpack
    if ( esType == NOEIGENSOLVER ) {
      esType = ARPACK;
    }


    // ==============
    //  Solver stuff
    // ==============

    // If no solver was specified use a direct one
    // Here we have not to determine, if we are faced with
    // a symmetric or non-symmetric problem.

    if ( sType == NOSOLVER ) {
      if( symmetricMat ) {
        sType = LDL_SOLVER;
      } else {
        sType = LU_SOLVER;
      }
    } else {
      // In this case, the user has specified a symmetric
      // solver explicitly, although, the system is non-symmetric.
      // In this case we issue a warning
      if( sType == LDL_SOLVER
          && !symmetricMat ) {
        warn = "Expert: Although the system matrix is non-symmetric ";
        warn += "a symmetric solver was chosen!";
        Info->Warning(warn);
      }
    }


    // =======================
    //  Preconditioner stuff
    // =======================

    // For direct solver we need no preconditioner (ID)
    if ( (sType == LDL_SOLVER
          || sType == LU_SOLVER
          || sType == LAPACK_LU )
         && !(pType == ID
              || pType == NOPRECOND) ) {
      warn = "Expert: Re-setting preconditioner type to 'NOPRECOND'";
      Info->Warning( warn );
      pType = NOPRECOND;
    }



    // ===============
    //  Matrix stuff
    // ===============

    // Lapack solvers want their own matrix format
    if ( sType == LAPACK_LU ) {
      if ( mType != BaseMatrix::LAPACK_GBMATRIX ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          warn = "Expert: Re-setting matrix storage type to ";
          warn += "'LAPACK_GBMATRIX'";
          Info->Warning( warn );
        }
        else {
          Info->PrintF( pdename,
                        "Expert: Using LAPACK_GBMATRIX as storage type\n" );
        }
        mType = BaseMatrix::LAPACK_GBMATRIX;
      }
    }

    // The direct solver LU_SOLVER expects a CRS matrix
    else if ( sType == LU_SOLVER ) {
      if ( mType != BaseMatrix::SPARSE_NONSYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          Enum2String( sType, tmp );

          (*warning) << "Expert: Changing matrix storage type from "
                     << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_NONSYM for "
                     << tmp << " solver";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for direct solver\n" );
        }
        mType = BaseMatrix::SPARSE_NONSYM;
      }
    }

    // The direct solver LDL_SOLVER expects an SCRS matrix
    else if ( sType == LDL_SOLVER ) {
      if ( mType != BaseMatrix::SPARSE_SYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          Enum2String( sType, tmp );

          (*warning) << "Expert: Changing matrix storage type from "
                     << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_SYM for "
                     << tmp << " solver";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_SYM as storage "
                        "type for directLDL solver\n" );
        }
        mType = BaseMatrix::SPARSE_SYM;
      }
    }


    // ILU-type preconditioners expect CRS matrices
    if ( pType == ILU0 || pType == ILUK ) {
      if ( mType != BaseMatrix::SPARSE_NONSYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          Enum2String( pType, tmp );

          (*warning) << "Expert: Changing matrix storage type from "
                     << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_NONSYM for "
                     << tmp << " preconditioner";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for preconditioner\n" );
        }
        mType = BaseMatrix::SPARSE_NONSYM;
      }
    }

    // ILDL-type preconditioners expect SCRS matrices
    if ( pType == ILDLK ) {
      if ( mType != BaseMatrix::SPARSE_SYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          Enum2String( pType, tmp );

          (*warning) << "Expert: Changing matrix storage type from "
                     << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_SYM for "
                     << tmp << " preconditioner";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_SYM as storage "
                        "type for preconditioner\n" );
        }
        mType = BaseMatrix::SPARSE_SYM;
      }
    }

    // The MG preconditioner requires a CRS matrix
    if ( pType == MG ) {
      if ( mType != BaseMatrix::SPARSE_NONSYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          Enum2String( pType, tmp );

          (*warning) << "Expert: Changing matrix storage type from "
                     << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_NONSYM for "
                     << tmp << " preconditioner";
          Warning( __FILE__, __LINE__ );
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for MG preconditioner\n" );
        }
        mType = BaseMatrix::SPARSE_NONSYM;
      }
    }

    // If no storage type was yet assigned use plain SCRS
    if ( mType == BaseMatrix::NOSTORAGETYPE ) {
      if ( symmetricMat ) {
        mType = BaseMatrix::SPARSE_SYM;
        Info->PrintF( pdename, "Expert: Using SPARSE_SYM as matrix "
                      "storage type\n" );
      }
      else {
        mType = BaseMatrix::SPARSE_NONSYM;
        Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as matrix "
                      "storage type\n" );
      }
    }

    // Check, if any entries

    // In case of a harmonic analysis or parameter identification
    // we assume complex matrix entries by default
    if ( analysisType == BasePDE::HARMONIC
         && eType != BaseMatrix::COMPLEX ) {
      eType = BaseMatrix::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for HARMONIC analysis\n" );
    }

    // Special treatment for parameter identification process
    std::string analysis;
    UInt seqStep = domain->GetSingleDriver()->GetActSequenceStep();
    param->Get("sequenceStep", "index", GenStr(seqStep))->Get("analysis")
      ->GetChild()->GetName();

    if ( analysis == "paramIdent" && eType != BaseMatrix::COMPLEX ) {
      eType = BaseMatrix::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for parameter identification\n" );
    }

    // If an eigenvalue problem is solved, the entrytype of the matrices has
    // always to be DOUBLE
    if ( analysisType == BasePDE::EIGENFREQUENCY &&
         eType != BaseMatrix::DOUBLE ) {
      eType = BaseMatrix::COMPLEX;
      Info->PrintF( pdename, "Expert: Using DOUBLE as matrix entry type "
                    "for eigenFrequency simulation\n" );
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
      if ( sType == LU_SOLVER || sType == LDL_SOLVER ) {
        if ( rType == NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'METIS'\n" );
          rType = METIS;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to METIS\n" );
          rType = METIS;
        }
      }
#else
      if ( sType == LU_SOLVER || sType == LDL_SOLVER ) {
        if ( rType == NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = SLOAN;
        }
      }
#endif

      // For the LAPACK solvers use SLOAN re-ordering
      if ( sType == LAPACK_LU || sType == LAPACK_LL ) {
        if ( rType == NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = SLOAN;
        }
      }

      // For Pardiso do not re-order (since Pardiso does this better itself)
      if ( sType == PARDISO && rType != NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = NOREORDERING;
      }


      // ================================
      //  Reordering for Preconditioners
      // ================================

      // For all advanced ILU type preconditioners use SLOAN re-ordering
      if ( pType == ILUK || pType == ILDLK ) {
        if ( rType == NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = SLOAN;
        }
      }

      // For ILU0 and JACOBI we do not use re-ordering.
      else if ( ( pType == ILU0 || pType == JACOBI )
                && rType != NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = NOREORDERING;
      }
    }
  }

}
