// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string.h>

#include <def_use_metis.hh>

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
    OLAS::SolverType sType;
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
    OLAS::String2Enum( sTypeString, sType );

    // Now determine the type of preconditioner for this PDE
    std::string pTypeString = "noPrecond";
    OLAS::PrecondType pType;
    
    if( solverNode )
      solverNode->Get("precond", pTypeString, false );

    OLAS::String2Enum( pTypeString, pType );

    // Now determine, if we are running an eigenfrequency analysis and if
    // yes, get the type of eigenvalue solver
    std::string esTypeString = "expertsChoice";
    OLAS::EigenSolverType esType = OLAS::NOEIGENSOLVER;
    
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

      OLAS::String2Enum( esTypeString, esType );
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
    OLAS:: MatrixEntryType eType = OLAS::DOUBLE;
    if( matrixNode ) {
      std::string eMatString =  matrixNode->Get("entry")->AsString();
      if ( eMatString == "double" ) {
        eType = OLAS::DOUBLE;
      }
      else if ( eMatString == "complex" ) {
        eType = OLAS::COMPLEX;
      }
    }

    // Following stuff only for required for standard systems
    OLAS::MatrixStorageType mType  = OLAS::NOSTORAGETYPE;
    OLAS::ReorderingType orderType = OLAS::NOREORDERING;
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
      OLAS::String2Enum( mMatString, mType );
      
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
      OLAS::String2Enum( orderString, orderType );
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
  }


  // *******************
  //   SetSolverParams
  // *******************
  void CFSOLASParams::SetSolverParams( std::string pdename,
                                       ParamNode *cfs,
                                       OLAS::OLAS_Params *olas,
                                       OLAS::SolverType sType ) {


    // check if base element for solver is present
    if( !cfs) return;
    
    ParamNode * bsNode = cfs->Get("solver", false );
    if ( !bsNode ) return;
    
    // parameter node for solver
    ParamNode * sNode = NULL;
    std::string val;
    
    switch( sType ) {
      
    case OLAS::DIAGSOLVER:
      break;
      
    case OLAS::CG:
      break;
      
    case OLAS::GMRES:
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
      
    case OLAS::MINRES:
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

    case OLAS::LAPACK_LU:
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
      
    case OLAS::LAPACK_LL:
      sNode = bsNode->Get("lapackLL", false);
      if(!sNode) return;
      
      if( sNode->Has( "logging") ) 
        olas->SetValue( "LAPACKLL_logging",
                        sNode->Get("logging")->AsBool() );
      break;
      
    case OLAS::LU_SOLVER:
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
      
    case OLAS::LDL_SOLVER:
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
      
    case OLAS::PARDISO:
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
        OLAS::ReorderingType ordering;
        OLAS::String2Enum( sNode->Get("ordering")->AsString(), ordering );
        olas->SetValue( "PARDISO_ordering", ordering );
      }

      if( sNode->Has("logging") ) 
        olas->SetValue("PARDISO_logging",
                       sNode->Get("logging")->AsBool() );  

      if( sNode->Has("stats") ) 
        olas->SetValue("PARDISO_stats",
                       sNode->Get("stats")->AsBool() );  
      break;

    case OLAS::ILUPACK_SOLVER:
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
    OLAS::StopCritType stopRule;
    OLAS::String2Enum( stopCrit, stopRule );
    olas->SetValue( "StoppingCriterion", stopRule );
  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::SetPrecondParams( std::string pdename,
                                        ParamNode *cfs,
                                        OLAS::OLAS_Params *olas,
                                        OLAS::PrecondType pType ) {
    


    // check for base element of solver
    if( !cfs) return;
    ParamNode * bsNode = cfs->Get("solver", false );
    if( !bsNode ) return;
    
    // parameter node for preconditioner
    ParamNode * pNode = NULL;

    switch( pType ) {

      // The easiest cases. No parameters exist.
    case OLAS::ID:
    case OLAS::NOPRECOND:
    case OLAS::JACOBI:
    case OLAS::ILU0:
    case OLAS::IC0:
      break;

    case OLAS::MG:

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
          olas->SetValue( "AMG_DirectSolver", OLAS::LAPACK_LU );
        else if( solverString == "directLDL" )
          olas->SetValue( "AMG_DirectSolver", OLAS::LDL_SOLVER );
        else if( solverString == "pardiso" )
          olas->SetValue( "AMG_DirectSolver", OLAS::PARDISO );
        else 
          olas->SetValue( "AMG_DirectSolver", OLAS::NOSOLVER );
      }
      
      if( pNode->Has("logging") )
        olas->SetValue("AMG_logging",
                       pNode->Get("logging")->AsBool() );
      break;
      

    case OLAS::ILUK:

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
      
    case OLAS::ILDLK:

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

    case OLAS::ILDLTP:

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

    case OLAS::ILDLCN:

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
                                            OLAS::OLAS_Params *olas, 
                                            OLAS::EigenSolverType sType ) {
    

    // check, if eigenfrequency solver is present at all
    if( !cfs) return;
    ParamNode * besNode = cfs->Get("eigenSolver", false );
    if( !besNode ) return;

    // Determine which parameters have been set by the user
    // and insert them into the olasParams object.
    ParamNode * esNode = NULL;

    switch (sType) {
    case OLAS::ARPACK:
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
      EXCEPTION( "Internal error. Maybe a missing implementation of a case?");
      break;
    }
  }
  
  
  // *******************
  //   Expert's Choice
  // *******************
  void CFSOLASParams::Expert( ParamNode *cfs,
                              std::string pdename,
                              OLAS::EigenSolverType &esType,
                              OLAS::SolverType &sType,
                              OLAS::PrecondType &pType,
                              OLAS::MatrixStorageType &mType,
                              OLAS::MatrixEntryType &eType,
                              OLAS::ReorderingType &rType,
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
    if ( esType == OLAS::NOEIGENSOLVER ) {
      esType = OLAS::ARPACK;
    }


    // ==============
    //  Solver stuff
    // ==============
    
    // If no solver was specified use a direct one
    // Here we have not to determine, if we are faced with
    // a symmetric or non-symmetric problem.

    if ( sType == OLAS::NOSOLVER ) {
      if( symmetricMat ) {
        sType = OLAS::LDL_SOLVER;
      } else {
        sType = OLAS::LU_SOLVER;
      }
    } else {
      // In this case, the user has specified a symmetric
      // solver explicitly, although, the system is non-symmetric.
      // In this case we issue a warning
      if( sType == OLAS::LDL_SOLVER
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
    if ( (sType == OLAS::LDL_SOLVER 
          || sType == OLAS::LU_SOLVER
          || sType == OLAS::LAPACK_LU )
         && !(pType == OLAS::ID
              || pType == OLAS::NOPRECOND) ) {
      warn = "Expert: Re-setting preconditioner type to 'NOPRECOND'";
      Info->Warning( warn );
      pType = OLAS::NOPRECOND;
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

    // Check, if any entries
     
    // In case of a harmonic analysis or parameter identification
    // we assume complex matrix entries by default
    if ( analysisType == BasePDE::HARMONIC
         && eType != OLAS::COMPLEX ) {
      eType = OLAS::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for HARMONIC analysis\n" );
    }
    
    // Special treatment for parameter identification process
    std::string analysis;
    UInt seqStep = domain->GetSingleDriver()->GetActSequenceStep();
    param->Get("sequenceStep", "index", GenStr(seqStep))->Get("analysis")
      ->GetChild()->GetName();
    
    if ( analysis == "paramIdent" && eType != OLAS::COMPLEX ) {
      eType = OLAS::COMPLEX;
      Info->PrintF( pdename, "Expert: Using COMPLEX as matrix entry type "
                    "for parameter identification\n" );
    }
    
    // If an eigenvalue problem is solved, the entrytype of the matrices has
    // always to be DOUBLE
    if ( analysisType == BasePDE::EIGENFREQUENCY && eType != OLAS::DOUBLE ) {
      eType = OLAS::COMPLEX;
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
      if ( pType == OLAS::ILUK || pType == OLAS::ILDLK ) {
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
      else if ( ( pType == OLAS::ILU0 || pType == OLAS::JACOBI )
                && rType != OLAS::NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = OLAS::NOREORDERING;
      }
    }
  }

}
