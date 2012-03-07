// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <ostream>

#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/domain.hh"
#include "Driver/assemble.hh"
#include "Driver/singleDriver.hh"
#include "General/Enum.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "def_use_arpack.hh"
#include "def_use_cholmod.hh"
#include "def_use_ilupack.hh"
#include "def_use_lapack.hh"
#include "def_use_metis.hh"
#include "def_use_pardiso.hh"
#include "def_use_pardiso64.hh"

namespace CoupledField {

  // *************
  //   SetParams
  // *************
  void CFSOLASParams::SetParams( std::string pdename,
                                 PtrParamNode cfs,
                                 BasePDE::AnalysisType analysisType,
                                 Assemble * assemble,
                                 bool overrideExpert ) {

    // Fetch setup node
    PtrParamNode setupNode  = cfs->Get("setup", ParamNode::INSERT );


    // First determine approach for handling inhomogeneous Dirichlet
    // boundary conditions
    // bool usingPenalty = true; // TODO: Unused variable usingPenalty
    // {
    //  std::string aux = "penalty";
    //  setupNode->GetValue("idbcHandling", aux, ParamNode::INSERT );
    //  usingPenalty = aux == "penalty" ? true : false;
    // }

    // Determine the type of solver for this PDE
    std::string sTypeString = "expertsChoice";
    BaseSolver::SolverType sType;
    PtrParamNode solverNode;

    solverNode = cfs->Get("solver", ParamNode::INSERT);
    solverNode->GetValue("type", sTypeString, ParamNode::INSERT);
    if ( sTypeString == "expertsChoice" ) {
      if ( overrideExpert ) {
        EXCEPTION( "You cannot specify expertsChoice as solver type "
                   << "and at the same time specify overrideExpert! "
                   << "This would leave the solver type undefined!" );
      }
      sTypeString = "noSolver";    
    }
    sType = BaseSolver::solverType.Parse(sTypeString);

    // Now determine the type of preconditioner for this PDE
    std::string pTypeString = "noPrecond";
    BasePrecond::PrecondType pType;

    solverNode->GetValue("precond", pTypeString, ParamNode::INSERT );
    pType = BasePrecond::precondType.Parse(pTypeString);

    // Now determine, if we are running an eigenfrequency analysis and if
    // yes, get the type of eigenvalue solver
    std::string esTypeString = "expertsChoice";
    BaseEigenSolver::EigenSolverType esType = BaseEigenSolver::NOEIGENSOLVER;
    
    if ( analysisType == BasePDE::EIGENFREQUENCY ) {
      PtrParamNode eigenSolverNode = cfs->Get("eigenSolver", ParamNode::INSERT );
      eigenSolverNode->GetValue("type", esTypeString, ParamNode::INSERT );
      
      if ( esTypeString == "expertsChoice" ) {
        if ( overrideExpert ) {
          EXCEPTION( "You cannot specify expertsChoice as EigenSolver type "
                     << "and at the same time specify overrideExpert! "
                     << "This would leave the solver type undefined!" );
        }
        esTypeString = "noEigenSolver";
      }

      esType = BaseEigenSolver::eigenSolverType.Parse(esTypeString);
    }

    // Before we try to determine matrix properties, we first should
    // find out, what type of linear system we are using. We do this
    // by trying to determine the symmetry attribute of the sbmMatrix
    // element
    bool sbmSymmetry = false;
    bool stdSystem = true;
    PtrParamNode matrixNode;

    if(cfs->Has("sbmMatrix") ) {
      matrixNode = cfs->Get("sbmMatrix");
      stdSystem = false;
      matrixNode->GetValue("symmetric", sbmSymmetry, ParamNode::INSERT );
    }

    if( !matrixNode ) {
      matrixNode = cfs->Get("matrix", ParamNode::INSERT);
    }

    // Determine matrix entry type
    BaseMatrix::EntryType eType;
    std::string eMatString = "double";
    matrixNode->GetValue("entry", eMatString, ParamNode::INSERT);
    eType = BaseMatrix::entryType.Parse(eMatString);  
    
    // Determine calculation of condition number
    bool calcCondition = false;
    matrixNode->GetValue("calcConditionNumber", calcCondition, ParamNode::INSERT );

    // Following stuff only for required for standard systems
    BaseMatrix::StorageType mType  = BaseMatrix::NOSTORAGETYPE;
    BaseOrdering::ReorderingType orderType = BaseOrdering::NOREORDERING;
    bool allowChangeOfReordering = false;


    if ( stdSystem == true ) {
      std::string mMatString = "expertsChoice";
      // Next determine the matrix type for this PDE
      matrixNode->GetValue( "storage", mMatString, ParamNode::INSERT );
      
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
      matrixNode->GetValue( "reordering", orderString, ParamNode::INSERT );
      
      if ( orderString == "expertsChoice" ) {
        if ( overrideExpert ) {
          EXCEPTION( "You cannot specify expertsChoice as re-ordering "
                     << "strategy and at the same time specify overrideExpert! "
                     << "This would leave the re-ordering strategy undefined!" );
        }
        orderString = "noReordering";
        allowChangeOfReordering = true;
      }
      orderType = BaseOrdering::reorderingType.Parse( orderString );
    }

    // Let expert module modify the settings
    if ( overrideExpert == false && stdSystem == true ) {
      CFSOLASParams::Expert( cfs, pdename, esType, sType, pType, mType, eType,
                             orderType, analysisType, assemble,
                             allowChangeOfReordering );
    }

    // Set special parameters for solver and preconditioner
    CFSOLASParams::SetSolverParams( pdename, cfs, sType );
    CFSOLASParams::SetPrecondParams( pdename, cfs, pType );
    CFSOLASParams::SetEigenSolverParams( pdename, cfs, esType );

  }


  // *******************
  //   SetSolverParams
  // *******************
  void CFSOLASParams::SetSolverParams( std::string pdename,
                                       PtrParamNode cfs,
                                       BaseSolver::SolverType sType ) {
    // Now set the stopping rule
    std::string stopCrit = "relNormRes0";
    PtrParamNode bsNode = cfs->Get("solver", ParamNode::INSERT );
    PtrParamNode stopRuleNode = bsNode->Get("stoppingRule", ParamNode::INSERT );
    stopRuleNode->GetValue("type", stopCrit, ParamNode::INSERT);
    StopCritType stopRule;
    String2Enum( stopCrit, stopRule );
  }


  // ********************
  //   SetPrecondParams
  // ********************
  void CFSOLASParams::SetPrecondParams( std::string pdename,
                                        PtrParamNode cfs,
                                        BasePrecond::PrecondType pType ) {
  }

  // *******************
  //   SetEigenSolverParams
  // *******************
  void CFSOLASParams::SetEigenSolverParams( std::string pdename,
                                            PtrParamNode cfs,
                                            BaseEigenSolver::EigenSolverType sType ) {
  }


  // *******************
  //   Expert's Choice
  // *******************
  void CFSOLASParams::Expert( PtrParamNode cfs,
                              std::string pdename,
                              BaseEigenSolver::EigenSolverType &esType,
                              BaseSolver::SolverType &sType,
                              BasePrecond::PrecondType &pType,
                              BaseMatrix::StorageType &mType,
                              BaseMatrix::EntryType &eType,
                              BaseOrdering::ReorderingType &rType,
                              BasePDE::AnalysisType analysisType,
                              Assemble * assemble,
                              bool allowChangeOfReordering ) {


    std::string warn;
    std::string tmp;

    // ===========================
    //  Determine Matrix symmetry
    // ===========================
    bool symmetricMat = assemble->IsFEMatSymmetric();

    // ==============
    //  EigenSolver stuff
    // ==============
    // If no eigenvalue solver was specified, use arpack
    if ( esType == BaseEigenSolver::NOEIGENSOLVER ) {
      esType = BaseEigenSolver::ARPACK;
    }
    if ( analysisType == BasePDE::EIGENFREQUENCY) {
      tmp = BaseEigenSolver::eigenSolverType.ToString(esType);
      cfs->Get("eigenSolver")->Get("type")->SetValue(tmp);
    }

    // ==============
    //  Solver stuff
    // ==============

    // If no solver was specified use a direct one
    // Here we have not to determine, if we are faced with
    // a symmetric or non-symmetric problem.

    if ( sType == BaseSolver::NOSOLVER ) {
      if( symmetricMat ) {
        sType = BaseSolver::LDL_SOLVER;
      } else {
        sType = BaseSolver::LU_SOLVER;
      }
    } else {
      // In this case, the user has specified a symmetric
      // solver explicitly, although, the system is non-symmetric.
      // In this case we issue a warning
      if( sType == BaseSolver::LDL_SOLVER
          && !symmetricMat ) {
        warn = "Expert: Although the system matrix is non-symmetric ";
        warn += "a symmetric solver was chosen!";
        WARN(warn);
      }
    }
    tmp = BaseSolver::solverType.ToString(sType);
    cfs->Get("solver")->Get("type")->SetValue(tmp);

    // =======================
    //  Preconditioner stuff
    // =======================

    // For direct solver we need no preconditioner (ID)
    // PARDISO is either a direct solver or uses its own preconditioners.
    if ( (sType == BaseSolver::LDL_SOLVER
          || sType == BaseSolver::LU_SOLVER
          || sType == BaseSolver::LAPACK_LU
          || sType == BaseSolver::LAPACK_LL
          || sType == BaseSolver::PARDISO )
         && !(pType == BasePrecond::ID
              || pType == BasePrecond::NOPRECOND) ) {
      std::string solverStr = BaseSolver::solverType.ToString(sType);
      pType = BasePrecond::NOPRECOND;
      std::string precondStr = BasePrecond::precondType.ToString(pType);

      WARN("Expert: Re-setting preconditioner type to '" << precondStr
           << "' for solver '" << solverStr << "'.");
    }
    tmp = BasePrecond::precondType.ToString(pType);
    cfs->Get("solver")->Get("precond")->SetValue(tmp);

    // ===============
    //  Matrix stuff
    // ===============

    // Lapack solvers want their own matrix format
    if ( sType == BaseSolver::LAPACK_LU ) {
      if ( mType != BaseMatrix::LAPACK_GBMATRIX ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          warn = "Expert: Re-setting matrix storage type to ";
          warn += "'LAPACK_GBMATRIX'";
          WARN( warn );
        }
        else {
          Info->PrintF( pdename,
                        "Expert: Using LAPACK_GBMATRIX as storage type\n" );
        }
        mType = BaseMatrix::LAPACK_GBMATRIX;
      }
    }

    // The direct solver LU_SOLVER expects a CRS matrix
    else if ( sType == BaseSolver::LU_SOLVER ) {
      if ( mType != BaseMatrix::SPARSE_NONSYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp = BaseSolver::solverType.ToString(sType);

          WARN("Expert: Changing matrix storage type from "
               << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_NONSYM for "
               << tmp << " solver");
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for direct solver\n" );
        }
        mType = BaseMatrix::SPARSE_NONSYM;
      }
    }

    // The direct solver LDL_SOLVER expects an SCRS matrix
    else if ( sType == BaseSolver::LDL_SOLVER ) {
      if ( mType != BaseMatrix::SPARSE_SYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp = BaseSolver::solverType.ToString(sType);

          WARN("Expert: Changing matrix storage type from "
               << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_SYM for "
               << tmp << " solver");
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_SYM as storage "
                        "type for directLDL solver\n" );
        }
        mType = BaseMatrix::SPARSE_SYM;
      }
    }


    // ILU-type preconditioners expect CRS matrices
    if ( pType == BasePrecond::ILU0 || pType == BasePrecond::ILUK ) {
      if ( mType != BaseMatrix::SPARSE_NONSYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          tmp = BasePrecond::precondType.ToString(pType);

          WARN("Expert: Changing matrix storage type from "
               << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_NONSYM for "
               << tmp << " preconditioner");
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_NONSYM as storage "
                        "type for preconditioner\n" );
        }
        mType = BaseMatrix::SPARSE_NONSYM;
      }
    }

    // ILDL-type preconditioners expect SCRS matrices
    if ( pType == BasePrecond::ILDLK ) {
      if ( mType != BaseMatrix::SPARSE_SYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          tmp = BasePrecond::precondType.ToString(pType);

          WARN("Expert: Changing matrix storage type from "
               << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_SYM for "
               << tmp << " preconditioner");
        }
        else {
          Info->PrintF( pdename, "Expert: Using SPARSE_SYM as storage "
                        "type for preconditioner\n" );
        }
        mType = BaseMatrix::SPARSE_SYM;
      }
    }

    // The MG preconditioner requires a CRS matrix
    if ( pType == BasePrecond::MG ) {
      if ( mType != BaseMatrix::SPARSE_NONSYM ) {
        if ( mType != BaseMatrix::NOSTORAGETYPE ) {
          std::string tmp;
          tmp = BasePrecond::precondType.ToString(pType);

          WARN("Expert: Changing matrix storage type from "
               << BaseMatrix::storageType.ToString( mType ) << " to SPARSE_NONSYM for "
               << tmp << " preconditioner");
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
    param->GetByVal("sequenceStep", std::string("index"), seqStep)->Get("analysis")
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
      if ( sType == BaseSolver::LU_SOLVER || sType == BaseSolver::LDL_SOLVER ) {
        if ( rType == BaseOrdering::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'METIS'\n" );
          rType = BaseOrdering::METIS;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to METIS\n" );
          rType = BaseOrdering::METIS;
        }
      }
#else
      if ( sType == BaseSolver::LU_SOLVER || sType == BaseSolver::LDL_SOLVER ) {
        if ( rType == BaseOrdering::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = BaseOrdering::SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = BaseOrdering::SLOAN;
        }
      }
#endif

      // For the LAPACK solvers use SLOAN re-ordering
      if ( sType == BaseSolver::LAPACK_LU || sType == BaseSolver::LAPACK_LL ) {
        if ( rType == BaseOrdering::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = BaseOrdering::SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = BaseOrdering::SLOAN;
        }
      }

      // For Pardiso do not re-order (since Pardiso does this better itself)
      if ( sType == BaseSolver::PARDISO && rType != BaseOrdering::NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = BaseOrdering::NOREORDERING;
      }


      // ================================
      //  Reordering for Preconditioners
      // ================================

      // For all advanced ILU type preconditioners use SLOAN re-ordering
      if ( pType == BasePrecond::ILUK || pType == BasePrecond::ILDLK ) {
        if ( rType == BaseOrdering::NOREORDERING ) {
          Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                        "'SLOAN'\n" );
          rType = BaseOrdering::SLOAN;
        }
        else {
          Info->PrintF( pdename, "Expert: Re-setting re-ordering strategy "
                        "to SLOAN\n" );
          rType = BaseOrdering::SLOAN;
        }
      }

      // For ILU0 and JACOBI we do not use re-ordering.
      else if ( ( pType == BasePrecond::ILU0 || pType == BasePrecond::JACOBI )
                && rType != BaseOrdering::NOREORDERING ) {
        Info->PrintF( pdename, "Expert: Setting re-ordering strategy to "
                      "'NOREORDERING'\n" );
        rType = BaseOrdering::NOREORDERING;
      }
    }
    
    // Now that we know what solver the user wants we check if our executable
    // has been built with support for it.
    
    switch(sType) {
    case BaseSolver::PARDISO:
#ifndef USE_PARDISO
      EXCEPTION("This executable has not been built with support for PARDISO!\n"
                "Please switch on USE_PARDISO and set CFS_PARDISO.");
#endif
      break;
    case BaseSolver::PARDISO64:
#ifndef USE_PARDISO64
      EXCEPTION("This executable has not been built with support for PARDISO64!\n"
                "Please switch on USE_PARDISO64 and set CFS_PARDISO.");
#endif
      break;
    case BaseSolver::ILUPACK:
#ifndef USE_ILUPACK
      EXCEPTION("This executable has not been built with support for ILUPACK!\n"
                "Please switch on USE_ILUPACK.");
#endif
      break;
    case BaseSolver::LAPACK_LL:
    case BaseSolver::LAPACK_LU:
#ifndef USE_LAPACK
      EXCEPTION("This executable has not been built with support for LAPACK!\n"
                "Please switch on USE_LAPACK.");
#endif
      break;
    case BaseSolver::CHOLMOD:
#ifndef USE_CHOLMOD
      EXCEPTION("This executable has not been built with support for CHOLMOD!\n"
                "Please switch on USE_CHOLMOD.");
#endif      
      break;      
    default:
      break;      
    }

    switch(eType) {
    case BaseEigenSolver::ARPACK:
#ifndef USE_ARPACK
      EXCEPTION("This executable has not been built with support for ARPACK!\n"
                "Please switch on USE_ARPACK.");
#endif
      break;
    default:
      break;      
    }
    
    tmp = BaseMatrix::storageType.ToString( mType );
    cfs->Get("matrix")->Get("storage")->SetValue(tmp);

    tmp = BaseMatrix::entryType.ToString( eType );
    cfs->Get("matrix")->Get("entry")->SetValue(tmp);
    
    tmp = BaseOrdering::reorderingType.ToString( rType );
    cfs->Get("matrix")->Get("reordering")->SetValue(tmp);

  }

}
