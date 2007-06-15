// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include <limits>

#include <def_use_lapack.hh>

#include "algsys/olasparams.hh"
#include "utils/utils.hh"

namespace OLAS {
  /**********************************************************/
  
  OLAS_Params::OLAS_Params() {
    SetDefaultParams();
  }
  
  /**********************************************************/

  void OLAS_Params::SetDefaultParams() {

    ENTER_FCN( "OLAS_Params::SetDefaultParams" );

    // Do we want parallel data structures?
    SetValue( "Parallel", false );

    // By default we assume no penalty formulation
    SetValue( "UsingPenaltyFormulation", false );
    SetValue( "RHSwithPenalty", false );

    // in general we want nonoverlapping partitionings
    SetValue( "PartitionOverlap", (Integer)0);
    
    // damping parameter for Richardson iteration
    SetValue( "R_Omega", (Double)1.0);
    
    // relaxation parameter for SSOR preconditioner
    SetValue( "SSOR_Omega", (Double)1.7);
        
    // AMG
    SetValue( "AMG_CoarseningAlpha",     (Double)0.15 );
    SetValue( "AMG_CoarseningGamma",     (Integer)(-1) );
    SetValue( "AMG_NumPreSmoothing",     (Integer)1 );
    SetValue( "AMG_NumPostSmoothing",    (Integer)1 );
    SetValue( "AMG_CycleParameter",      (Integer)1 );
    SetValue( "AMG_NumBadCoarsenings",   (Integer)5 );
    SetValue( "AMG_MinSystemSize",       (Integer)20 );
    SetValue( "AMG_MaxCoarseDependency", (Integer)10 );
    SetValue( "AMG_GaussSeidelOmega",    (Double)1.0 );
    SetValue( "AMG_JacobiOmega",         (Double)0.8 );
    SetValue( "AMG_StrongDiagRatio",     (Double)0.0 );
    SetValue( "AMG_ForceFineRatio",      (Double)1e-10 );
    SetValue( "AMG_BadCoarseningRatio",  (Double)0.6 );
    SetValue( "AMG_InterpolationType",    AMG_INTERPOLATION_CONSTANT );
    SetValue( "AMG_SmootherType",         AMG_SMOOTHER_GAUSSSEIDEL );
    SetValue( "AMG_logging",              false );
    // parallel AMG
    SetValue( "AMG_BoundaryCoarsening",   false );
#ifdef USE_LAPACK
    SetValue( "AMG_DirectSolver", (SolverType)LAPACK_LU );
#else
    SetValue( "AMG_DirectSolver", (SolverType)NOSOLVER );
#endif

    // LAPACK_LU (NOTE: Changes here must be documented in LAPACK_LU!)
    SetValue( "LAPACKLU_tryScaling", true );
    SetValue( "LAPACKLU_refineSol" , true );
    SetValue( "LAPACKLU_logging"   , true );

    // LAPACK_LL (NOTE: Changes here must be documented in LAPACK_LL!)
    SetValue( "LAPACKLL_logging"         , true  );
    SetValue( "LAPACKLL_newMatrixPattern", false );

    // CG (NOTE: Changes here must be documented in CGSolver!)
    SetValue( "CG_epsilon", 1e-6 );
    SetValue( "CG_epsmach", std::numeric_limits<double>::epsilon() );
    SetValue( "CG_maxIter",  50  );
    SetValue( "CG_logging", true );
    SetValue( "CG_resDirectly", 50 );

    // GMRES (NOTE: Changes here must be documented in GMRESSolver!)
    SetValue( "GMRES_maxKrylovDim",  50 );
    SetValue( "GMRES_maxIter",        1 );
    SetValue( "GMRES_epsilon",     1e-6 );
    SetValue( "GMRES_lanczos",    false );
    SetValue( "GMRES_logging",     true );

    // MINRES (NOTE: Changes here must be documented in MINRESSolver!)
    SetValue( "MINRES_maxIter",      50 );
    SetValue( "MINRES_epsilon",    1e-6 );
    SetValue( "MINRES_logging",    true );

    // PARDISO (NOTE: Changes here must be documented in PardisoSolver!)
    SetValue( "PARDISO_posDef"      , false );
    SetValue( "PARDISO_hermitian"   , false );
    SetValue( "PARDISO_symStructure", false ); 
    SetValue( "PARDISO_ordering"    , NESTED_DISSECTION );
    SetValue( "PARDISO_logging"     , false );
    SetValue( "PARDISO_stats"       , false );
    SetValue( "PARDISO_pivoting"    , 1     );

    // MINRES (NOTE: Changes here must be documented in MINRESSolver!)
    SetValue( "MINRES_maxIter",      50 );
    SetValue( "MINRES_epsilon",    1e-6 );
    SetValue( "MINRES_logging",    true );

    // ILUTP (NOTE: Changes here must be documented in ILUTP_Precond!)
    SetValue( "ILUTP_tau"    , 0.01 );
    SetValue( "ILUTP_fillVal",   -2 );

    // ILUK (NOTE: Changes here must be documented in ILUK_Precond!)
    SetValue( "ILUK_level", 1 );
    SetValue( "ILUK_logging", true );

    // LDL_SOLVER (NOTE: Changes here must be documented in LDLSolver!)
    SetValue( "LDLSOLVER_logging"        , true  );
    SetValue( "LDLSOLVER_saveFacToFile"  , false );
    SetValue( "LDLSOLVER_savePatternOnly", false );
    SetValue( "LDLSOLVER_itRefSteps"     ,     2 );
    SetValue( "LDLSOLVER_itRefVerbosity" ,     2 );

    // LU_SOLVER (NOTE: Changes here must be documented in LUSolver!)
    SetValue( "LUSOLVER_logging"        , true );
    SetValue( "LUSOLVER_itRefSteps"     ,    2 );
    SetValue( "LUSOLVER_itRefVerbosity" ,    2 );

    // ILDL preconditioners
    // (NOTE: Changes here must be documented in ILDLPrecond!)
    SetValue( "ILDLPRECOND_logging"        ,     1 );
    SetValue( "ILDLPRECOND_saveFacToFile"  , false );
    SetValue( "ILDLPRECOND_savePatternOnly", false );
    SetValue( "ILDLPRECOND_saveLevels"     , false );
    SetValue( "ILDLPRECOND_level"          , 1     );
    SetValue( "ILDLPRECOND_tau"            , 0.01  );
    SetValue( "ILDLPRECOND_fillVal"        , 2     );

    // ARPACK Eigenvalue solver
    // (NOTE: Since arpack has internally some complicated mechanisms
    // to determine default values for the number of vectors and the tolerance,
    // we define here either -1 or -1.0 as defualt value to indicate that
    // we want arpack to choose the default value!)
    SetValue( "ARPACK_tolerance", (Double) -1.0);
    SetValue( "ARPACK_maxIt"    , (Integer) -1);
    SetValue( "ARPACK_numVec"   , (Integer) -1);
    SetValue( "ARPACK_which"    , std::string("LM") );
    SetValue( "ARPACK_logging"  , false);
    
    // By default we do not export the linear system to a file
    SetValue( "exportLinSys", false );

    // By default we use a relative stopping criterion based
    // on the norm of the right-hand side of the problem
    SetValue( "StoppingCriterion", RELNORM_RHS );

    // By default we do not export the factorisation matrix
    // computed by a Crout object (used in an ILU preconditioner
    // or the LUSolver) to a file
    SetValue( "CROUT_saveFacToFile", false );

    // By default we do not try to re-order the graph
    SetValue( "GRAPH_reordering", NOREORDERING );

    // By default we assume that the matrix sparsity pattern
    // remains constant across all solves
    SetValue( "newMatrixPattern", false );

    // When generating a SBM_System we assume by default that
    // the Finite-Element matrices will not be symmetric
    SetValue( "SBM_Symmetry", false );
  }

  /**********************************************************/
}
