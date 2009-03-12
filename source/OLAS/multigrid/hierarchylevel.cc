// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

%/* $Id$ */

#include <olas_use_lapack.hh>
#include <olas_use_pardiso.hh>

#include "multigrid/hierarchylevel.hh"
#include "multigrid/gaussseidel.hh"
#include "multigrid/jacobi.hh"
#include "algsys/olasparams.hh"

#include "solver/generatesolver.hh"
#include "solver/basesolver.hh"
#include "precond/idprecond.hh"

#ifdef USE_LAPACK
#include "external/lapack/lapack.hh"
#endif

// avoid linking problems with icc and old g++
#if defined(__INTEL_COMPILER)||(!__GNUC_PREREQ(3,3) )
#include "matvec/crs_matrix.cc"
#include "matvec/vector.cc"
#endif

/**********************************************************
 * dependend settings
 **********************************************************/
#ifdef  AMG_EXPORT_HIERARCHY
#include "exporttools.hh"

//// In case we want to export hierarchy components we define
//// a macro, that exports the the components by simply calling
//// the export method, and a macro, that exports a matrix
#ifdef  AMG_EXPORT_COMPONENTS
//// macro, that simply calls HierarchyLevel::ExportComponents
#define HIERARCHY_EXPORT_COMPONENTS ExportComponents()
//// macro, that exports a matrix (calls macro in "exporttools.hh"
#define HIERARCHY_EXPORT_MATRIX( MATRIX, NAMEFORMAT ) { \
    Integer exportformatid = MATRIX.GetNumRows(); \
    if( exportformatid < MATRIX.GetNumCols() )  exportformatid = MATRIX.GetNumCols(); \
    EXPORT_MATRIX_FN1( MATRIX, filename, AMG_EXPORT_PREFIX NAMEFORMAT, exportformatid ) \
}
#else
//// empty definitions
#define HIERARCHY_EXPORT_COMPONENTS
#define HIERARCHY_EXPORT_MATRIX( MATRIX, NAMEFORMAT )
#endif

//// In case we want to export the intermediate result vectors
//// during an AMG cycle we define a macro, that exports a vector.
#ifdef  AMG_EXPORT_CYCLE
/// macro, that exports a vector
#define HIERARCHY_EXPORT_VECTOR( VECTOR, NAMEFORMAT ) \
    EXPORT_VECTOR_FN1( VECTOR, filename, AMG_EXPORT_PREFIX NAMEFORMAT, VECTOR.GetSize() )
#else // AMG_EXPORT_CYCLE
#define HIERARCHY_EXPORT_VECTOR( VECTOR, NAMEFORMAT )
#endif // AMG_EXPORT_CYCLE

#else // EXPORT_HIERARCHY

//// define empty macros, if we don not export anything
#define HIERARCHY_EXPORT_COMPONENTS
#define HIERARCHY_EXPORT_VECTOR( VECTOR, NAMEFORMAT )
#define HIERARCHY_EXPORT_MATRIX( MATRIX, NAMEFORMAT )

#endif // EXPORT_HIERARCHY
/**********************************************************/
#ifdef AMG_TRACE_CF_SPLITTING
#define  HIERARCHY_TRACE_CF_SPLITTING \
         if( !GetLevelID() ) TraceCFSplitting();
#else
#define  HIERARCHY_TRACE_CF_SPLITTING
#endif
/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_HIERARCHYLEVEL
#define DEBUG_HIERARCHYLEVEL
#endif // DEBUG_HIERARCHYLEVEL
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR

#ifdef PROFILING
#ifdef PROFILE_HIERARCHYLEVEL
#include "utils/utils.hh"
#endif
#endif

#ifndef  LOGLINE
#define  LOGLINE  " --------------------------------------\n"
#endif
/**********************************************************/

namespace CoupledField {
/**********************************************************/

template <typename T>
HierarchyLevel<T>::HierarchyLevel()
                 : LevelID_( 0 ),
                   SysMatrix_( NULL ),
                   AuxMatrix_( NULL ),
                   DirSysMatrix_( NULL ),
                   dirSysMatPPool_( NULL ),
                   Topology_( NULL ),
                   Transfer_( NULL ),
                   Smoother_( NULL ),
                   directSolver_( NULL ),
                   directSolverParams_( NULL ),
                   vh1_( NULL ),
                   vH1_( NULL ),
                   vH2_( NULL ),
                   nextLevel_( NULL ),
#ifdef PROFILE_HIERARCHYLEVEL
                   logging_( true ),
#else
                   logging_( false ),
#endif
                   deleteDirSysMatrix_( true )
{
}

/**********************************************************/

template <typename T>
HierarchyLevel<T>::HierarchyLevel( const Integer level_id )
                 : LevelID_( level_id ),
                   SysMatrix_( NULL ),
                   AuxMatrix_( NULL ),
                   DirSysMatrix_( NULL ),
                   dirSysMatPPool_( NULL ),
                   Topology_( NULL ),
                   Transfer_( NULL ),
                   Smoother_( NULL ),
                   directSolver_( NULL ),
                   directSolverParams_( NULL ),
                   vh1_( NULL ),
                   vH1_( NULL ),
                   vH2_( NULL ),
                   nextLevel_( NULL ),
#ifdef PROFILE_HIERARCHYLEVEL
                   logging_( true ),
#else
                   logging_( false ),
#endif
                   deleteDirSysMatrix_( true )
{
}

/**********************************************************/

template <typename T>
HierarchyLevel<T>::~HierarchyLevel()
{
    Reset();
}

/**********************************************************/

template <typename T>
inline const HierarchyLevel<T>* HierarchyLevel<T>::
GetNextLevel() const
{
    return nextLevel_;
}

/**********************************************************/

template <typename T>
inline Integer HierarchyLevel<T>::GetLevelID() const
{
    return LevelID_;
}

/**********************************************************/

template <typename T>
void HierarchyLevel<T>::Reset()
{

    if( !GetLevelID() && (SysMatrix_ || AuxMatrix_) ) {
        std::cerr
        << "called \033[31mHierarchyLevel::Reset\033[0m"
        << std::endl << "on level "<<GetLevelID()<<" with still"
           " present system or auxiliary matrix" << std::endl
        << "might not be intended -> eventually crash"
        << std::endl;
    }

    delete SysMatrix_;     SysMatrix_ = NULL;
    delete AuxMatrix_;     AuxMatrix_ = NULL;

    if( deleteDirSysMatrix_ ) {
        delete DirSysMatrix_;
        DirSysMatrix_ = NULL;
    }
    deleteDirSysMatrix_ = true;
    // It is crucial to delete the pattern pool AFTER
    // DirSysMatrix_, since the destructor of DirSysMatrix_
    // deregisters itself at the pool. Otherwise the pool
    // would not free the pattern memory.
    delete dirSysMatPPool_;  dirSysMatPPool_ = NULL;

    delete Topology_;      Topology_     = NULL;
    delete Transfer_;      Transfer_     = NULL;
    delete Smoother_;      Smoother_     = NULL;
    delete directSolver_;  directSolver_ = NULL;
    delete directSolverParams_;  directSolverParams_ = NULL;

    // delete temporary vectors
    delete vh1_;  vh1_ = NULL;
    delete vH1_;  vH1_ = NULL;
    delete vH2_;  vH2_ = NULL;

    // remove all coarser levels
    delete nextLevel_;  nextLevel_ = NULL;
}

/**********************************************************/

template <typename T>
void HierarchyLevel<T>::
InsertSysMatrix( const CRS_Matrix<T>* const sysmat )
{

    delete SysMatrix_;
    SysMatrix_ = (CRS_Matrix<T>*)sysmat;
}


template <typename T>
void HierarchyLevel<T>::
InsertAuxMatrix( const CRS_Matrix<T>* const auxmat )
{
    
    delete AuxMatrix_;
    AuxMatrix_ = (CRS_Matrix<T>*)auxmat;
}

/**********************************************************/

template <typename T>
const CRS_Matrix<T>* HierarchyLevel<T>::GetSysMatrixPtr() const
{
    
    return (const CRS_Matrix<T>*) SysMatrix_;
}


template <typename T>
const CRS_Matrix<T>* HierarchyLevel<T>::GetAuxMatrixPtr() const
{
    
    return (const CRS_Matrix<T>*) AuxMatrix_;
}

/**********************************************************/

template <typename T>
CRS_Matrix<T>* HierarchyLevel<T>::UnhookSysMatrix()
{
    
    CRS_Matrix<T> *temp = SysMatrix_;
    SysMatrix_ = NULL;
    return temp;
}


template <typename T>
CRS_Matrix<T>* HierarchyLevel<T>::UnhookAuxMatrix()
{
    
    CRS_Matrix<T> *temp = AuxMatrix_;
    AuxMatrix_ = NULL;
    return temp;
}

/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::Setup( Settings* const settings,
                               Integer         numBadCoarsenings )
{

    bool *penaltyFlags = NULL;

    if( SysMatrix_ == NULL )  return false;
    // (de)activate logging
    SetLogging( settings->logging );

#ifdef PROFILE_HIERARCHYLEVEL
    Double t1 = AMG_GET_REAL_TIME
#endif // PROFILE_HIERARCHYLEVEL

 //////////////////// 
 // setup THIS level

    // some logging information
    if( logging_ ) {
        if( ! GetLevelID() ) {
            (*cla) <<LOGLINE<<" AMG: SETTINGS:\n";
            settings->Print( *cla );
            (*cla) << LOGLINE;
        }
        (*cla) << " AMG: starting setup on level ["<<GetLevelID()
               << "]"<<std::endl;
    }

#ifdef AMG_DIRICHLET_MIXED_SMOOTHING
    if( !GetLevelID() ) {
        NEWARRAY( penaltyFlags, bool, SysMatrix_->GetNumRows() ); }
#endif // AMG_DIRICHLET_MIXED_SMOOTHING

    // setup the topology (also creates the dependency graphs)
    if( ! SetupTopology(settings->alpha,
                        settings->strongDiagRatio,
                        settings->forceFineRatio,
                        penaltyFlags) ) {
        DELETEARRAY( penaltyFlags );
        return false;
    }

    // compute the coarse-fine splitting.
    const Integer SizeH =
    Topology_->CalcCoarseFineSplitting( settings->maxCoarseDepend,
                                        settings->gamma );
    if( logging_ ) {
        (*cla)
        << " AMG: coarsening nh:"<<Topology_->GetSizeh()
        << " -> nH:"<<SizeH<<std::endl;
    }

    // create auxiliary data. Here we need the sizes of both, the
    // fine and the coarse system (this is the reason, why we stored it)
    // In case SizeH > settings->minSystemSize, this level will stay
    // the bottom level, and we do not need the temporary vectors of
    // size nH. We prevent their creation by passing 0 as second
    // parameter to SetupAuxData.
    if( ! SetupAuxData(
              Topology_->GetSizeh(),
              SizeH < settings->minSystemSize ? 0 : SizeH) ) {
        DELETEARRAY( penaltyFlags );
        return false;
    }
    // check, if the coarsenig gets slow
    if( settings->badCoarseningRatio * Topology_->GetSizeh() <= SizeH ) {
        numBadCoarsenings++;
    } else {
        numBadCoarsenings = 0;
    }
    // If the number of coarse unknowns falls below the minimum, create
    // a direct solver and interrupt the setup. (note, that we still
    // keep the temporary vectors of size h).
    if( settings->minSystemSize     >  SizeH             ||
        settings->numBadCoarsenings <  numBadCoarsenings ||  
        Topology_->GetSizeh()       == SizeH                ) {
        // If logging is switched on, write some information to
        // the (*cla) stream. If profiling is switched on, write
        // additionally the profiling result. Note that profiling
        // implicitly switches on the logging.
        if( logging_ ) {
            if( settings->minSystemSize > SizeH ) {
                (*cla)
                 << " AMG: minimal system size "
                 << settings->minSystemSize<<std::endl
                 << " AMG:  => reached coarsest level of size "
                 << Topology_->GetSizeh()<<std::endl;
            // settings->numBadCoarsenings < numBadCoarsenings
            } else if( settings->numBadCoarsenings < numBadCoarsenings ) {
                (*cla)
                 << " AMG: had too many ("<<numBadCoarsenings<<") "
                    "\"bad coarsenings\""<<std::endl
                 << " AMG:  => reached coarsest level of size "
                 << Topology_->GetSizeh()<<std::endl;
            // Topology_->GetSizeh() == SizeH
            } else {
                (*cla)
                 << " AMG: nh = "<<SizeH<<" could not be reduced\n"
                 << " AMG:  => reached coarsest level\n";
            }
        }

        // remove topology object and eventually the penalty flags
        delete Topology_;  Topology_ = NULL;
        DELETEARRAY( penaltyFlags );
        // setup direct solver
        const bool rValue = SetupDirectSolver( settings );

        if( logging_ ) {
            (*cla)
            << " AMG: finished setup on level ["<<GetLevelID()<<"]";
#ifdef PROFILE_HIERARCHYLEVEL
            const Double t2 = AMG_GET_REAL_TIME
            (*cla) << ": " <<(t2-t1)<<" seconds";
#endif
            (*cla) <<std::endl<<LOGLINE;
        }
        // eventualy export the system matrix
        HIERARCHY_EXPORT_COMPONENTS;
        return rValue;
    }

    // setup the smoother
    if( ! SetupSmoother(settings, penaltyFlags) ) {
        DELETEARRAY( penaltyFlags );
        return false;
    }
    // create transfer operator
    if( NULL  == (Transfer_ = New TransferOperator<T>()) )  return false;
    if( ! Transfer_->CreateOperators(*SysMatrix_, *Topology_,
                                     settings->InterpolationType,
                                     true) ) {
        return false;
    }

    // create coarse matrix
    CRS_Matrix<T> *coarseMatrix = NULL;

#ifdef HIERARCHYLEVEL_CALC_COARSE_GRAPHS
    // create the graphs for the galerkin product
    DependencyGraph<T> *graph_AHT = New DependencyGraph<T>(),
                       *graph_VT  = New DependencyGraph<T>();
    Topology_->CalcGalerkinGraphs( *graph_AHT, *graph_VT );

    // now we do not need the topology any more
    if( false == settings->keepTopology ) {
        delete Topology_;  Topology_ = NULL; }

    // calculate the galerkin product directly
    Transfer_->GalerkinProduct( &coarseMatrix, *SysMatrix_,
                                *graph_AHT   , *graph_VT    );

    // now we do not need the graphs any more
    delete graph_AHT;  graph_AHT = NULL;
    delete graph_VT;   graph_VT  = NULL;
#else
    // now we do not need the topology any more
    if( false == settings->keepTopology ) {
        delete Topology_;  Topology_ = NULL; }

    // calc the Galerkin product using a temporary PreMatrix
    Transfer_->GalerkinProduct( &coarseMatrix, *SysMatrix_ );
#endif

    // some logging and profiling
#ifdef PROFILE_HIERARCHYLEVEL
    Double t2 = AMG_GET_REAL_TIME
#endif
    if( logging_ ) {
        (*cla) << " AMG: finished setup on level ["<<GetLevelID()<<"]:\n"
#ifdef PROFILE_HIERARCHYLEVEL
               << " AMG: setup on level ["<<GetLevelID()<<"] "
                  "in "<<(t2-t1)<<" seconds\n"
#endif
               << LOGLINE;
    }

 ////////////////////////////////////////////////////
 // continue setup process recursively on next level
    // create level object for the next level
    delete nextLevel_;
    nextLevel_ = New HierarchyLevel( GetLevelID()+1 );
    nextLevel_->SetLogging( logging_ );
    nextLevel_->InsertSysMatrix( coarseMatrix );
    const bool setupResult = nextLevel_->Setup( settings,
                                                numBadCoarsenings );
    // eventual export of the set up components
    HIERARCHY_EXPORT_COMPONENTS;
    // eventual trace CF-splitting
    HIERARCHY_TRACE_CF_SPLITTING;

    return setupResult;
}

/**********************************************************/

template <typename T>
void HierarchyLevel<T>::
Cycle( const Vector<T>& rhs,
             Vector<T>& sol,
       const Integer    num_presmoothings,
       const Integer    num_postsmoothings,
       const Integer    cycle_parameter     )
{

#ifdef AMG_EXPORT_CYCLE
    char filename[256];
#endif
    const Integer maxNumSmoothingSteps = 100;

// eventual export
HIERARCHY_EXPORT_VECTOR( rhs, "b_%d.00" );
if( nextLevel_ )  HIERARCHY_EXPORT_VECTOR( sol, "u_%d.00" );

    // if there exists a coarser level, apply the usual
    // AMG-circus, i.e. smoothing, restriction, ...
    if( nextLevel_ ) {
        // If this is the top level, we use the passed initial
        // solution. On all coarser levels we start with the zero
        // initial guess.
        if( GetLevelID() > 0 ) {
            sol.Init();
// eventual export
HIERARCHY_EXPORT_VECTOR( sol, "u_%d.00" );
        }

        // use the cycle parameter only on coarser levels
        const Integer lambda = GetLevelID() == 0 ? 1 : cycle_parameter;
        for( Integer i = 0; i < lambda; i++ ) {
            // v1 presmoothing steps (forward, suppress new setup)
            for( Integer v1 = 0; v1 < num_presmoothings; v1++ ) {
                Smoother_->Step( *SysMatrix_, rhs, sol,
                                 Smoother<T>::FORWARD, false );
            }
            // calculate residual (vh1_ = r = b - Ax) and restrict it
            SysMatrix_->CompRes( *vh1_, sol, rhs );

// eventual export
HIERARCHY_EXPORT_VECTOR( sol,    "u_%d.01" );
HIERARCHY_EXPORT_VECTOR((*vh1_), "r_%d.00" );

            Transfer_->Restrict( *vh1_, *vH1_, false );
            // AMG cycle on coarser level
            nextLevel_->Cycle( *vH1_, *vH2_,
                               num_presmoothings, num_postsmoothings,
                               cycle_parameter );
            // interpolate coarse correction, including addition
            Transfer_->Prolongate( *vH2_, sol, true );

// eventual export
HIERARCHY_EXPORT_VECTOR( sol, "u_%d.02" );

            // v2 postsmoothing steps (backward, suppress new setup)
            for( Integer v2 = 0; v2 < num_postsmoothings; v2++ ) {
                Smoother_->Step( *SysMatrix_, rhs, sol,
                                 Smoother<T>::BACKWARD, false );
            }
            
// eventual export
HIERARCHY_EXPORT_VECTOR( sol, "u_%d.04" );

        }
    // if there does not exist a coarser level, we must
    // solve exactly on this level
    } else {
        if( directSolver_ ) {
            IdPrecond idprecond;
            directSolver_->Solve( *DirSysMatrix_, idprecond, rhs, sol );
        } else {
            if( logging_ ) {
                (*cla) << " AMG: using smoother as direct solver, "
                       << "maximal "<<maxNumSmoothingSteps<<" steps"
                       << std::endl;
            }
            Integer nExactSmoothingSteps = 0;
            Double residualL2Norm = 0;
            typename Smoother<T>::Direction dir = Smoother<T>::FORWARD;
            sol.Init();
            do{ if( ++nExactSmoothingSteps > maxNumSmoothingSteps ) {
                    Warning( __FILE__, __LINE__,
                             "AMG: Used smoother as exact solver on coarsest "
                             "level, but did not succeed within %d steps "
                             "-> aborting exact solving, choose a better exact "
                             "solver by setting AMG_DirectSolver",
                             maxNumSmoothingSteps );
                    break;
                }
                Smoother_->Step( *SysMatrix_, rhs, sol, dir, false );
                // flip direction
                dir = (dir == Smoother<T>::FORWARD ?
                       Smoother<T>::BACKWARD : Smoother<T>::FORWARD);
                SysMatrix_->CompRes( *vh1_, sol, rhs );
                residualL2Norm = vh1_->NormL2();
                if( logging_ ) {
                    (*cla) << " AMG: smoother: step "<<nExactSmoothingSteps
                           << " : ||r|| = "<<std::scientific<<residualL2Norm
                           << std::endl;
                }
            } while( residualL2Norm > 1e-15 );
        }

// eventual export of exact coarse solution
HIERARCHY_EXPORT_VECTOR( sol, "u_%d.04" );

    }
}

/**********************************************************/

template <typename T>
Double HierarchyLevel<T>::GetOperatorComplexity() const
{

    const HierarchyLevel<T>* level      = this;
    Double                   complexity = 0.0;
    
    while( level ) {
        complexity += level->SysMatrix_->GetNnz();
        level = level->GetNextLevel();
    }
    complexity /= SysMatrix_->GetNnz();

    return complexity;
}

/**********************************************************/

template <typename T>
std::ostream& HierarchyLevel<T>::Print( std::ostream& out,
                                        const bool color ) const
{

    if( SysMatrix_ == NULL ) {
        out << "no matrix present in hierarchy level "
            << LevelID_ << std::endl;
        return out;
    }
    
    out
    << (color ? "\033[0;1m" : "")
    << "hierarchy level "<<LevelID_
    << (color ? "\033[0m" : "")<<std::endl
    << "  system matrix : size : " << SysMatrix_->GetNumRows() << 'x'
                                   << SysMatrix_->GetNumCols() << std::endl
    << "                : nnz  : " << SysMatrix_->GetNnz()   << std::endl
    << "                : dens : " << (((double)SysMatrix_->GetNnz()  /
                                        (double)SysMatrix_->GetNumRows() )/
                                        (double)SysMatrix_->GetNumRows()   )
    << std::endl
    << "  complexity    : " << GetOperatorComplexity() << std::endl;

    if( nextLevel_ )  nextLevel_->Print( out );

    return out;
}

/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::
SetupTopology( const Double alpha,
               const Double strongDiagRatio,
               const Double forceFineRatio,
               bool *const  penaltyFlags )
{

    if( logging_ )  (*cla) << " AMG: setup topology\n";

    // Create topology object. The constructor, we use here
    // also creates the graphs of strong dependencies S and ST.
    if( Topology_ ) { delete Topology_;  Topology_ = NULL; }
    
    // create new topology object
    if( NULL == (Topology_ = New Topology<T>) ) {
        (*cla) << " AMG: creation of topology object failed\n";
        return false;
    }

    // setup dependency graphs
    if( 0 > Topology_->CreateDependencyGraphs(
                 *SysMatrix_,       // system matrix
                  alpha,            // alpha for coarsening
                  strongDiagRatio, 
                  forceFineRatio
#ifdef AMG_DIRICHLET_MIXED_SMOOTHING
                , penaltyFlags
#endif
    ) ) {
        (*cla) << " AMG: topology setup failed\n";
        return false;
    }

    return true;
}

/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::
SetupSmoother( const Settings *const settings,
               bool *const penaltyFlags )
{

    if( logging_ ) (*cla) << " AMG: setup smoother\n";

    // create smoother object
    if( Smoother_ ) { delete Smoother_; Smoother_ = 0x0; }
    switch( settings->SmootherType ) {
        case AMG_SMOOTHER_GAUSSSEIDEL: {
            GaussSeidel<T> *gss = New GaussSeidel<T>;
            if( gss ) {
                gss->SetOmega( settings->smootherOmega );
                Smoother_ = gss;
            }
            break;
        }
        case AMG_SMOOTHER_DAMPED_JACOBI: {
            Jacobi<T> *js = New Jacobi<T>;
            if( js ) {
                js->SetOmega( settings->smootherOmega );
                Smoother_ = js;
            }
            break;
        }
        default: {
            (*cla)
            << " AMG: not supported smoother type ID "
            << settings->SmootherType<<std::endl;
        }
    }

    if( ! Smoother_ ) {
        (*cla)
        << " AMG: creation of smoother object failed"<<std::endl;
        return false;
    } else {
        // setup the smoother
        if( ! Smoother_->Setup(*SysMatrix_) ) {
            (*cla) << " AMG: setup of smoother failed\n";
            delete Smoother_;
            Smoother_ = 0x0;
            return false;
        }
    }    
    
    return true;
}

/**********************************************************/
// macro for the switch cases in method SetupDirectSolver
// We need all these cases explicitly, because
//  -> GenerateStdMatrixObject returns a BaseMatrix pointer
//  -> BaseMatrix has not got a method Convert
//  -> the template parameter for the cast to LapackGBMatrix
//     cannot be specified dynamically
#define CASE_CONVERT_CRS_TO_LAPACKGB( FTypeID, entryF ) \
case FTypeID: { \
    LapackGBMatrix< entryF, T_Mtype >* lapackA = \
        dynamic_cast<LapackGBMatrix< entryF, T_Mtype >*>(DirSysMatrix_); \
    if( lapackA == NULL ) { \
        Error( "HierarchyLevel<T>::SetupDirectSolver: dynamic cast " \
               "failed", __FILE__, __LINE__ ); \
    } else { \
        lapackA->Convert( *SysMatrix_ ); \
    } \
    break; \
}
/****************************/

template <typename T>
bool HierarchyLevel<T>::
SetupDirectSolver( const Settings* const settings )
{

    if( logging_ ) (*cla)<<" AMG: setup exact solver\n";

    switch( settings->directSolver ) {

        //////////////////////////////////////////
        // NOSOLVER as setting for the direct solver in AMG means to
        // use the smoother as direct solver on the coarsest level
        //////////////////////////////////////////

        case NOSOLVER: {
            if( ! SetupSmoother(settings) )  return false;
            break;
        }

        //////////////////////////////////////////
        // use the Lapack LU decomposition
        //////////////////////////////////////////

        case LAPACK_LU: {
#ifdef USE_LAPACK
            // first sort the CRS_Matrix
            SysMatrix_->SortConformingToLayout();
            // create the lapack matrix
            DirSysMatrix_ = GenerateStdMatrixObject(
                                SysMatrix_->GetEntryType(),
                                LAPACK_GBMATRIX, 1,
                                SysMatrix_->GetNumRows(),
                                SysMatrix_->GetNumCols(),
                                SysMatrix_->GetNnz()    );

            // only LapackGBMatrix knows the Convert method
            switch( DirSysMatrix_->GetEntryType() ) {
                CASE_CONVERT_CRS_TO_LAPACKGB( F77REAL8, F77real8 );
                // CASE_CONVERT_CRS_TO_LAPACKGB( F77COMPLEX16, F77complex16 );
                default:
                    Error( "No case-statement for FORTRAN entry type",
                           __FILE__, __LINE__ );
                break;
            }

            // Create parameter object and set parameters
            // for direct solution with Lapack_LU solver.
            // The solution on the coarsest grid need not
            // be too exact so we go for a fast approach
            if ( directSolverParams_ == NULL ) {
                directSolverParams_ = New OLAS_Params;
            }
            directSolverParams_->SetValue( "LAPACKLU_logging", logging_ );

            // Now generate solver using the factory and trigger factorisation
            directSolver_ = GenerateSolverObject( *DirSysMatrix_, LAPACK_LU,
                                                   directSolverParams_, NULL );
            directSolver_->Setup( *DirSysMatrix_ );
#else // USE_LAPACK
            Error( "To use the LAPACK LU decomposition as direct solver "
                   "on the coarsest AMG level the code must be compiled "
                   "with preprocessor flag USE_LAPACK", __FILE__, __LINE__ );
#endif // USE_LAPACK
            break;
        }

        //////////////////////////////////////////
        // use the PARDISO LU decomposition
        //////////////////////////////////////////

        case PARDISO: {
#ifdef USE_PARDISO
            // sort the serial CRS matrix
            SysMatrix_->SortConformingToLayout();
            // The PARDISO solver works with CRS matrices. Therefore
            // we do not need to convert the matrix, but only point
            // to it with the pointer to the general direct solver
            // matrix.
            DirSysMatrix_ = SysMatrix_;
            // remember that DirSysMatrix_ must not be deleted
            // separately
            deleteDirSysMatrix_ = false;
            // create new parameter object for the direct solver
            if ( directSolverParams_ == NULL ) {
                directSolverParams_ = New OLAS_Params;
            }
            // >> this would be the right place to modify the
            // >> settings for the PARDISO solver's behaviour
            // generate solver object
            directSolver_ =
            GenerateSolverObject( *DirSysMatrix_, PARDISO,
                                  directSolverParams_, NULL );
            directSolver_->Setup( *DirSysMatrix_ );
#else // USE_PARDISO
            Error( "To use PARDISO as direct solver on the coarsest "
                   "AMG level the code must be compiles with the "
                   "preprocessor flags USE_PARDISO",
                   __FILE__, __LINE__ );
#endif // USE_PARDISO
            break;
        }

        //////////////////////////////////////////
        // use the LDL solver
        //////////////////////////////////////////

        case LDL_SOLVER: {
            // sort the serial CRS matrix
            SysMatrix_->SortConformingToLayout();
            // The LDL solver works with SCRS matrices.
            DirSysMatrix_ = 
            GenerateStdMatrixObject( SysMatrix_->GetEntryType(),
                                     SPARSE_SYM,
                                     BlockSize<T_Mtype>::size,
                                     SysMatrix_->GetNumRows(),
                                     SysMatrix_->GetNumCols(),
                                     SysMatrix_->GetNnz() );
            // cast matrix pointer
            PTRCAST( DirSysMatrix_, SCRS_Matrix<T>, scrs );
            // create sparsity pattern for the SCRS_Matrix
            // NOTE that the number of nonzeros is only the one of the
            //      upper diagonal part, including the diagonal.
            SCRS_Pattern *pattern = New SCRS_Pattern;
            NEWARRAY( pattern->rptr_, Integer, SysMatrix_->GetNumRows()+1 );
            NEWARRAY( pattern->cidx_, Integer,
                      (SysMatrix_->GetNnz() + SysMatrix_->GetNumRows())/2 );
            // insert the pattern into the pattern pool and set
            // the pool's pattern in the SCRS_Matrix
            dirSysMatPPool_ = New PatternPool;
            scrs->SetSparsityPattern( dirSysMatPPool_,
                                      dirSysMatPPool_->InsertPattern(pattern) );
            // remember that DirSysMatrix_ must be deleted separately
            deleteDirSysMatrix_ = true;

            // convert the system matrix into the symmetrix sparse format
            Integer *const crsRow  = SysMatrix_->GetRowPointer();
            Integer *const crsCol  = SysMatrix_->GetColPointer();
            T_Mtype *const crsDat  = SysMatrix_->GetDataPointer();
            Integer *const scrsRow = scrs->GetRowPointer();
            Integer *const scrsCol = scrs->GetColPointer();
            T_Mtype *const scrsDat = scrs->GetDataPointer();
            Integer  crs_ij = 1,
                    scrs_ij = 1;
            scrsRow[1] = scrs_ij;

            // fill in the values
            for( Integer i = 1; i <= SysMatrix_->GetNumRows(); i++ ) {
                // set diagonal entry
                scrsCol[scrs_ij  ] = i;
                scrsDat[scrs_ij++] = crsDat[crs_ij = crsRow[i]];
                // skip entries in lower diagonal part. NOTE that we
                // assume the diagonal element to be present in each
                // line.
                while( ++crs_ij < crsRow[i+1] && crsCol[crs_ij] < i );
                // fill in entries in upper diagonal part
                while( crs_ij < crsRow[i+1] ) {
                    scrsCol[scrs_ij  ] = crsCol[crs_ij  ];
                    scrsDat[scrs_ij++] = crsDat[crs_ij++];
                }
                scrsRow[i+1] = scrs_ij;
            }

            // create new parameter object for the direct solver
            if ( directSolverParams_ == NULL ) {
                directSolverParams_ = New OLAS_Params;
            }
            // >> this would be the right place to modify the
            // >> settings for the LDL solver's behaviour
            // generate solver object
            directSolver_ =
            GenerateSolverObject( *DirSysMatrix_, LDL_SOLVER,
                                  directSolverParams_, NULL );
            directSolver_->Setup( *DirSysMatrix_ );

            break;
        }

        //////////////////////////////////////////
        // not supported solver type
        //////////////////////////////////////////

        default: {
              Error( "Solver type not supported as direct solver in AMG",
                     __FILE__, __LINE__ );
              break;
        }
    }
    
    return true;
}

/**********************************************************/
#undef CASE_CONVERT_CRS_TO_LAPACKGB
/**********************************************************/

template <typename T>
bool HierarchyLevel<T>::SetupAuxData( const Integer sizeh,
                                      const Integer sizeH )
{

    vh1_ = New Vector<T>( sizeh );

    if( sizeH ) {
        vH1_ = New Vector<T>( sizeH );
        vH2_ = New Vector<T>( sizeH );
    }

    return true;
}

/**********************************************************/
#ifdef AMG_EXPORT_COMPONENTS

template <typename T>
void HierarchyLevel<T>::ExportComponents()
{
    char filename[500];
    HIERARCHY_EXPORT_MATRIX( (*SysMatrix_), "A_%d" );
    if( nextLevel_ ) {
        HIERARCHY_EXPORT_MATRIX( (*Transfer_->GetProlongation()), "P_%d" );
        HIERARCHY_EXPORT_MATRIX( (*Transfer_->GetRestriction()),  "R_%d" );
    }
}

#endif // AMG_EXPORT_COMPONENTS
/**********************************************************/
#ifdef AMG_TRACE_CF_SPLITTING

template <typename T>
bool HierarchyLevel<T>::
TraceCFSplitting( AMG_CF_TRACE_TYPE *trace,
                  Integer           *indexMap ) const
{
    bool  traceResult   = true;
    const Integer sizeh = SysMatrix_->GetNumRows();
    
    // On top level create and initialize arrays for splitting
    // trace and for index map [index] -> [index on fines level].
    if( GetLevelID() == 0 ) {
        NEWARRAY( trace,    AMG_CF_TRACE_TYPE, sizeh );
        NEWARRAY( indexMap, Integer,           sizeh );
        // Initializing trace with -1 allows succesive increment
        // of the values in trace instead of assigning the level ID.
        // This enables correctness check.
        for( Integer i = 1; i <= sizeh; i++ )  trace[i]    = -1;
        // initialize index map on top level
        for( Integer i = 1; i <= sizeh; i++ )  indexMap[i] =  i;
    } else {
        if( !trace || !indexMap) {
            Warning( __FILE__, __LINE__,
                     "HierarchyLevel::TraceCFSplitting: (level %d) this "
                     " function should be called on coarser levels only "
                     "recursively by the next finer level", GetLevelID() );
            return false;
        }
    }

    // Increment the ID of all top level point, that are
    // also part of this level.
    for( Integer ih = 1; ih <= sizeh; ih++ ) {
        if( ++trace[indexMap[ih]] != GetLevelID() ) {
            Warning( __FILE__, __LINE__,
                     "HierarchyLevel::TraceCFSplitting: error on "
                     "level %d", GetLevelID() );
            traceResult = false;
            break;
        }
    }
    
    if( traceResult && GetNextLevel() ) {
        if( Topology_ ) {
            Integer *indexMapH = NULL;
            NEWARRAY( indexMapH, Integer, Topology_->GetSizeH() );
            for( Integer ih = 1; ih <= sizeh; ih++ ) {
                if( Topology_->IsCPoint(ih) ) {
                    indexMapH[Topology_->GetCoarseIndex(ih)] = indexMap[ih];
                }
            }
            traceResult &=
            GetNextLevel()->TraceCFSplitting( trace, indexMapH );

            DELETEARRAY( indexMapH );
        // Topology_ == NULL
        } else {
            Warning( __FILE__, __LINE__,
                     "HierarchyLevel::TraceCFSplitting: no topology "
                     "object on level %d. Cannot trace coarsening",
                     GetLevelID() );
            traceResult = false;
        }
    }

    // write splitting trace
    if( GetLevelID() == 0 ) {
        if( traceResult ) {
            const char* const filename =
                AMG_EXPORT_PREFIX   AMG_CF_SPLITTING_TRACE_FILE_NAME;
            FILE *file = fopen( filename, "wt" );
            if( file ) {
                std::cout
                << "\033[37m exporting C-F-splitting trace to "
                   "\033[0;1m\""<<filename<<"\"\033[0m\n";
#ifdef AMG_CF_TRACE_BINARY
                fwrite( trace+1, sizeof(AMG_CF_TRACE_TYPE),
                        Topology_->GetSizeh(), file );
#else
                for( Integer i = 1; i <= Topology_->GetSizeh(); i++ ) {
                    fprintf( file, "%d\n", trace[i] );
                }
#endif
                fclose( file );
            } else {
                Warning( __FILE__, __LINE__,
                         "HierarchyLevel::TraceCFSplitting: could not open "
                         "file \033[1\"%s\"\033[0m. Trace file not written",
                         filename );
                traceResult = false;
            }
        } else {
            Warning( "HierarchyLevel::TraceCFSplitting:\n since there were "
                     "errors on at least one level, no splitting trace file "
                     "is written.", __FILE__, __LINE__ );
        }

        DELETEARRAY( trace );
        DELETEARRAY( indexMap );
    }

    return traceResult;
}

#endif // AMG_TRACE_CF_SPLITTING
/**********************************************************/

/**********************************************************
 * HierarchyLevel<T>::Settings
 **********************************************************/

template <typename T>
HierarchyLevel<T>::Settings::Settings( OLAS_Params* const params )
                 : maxCoarseDepend( 10 ),
                   minSystemSize( 20 ),
                   numPreSmoothings( 1 ),
                   numPostSmoothings( 1 ),
                   CycleParameter( 1 ),
                   numBadCoarsenings( 5 ),
                   gamma( -1 ),
                   alpha( 0.15 ),
                   strongDiagRatio( 0.0 ),
                   forceFineRatio( 1e-10 ),
                   badCoarseningRatio( 0.6 ),
                   keepTopology( false ),
                   logging( false ),
                   directSolver( NOSOLVER ),
                   InterpolationType( AMG_INTERPOLATION_CONSTANT ),
                   SmootherType( AMG_SMOOTHER_GAUSSSEIDEL ),
                   smootherOmega( 1.0 )
{
    if( params ) {
        maxCoarseDepend    = params->GetIntValue( "AMG_MaxCoarseDependency" );
        minSystemSize      = params->GetIntValue( "AMG_MinSystemSize" );
        numPreSmoothings   = params->GetIntValue( "AMG_NumPreSmoothing" );
        numPostSmoothings  = params->GetIntValue( "AMG_NumPostSmoothing" );
        CycleParameter     = params->GetIntValue( "AMG_CycleParameter" );
        numBadCoarsenings  = params->GetIntValue( "AMG_NumBadCoarsenings" );
        gamma              = params->GetIntValue( "AMG_CoarseningGamma" );
        alpha              = params->GetDoubleValue( "AMG_CoarseningAlpha" );
        strongDiagRatio    = params->GetDoubleValue( "AMG_StrongDiagRatio" );
        forceFineRatio     = params->GetDoubleValue( "AMG_ForceFineRatio" );
        badCoarseningRatio = params->GetDoubleValue( "AMG_BadCoarseningRatio" );
        logging            = params->GetBoolValue( "AMG_logging" );
        params->GetEnumValue( "AMG_DirectSolver", directSolver );
        params->GetEnumValue( "AMG_InterpolationType", InterpolationType );
        params->GetEnumValue( "AMG_SmootherType", SmootherType );
        if( SmootherType == AMG_SMOOTHER_GAUSSSEIDEL ) {
            smootherOmega  = params->GetDoubleValue( "AMG_GaussSeidelOmega" );
        } else {
            smootherOmega  = params->GetDoubleValue( "AMG_JacobiOmega" );
        }
    }
}

/**********************************************************/

template <typename T>
std::ostream& HierarchyLevel<T>::Settings::
Print( std::ostream& out ) const
{
    out
    << " AMG:    MaxCoarseDependency : "<<maxCoarseDepend     <<std::endl
    << " AMG:    MinSystemSize       : "<<minSystemSize       <<std::endl
    << " AMG:    NumPreSmoothing     : "<<numPreSmoothings    <<std::endl
    << " AMG:    NumPostSmoothing    : "<<numPostSmoothings   <<std::endl
    << " AMG:    CycleParameter      : "<<CycleParameter      <<std::endl
    << " AMG:    NumBadCoarsenings   : "<<numBadCoarsenings   <<std::endl
    << " AMG:    CoarseningAlpha     : "<<alpha               <<std::endl
    << " AMG:    CoarseningGamma     : "<<gamma               <<std::endl
    << " AMG:    StrongDiagRatio     : "<<strongDiagRatio     <<std::endl
    << " AMG:    ForceFineRatio      : "<<forceFineRatio      <<std::endl
    << " AMG:    BadCoarseningRatio  : "<<badCoarseningRatio  <<std::endl
    << " AMG:    logging             : "<<(logging?"yes":"no")<<std::endl
    << " AMG:    DirectSolver        : "
    << Enum2String<SolverType>(directSolver)                  <<std::endl
    << " AMG:    InterpolationType   : "
    << Enum2String<AMGInterpolationType>(InterpolationType)   <<std::endl
    << " AMG:    SmootherType        : "
    << Enum2String<AMGSmootherType>(SmootherType)             <<std::endl
    << " AMG:    SmootherOmega       : "<<smootherOmega
    << std::endl;

    return out;
}

/**********************************************************/
} // namespace CoupledField

/**********************************************************
 * print-out operators
 **********************************************************/

namespace std {
template <typename T> std::ostream& operator <<
( std::ostream& out, const OLAS::HierarchyLevel<T>& level ) {
    return level.Print( out ); }
}

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
