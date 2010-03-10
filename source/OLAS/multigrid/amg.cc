// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#include "multigrid/amg.hh"

#if defined(__INTEL_COMPILER)
#include "algsys/olascomm_impl.hh"
#endif
/**********************************************************/
#ifdef AMG_TRACE_CF_SPLITTING
#define  AMG_SETTING_SET_KEEP_TOPOLOGY \
         settings.keepTopology = true;
#else
#define  AMG_SETTING_SET_KEEP_TOPOLOGY
#endif
/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_AMG
#define DEBUG_AMG
#endif // DEBUG_AMG
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR

#ifndef  LOGLINE
#define  LOGLINE  " --------------------------------------\n"
#endif
/**********************************************************/

namespace CoupledField {
/**********************************************************/

template <typename T>
AMGSolver<T>::AMGSolver( OLAS_Params* params )
            : Parameters_( params ),
              Hierarchy_( NULL ),
              prepared_( false ),
              logging_( params == NULL ? false :
                        params->GetBoolValue("AMG_logging") )
{
}

/**********************************************************/

template <typename T>
AMGSolver<T>::~AMGSolver()
{

    Reset();
}

/**********************************************************/

template <typename T>
bool AMGSolver<T>::Setup( const CRS_Matrix<T>* const sys_matrix,
                          const CRS_Matrix<T>* const aux_matrix  )
{

    if( aux_matrix ) {
        WARN( "AMGSolver::Setup: please note that the auxiliary"
                 " matrix is not yet used\n", __FILE__, __LINE__ );
    }

    // reset the whole object
    Reset();
    Hierarchy_ = New HierarchyLevel<T>();

    // set logging state
    Hierarchy_->SetLogging( logging_ );
    // insert the matrices into the top level
    Hierarchy_->InsertSysMatrix( sys_matrix );
    Hierarchy_->InsertAuxMatrix( aux_matrix );

    // create a settings object for the hierarchy
    typename HierarchyLevel<T>::Settings settings( Parameters_ );
    // internal, eventual developing settings
    AMG_SETTING_SET_KEEP_TOPOLOGY

    // setup
    if( false == Hierarchy_->Setup(&settings) ) {
        Reset();
        return false;
    }

    return  prepared_ = true;
}

/**********************************************************/

template <typename T>
void AMGSolver<T>::Reset()
{

    // Remove system and auxiliary matrix from top level
    // of hierarchy, since they are not intended to be owned
    // by this object.
    if( Hierarchy_ ) {
        Hierarchy_->UnhookSysMatrix();
        Hierarchy_->UnhookAuxMatrix();
    }
    
    delete Hierarchy_;  Hierarchy_ = NULL;

    prepared_ = false;
}

/**********************************************************/

template <typename T>
void AMGSolver<T>::Cycle( const Vector<T>& rhs,
                                Vector<T>& sol ) const
{

    if( !prepared_ ) {
        WARN( "called AMGSolver::Cycle at a non prepared AMG" );
        return;
    }

    typename HierarchyLevel<T>::Settings settings( Parameters_ );

    if( logging_ ) (*cla)<<LOGLINE<<" AMG: starting cycle\n";
    // We need this dirty cast to a non-constant HierarchyLevel,
    // since the cycle will change the HierarchyLevel's temporary
    // vectors.
    static_cast<HierarchyLevel<T>*>(Hierarchy_)
                ->Cycle( rhs, sol,
                         settings.numPreSmoothings,
                         settings.numPostSmoothings,
                         settings.CycleParameter );

    if( logging_ )  (*cla) <<" AMG: finished cycle\n";
    if( logging_ )  (*cla) << LOGLINE;
}

/**********************************************************/

template <typename T>
std::ostream& AMGSolver<T>::Print( std::ostream& out ) const
{
    out << LOGLINE" AMG solver:"<< std::endl;

    if( prepared_ ) {
        Hierarchy_->Print( out );
    } else {
        out << " AMG is not prepared" << std::endl;
    }

    out << LOGLINE;
    return out;
}

/**********************************************************/
} // namespace CoupledField

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
