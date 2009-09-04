// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#define  INCLUDE_MULTIGRID_CC_FILES
// #include "multigrid/multigrid.hh"
#include "precond/mgprecond.hh"

namespace CoupledField {
/**********************************************************/

template <typename T>
MGPrecond<T>::MGPrecond( OLAS_Params *params )
    : params_( params ),
      report_( NULL ),
      AMG_( NULL )
{
}    

/**********************************************************/

template <typename T>
MGPrecond<T>::MGPrecond( const StdMatrix   &matrix,
                         OLAS_Params *params,
                         OLAS_Report *report )
    : params_( params ),
      report_( report ),
      AMG_( NULL )
{
}

/**********************************************************/

template <typename T>
MGPrecond<T>::~MGPrecond()
{

    delete AMG_;
}

/**********************************************************/

template <typename T>
void MGPrecond<T>::Setup( StdMatrix& sysmatrix )
{

    if( ! AMG_ ) {
      if( NULL == (AMG_ = New AMGSolver<T>(params_)) ) {
	Error( "MGPrecond::Setup: could not create serial "
	       "AMG object", __FILE__, __LINE__ );
      } else {
	AMG_->Reset();
      }
    }
    // cast the matrix to a serial CRS matrix and call Setup

    CRS_Matrix<T>& crsSysMatrix = dynamic_cast<CRS_Matrix<T>&>(sysmatrix);
      // eventually change the layout of the matrix
      if( CRS_Matrix<T>::LEX_DIAG_FIRST != crsSysMatrix.GetCurrentLayout() ) {
	(*cla)
	  << " --------------------------------------\n"
	  << " MGPrecond::Setup( StdMatrix& ):\n"
	  "   changing layout of  passed serial system matrix:"
	  << std::endl;
	crsSysMatrix.ChangeLayout( CRS_Matrix<T>::LEX_DIAG_FIRST );
	(*cla) << " --------------------------------------"<<std::endl;
      }
    // setup the serial AMG object
    if( ! AMG_->Setup(&crsSysMatrix) ) {
      Error( "MGPrecond::Setup: could not set up the AMG "
	     "preconditioner\n", __FILE__, __LINE__ );
    }
    this->readyToUse_ = true;  
}

/**********************************************************/

template <typename T>
void MGPrecond<T>::Setup( const StdMatrix& sysmatrix,
                          const StdMatrix& auxmatrix )
{

    ///////////////////////////////////////////////////////////
    Error( "MGPrecond::Setup(StdMatrix,StdMatrix) not yet "
           "implemented\n", __FILE__, __LINE__ );
    ///////////////////////////////////////////////////////////

    // cast the standart matrix objects (AMG works with
    // CRS matrices only)
    CRS_Matrix<T>* const pSysMatrix = dynamic_cast<CRS_Matrix<T>*>(&sysmatrix);
    CRS_Matrix<T>* const pAuxMatrix = dynamic_cast<CRS_Matrix<T>*>(&auxmatrix);

    if( ! pSysMatrix && ! pAuxMatrix ) {
        Error( "MGPrecond::Setup: AMG works with CRS matrices "
               "only\n", __FILE__, __LINE__ );
    }

    // setup the serial matrix object
    if( !AMG_.Setup(&sysmatrix, &auxmatrix) ) {
        Error( "could not set up the AMG preconditioner\n",
               __FILE__, __LINE__ );
    }
}

/**********************************************************/

template <typename T>
void MGPrecond<T>::Apply( const StdMatrix& sysmatrix,
    const SingleVector& rhs,
    SingleVector& sol ) const
    {

  if( this->readyToUse_ ) {
    sol.Init();
    // We have set up a serial AMG. This might also have
    // been the decition due to a single process run of
    // a parallel program. Also in this case the casts to
    // (const) Vector<T> are OK, because Vector<T> is a
    // base class of ParVector<T>
    if( AMG_ ) {
      const Vector<T>& vectorRhs = dynamic_cast<const Vector<T>&>(rhs);
      Vector<T>& vectorSol = dynamic_cast<Vector<T>&>(sol);

      // AMG cycle
      AMG_->Cycle( vectorRhs, vectorSol );
    } else {
      Error( "AMG preconditioner used in undefined "
          "state. It is a bug that you could reach "
          "this line.", __FILE__, __LINE__ );
    }
  } else {
    Error( "AMG preconditioner used before setup",
        __FILE__, __LINE__ );
  }
    }

/**********************************************************/
} // namespace CoupledField
