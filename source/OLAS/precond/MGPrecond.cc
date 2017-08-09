//#define  INCLUDE_MULTIGRID_CC_FILES
#include <def_expl_templ_inst.hh>

#include "OLAS/precond/MGPrecond.hh"
#include "DataInOut/ProgramOptions.hh"

namespace CoupledField {
/**********************************************************/
/*
template <typename T>
MGPrecond<T>::MGPrecond( PtrParamNode params )
    : params_( params ),
      report_( NULL ),
      AMG_( NULL )
{
  // The MG preconditioner is currently not working
  WARN("The MG-Preconditioner is not yet ported from 1 to 0 based"
      "numbering, i.e. it will NOT work correctly!" );
}    
*/
/**********************************************************/

template <typename T>
MGPrecond<T>::MGPrecond( const StdMatrix   &matrix,
                        PtrParamNode precondNode,
                        PtrParamNode olasInfo )
    : params_( precondNode ),
      report_( olasInfo ),
      AMG_( NULL )
{

  // Set pointers to communication objects
  this->xml_ = precondNode;
  this->infoNode_ = olasInfo->Get("MG", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);

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

  // The new implementation only works with the auxiliary-
  // matrix approach
  EXCEPTION("MGPrecond::Setup : wrong constructor was called!");

}

/**********************************************************/

template <typename T>
void MGPrecond<T>::SetupMG(StdMatrix& sysmatrix,
                           StdMatrix& auxmatrix,
                           const AMGType amgType,
                           const StdVector< StdVector< Integer> >& edgeIndNode,
                           const StdVector<Integer>& nodeNumIndex)
{
    // cast the standard matrix objects (AMG works with
    // CRS matrices only)
    CRS_Matrix<T>& pSysMatrix = dynamic_cast<CRS_Matrix<T>&>(sysmatrix);
    CRS_Matrix<T>& pAuxMatrix = dynamic_cast<CRS_Matrix<T>&>(auxmatrix);

    AMG_ = new AMGSolver<T>(params_, report_);

    // setup the serial matrix object
    if( !AMG_->Setup(&pSysMatrix, &pAuxMatrix, amgType, edgeIndNode, nodeNumIndex) ) {
      EXCEPTION( "could not set up the AMG preconditioner");
    }

AMG_->Print(std::cout);
    this->readyToUse_ = true;
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
      EXCEPTION( "AMG preconditioner used in undefined "
          "state. It is a bug that you could reach this line.");
    }
  } else {
    EXCEPTION( "AMG preconditioner used before setup");
  }


    }

/**********************************************************/

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class MGPrecond<Double>;
//template class MGPrecond<Complex>;
#endif

} // namespace CoupledField
