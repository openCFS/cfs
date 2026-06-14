#include "MatVec/BaseMatrix.hh"
#include "MatVec/StdMatrix.hh"
#include "DataInOut/ProgramOptions.hh"

#include "def_use_pardiso.hh"
#ifdef USE_PARDISO
#include "OLAS/external/pardiso/PardisoSolver.hh"
#endif

#include "OLAS/precond/LDLPrecond.hh"

namespace CoupledField {

  // Constructor
  template<typename T>
  LDLPrecond<T>::LDLPrecond( PtrParamNode precondNode, PtrParamNode olasInfo ) {
    this->xml_      = precondNode;
    this->infoNode_ = olasInfo->Get("ldl", progOpts->DoDetailedInfo() ? ParamNode::APPEND : ParamNode::INSERT);

    // adaptive refresh threshold
    precondNode->GetValue("maxApplyIter", maxApplyIter_, ParamNode::INSERT);

#ifdef USE_PARDISO
    // PARDISO reads its own params (ordering, posDef, ...) from the same <ldl>
    // node; absent values fall back to PARDISO defaults. Keep symStruct off so
    // complex-symmetric systems get mType 6, and IterRefineSteps stays 0 (no
    // refinement against the mismatched frozen matrix).
    pardiso_ = shared_ptr<PardisoSolver<T>>( new PardisoSolver<T>(precondNode, olasInfo) );
#else
    EXCEPTION("LDLPrecond requires building openCFS with USE_PARDISO.");
#endif

    this->infoNode_->Get("maxApplyIter")->SetValue(maxApplyIter_);
  }

  // Destructor: shared_ptr releases the PardisoSolver (and its factors) here
  template<typename T>
  LDLPrecond<T>::~LDLPrecond() {}

  // Setup: freeze/refresh decision, called once per wavenumber
  template<typename T>
  void LDLPrecond<T>::Setup( BaseMatrix& sysmat ) {
    const bool refresh = !factored_ || RefreshNeeded();

    if ( refresh ) {
      // First call: symbolic + numeric. Later: numeric only, since the
      // sparsity pattern is identical across wavenumbers.
      pardiso_->Setup( sysmat );
      factored_ = true;
    }
    // else: keep frozen factors, skip factorization

    lastApplyCount_ = 0;        // reset per-solve iteration proxy
    this->readyToUse_ = true;
  }

  // Apply: z = M^-1 r via PARDISO phase 33 with the frozen factors
  template<typename T>
  void LDLPrecond<T>::Apply( const BaseMatrix& sysmat, const BaseVector& r,
                             BaseVector& z ) {
    pardiso_->Solve( sysmat, r, z );
    ++lastApplyCount_;
  }

  // RefreshNeeded: refactorize if the previous COCR solve became too slow
  template<typename T>
  bool LDLPrecond<T>::RefreshNeeded() const {
    return ( maxApplyIter_ > 0 && lastApplyCount_ > maxApplyIter_ );
  }

  // Explicit template instantiation
  template class LDLPrecond<Double>;
  template class LDLPrecond<Complex>;

}