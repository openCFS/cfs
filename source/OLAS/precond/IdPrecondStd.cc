#include "MatVec/StdMatrix.hh"
#include "MatVec/SBM_Matrix.hh"

#include "MatVec/SBM_Vector.hh"
#include "MatVec/SingleVector.hh"

#include "IdPrecondStd.hh"

namespace CoupledField {

  IdPrecondStd::IdPrecondStd(PtrParamNode xml, PtrParamNode olasInfo ) {
    this->infoNode_ = olasInfo->Get("idPrecond");
  }
  
  void IdPrecondStd::Apply( const BaseMatrix &sysmat, const BaseVector &rhs,
                         BaseVector &sol ) {
    sol = rhs;
  }
  
//  void IdPrecondSBM::Apply( const SBM_Matrix &sysmat, const SBM_Vector &rhs,
//                            SBM_Vector &sol ) const {
//    sol = rhs;
//  }  
}
