// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/stdmatrix.hh"
#include "MatVec/sbmmatrix.hh"

#include "MatVec/basevector.hh"
#include "MatVec/sbmvector.hh"

#include "baseprecond.hh"

namespace CoupledField {

  void BaseStdPrecond::Apply( const BaseMatrix &sysmat, const BaseVector &r, 
                           BaseVector &z ) const {
    TRY_CAST {
      CONSTREFCAST(sysmat,StdMatrix,stdsysmat);
      CONSTREFCAST(r,SingleVector,stdr);
      REFCAST(z,SingleVector,stdz);
      Apply(stdsysmat,stdr,stdz);
    } CATCH_CAST;
  }
  
  void BaseStdPrecond::Setup( BaseMatrix &sysMat ) {
    TRY_CAST {
      REFCAST( sysMat, StdMatrix, stdMat );
      Setup( stdMat );
    } CATCH_CAST;
  }

  void BaseSBMPrecond::Apply( const BaseMatrix& sysmat, const BaseVector& r, 
                              BaseVector& z ) const {
      const SBM_Matrix& sbmsysmat = dynamic_cast<const SBM_Matrix&>(sysmat);
      const SBM_Vector& sbmr      = dynamic_cast<const SBM_Vector&>(r);
      SBM_Vector& sbmz            = dynamic_cast<SBM_Vector&>(z);
      Apply(sbmsysmat,sbmr,sbmz);
    }

  void BaseSBMPrecond::Setup( BaseMatrix &A ) {
    SBM_Matrix& sbmA = dynamic_cast<SBM_Matrix&>(A);
    Setup(sbmA);
  }

}
