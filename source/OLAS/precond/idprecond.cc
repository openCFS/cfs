// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/SingleVector.hh"
#include "MatVec/sbmvector.hh"
#include "idprecond.hh"

namespace CoupledField {
class SBM_Matrix;
class StdMatrix;

  void IdPrecondStd::Apply( const StdMatrix &sysmat, const SingleVector &rhs,
                         SingleVector &sol ) const {
    sol = rhs;
  }
  
  void IdPrecondSBM::Apply( const SBM_Matrix &sysmat, const SBM_Vector &rhs,
                            SBM_Vector &sol ) const {
    sol = rhs;
  }  
}
