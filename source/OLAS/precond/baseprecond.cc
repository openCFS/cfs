// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/SingleVector.hh"
#include "MatVec/basematrix.hh"
#include "MatVec/basevector.hh"
#include "MatVec/sbmmatrix.hh"
#include "MatVec/sbmvector.hh"
#include "MatVec/stdmatrix.hh"
#include "baseprecond.hh"

namespace CoupledField {

  static EnumTuple precondTypeTuples[] = 
  {
    EnumTuple( BasePrecond::NOPRECOND, "noPrecond" ),
    EnumTuple( BasePrecond::ID, "Id" ),
    EnumTuple( BasePrecond::MG, "MG"),
    EnumTuple( BasePrecond::JACOBI, "Jacobi"),
    EnumTuple( BasePrecond::SSOR, "SSOR" ),
    EnumTuple( BasePrecond::ILU0, "ILU0" ),
    EnumTuple( BasePrecond::ILUTP, "ILUTP"),
    EnumTuple( BasePrecond::ILUK, "ILUK"),
    EnumTuple( BasePrecond::ILDL0, "ILDL0" ),
    EnumTuple( BasePrecond::ILDLK, "ILDLK" ),
    EnumTuple( BasePrecond::ILDLTP, "ILDLTP"),
    EnumTuple( BasePrecond::ILDLCN, "ILDLCN"),
    EnumTuple( BasePrecond::IC0, "IC0" ),
  };

  Enum<BasePrecond::PrecondType> BasePrecond::precondType = \
  Enum<BasePrecond::PrecondType>("Preconditioner Types",
      sizeof(precondTypeTuples) / sizeof(EnumTuple),
      precondTypeTuples); 


  void BaseStdPrecond::Apply( const BaseMatrix &sysmat, const BaseVector &r, 
                           BaseVector &z ) const {
    const StdMatrix& stdsysmat = dynamic_cast<const StdMatrix&>(sysmat);
    const SingleVector& stdr = dynamic_cast<const SingleVector&>(r);
    SingleVector& stdz = dynamic_cast<SingleVector&>(z);

    Apply(stdsysmat,stdr,stdz);
  }
  
  void BaseStdPrecond::Setup( BaseMatrix &sysMat ) {
    Setup( dynamic_cast<StdMatrix&>(sysMat) );
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
