// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_IDSTDPRECOND_HH
#define OLAS_IDSTDPRECOND_HH

#include "MatVec/sbmvector.hh"
#include "baseprecond.hh"

namespace CoupledField {

  //! Identity Preconditioner for Standard Matrices

  //! This class implements an identity preconditioner for problems involving
  //! standard matrices. Since preconditioning with the identity matrix does
  //! not really change anything this is a quite simple task. The Setup routine
  //! does not really do anything and the Apply routine simply returns the rhs
  //! vector unchanged as sol vector. The class is provided in order to render
  //! constant case distinctions unneccessary.
class SBM_Matrix;
class SingleVector;
class StdMatrix;

  class IdPrecondStd : public BaseStdPrecond {

  public:

    using BasePrecond::Apply;
    using BasePrecond::Setup;

    //! Default Constructor
    IdPrecondStd() {
    };

    //! Default Destructor
    ~IdPrecondStd() {
    };

    //! Application of preconditioner

    //! This method applies the identity preconditioner to solve Id * sol = rhs
    //! for sol. As a consequence it simply copies the rhs entries into the sol
    //! vector.
    void Apply( const StdMatrix &sysmat, const SingleVector &rhs,
                SingleVector &sol ) const;

    //! Triggers setup of the identity preconditioner

    //! Since the application of the identity precondtioner requires nothing
    //! but to copy some vector entries, the setup method does not do anything.
    void Setup( StdMatrix &sysmat ) {
    }

    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is ID.
    PrecondType GetPrecondType() const {
      return ID;
    }

  };
  
  //! Identity Preconditioner for SBM Matrices

  //! This class implements an identity preconditioner for problems involving
  //! SBM matrices. Since preconditioning with the identity matrix does
  //! not really change anything this is a quite simple task. The Setup routine
  //! does not really do anything and the Apply routine simply returns the rhs
  //! vector unchanged as sol vector. The class is provided in order to render
  //! constant case distinctions unneccessary.
  class IdPrecondSBM : public BaseSBMPrecond {

  public:

    using BasePrecond::Apply;
    using BasePrecond::Setup;

    //! Default Constructor
    IdPrecondSBM() {
    };

    //! Default Destructor
    ~IdPrecondSBM() {
    };

    //! Application of preconditioner

    //! This method applies the identity preconditioner to solve Id * sol = rhs
    //! for sol. As a consequence it simply copies the rhs entries into the sol
    //! vector.
    void Apply( const SBM_Matrix &sysmat, const SBM_Vector &rhs,
                SBM_Vector &sol ) const;

    //! Triggers setup of the identity preconditioner

    //! Since the application of the identity precondtioner requires nothing
    //! but to copy some vector entries, the setup method does not do anything.
    void Setup( SBM_Matrix &sysmat ) {
    }

    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is ID.
    PrecondType GetPrecondType() const {
      return ID;
    }

  };

}

#endif
