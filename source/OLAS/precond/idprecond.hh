// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_IDSTDPRECOND_HH
#define OLAS_IDSTDPRECOND_HH

#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "precond/baseprecond.hh"

namespace OLAS {

  //! Identity Preconditioner for Standard Matrices

  //! This class implements an identity preconditioner for problems involving
  //! standard matrices. Since preconditioning with the identity matrix does
  //! not really change anything this is a quite simple task. The Setup routine
  //! does not really do anything and the Apply routine simply returns the rhs
  //! vector unchanged as sol vector. The class is provided in order to render
  //! constant case distinctions unneccessary.
  class IdPrecond : public BaseStdPrecond {

  public:

    //! Default Constructor
    IdPrecond() {
    };

    //! Default Destructor
    ~IdPrecond() {
    };

    //! Application of preconditioner

    //! This method applies the identity preconditioner to solve Id * sol = rhs
    //! for sol. As a consequence it simply copies the rhs entries into the sol
    //! vector.
    void Apply( const StdMatrix &sysmat, const SparseVector &rhs,
		SparseVector &sol ) const {
      sol = rhs;
    }

    //! Triggers setup of the identity preconditioner

    //! Since the application of the identity precondtioner requires nothing
    //! but to copy some vector entries, the setup method does not do anything.
    /// diagonal entries of the system matrix
    void Setup( StdMatrix &sysmat ) {
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
