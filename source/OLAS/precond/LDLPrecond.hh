#ifndef OLAS_LDLPRECOND_HH
#define OLAS_LDLPRECOND_HH

#include "BasePrecond.hh"

namespace CoupledField {

  // Forward declaration; the concrete type is only needed in the .cc
  template<typename> class PardisoSolver;

  //! Frozen-factorization preconditioner backed by PARDISO.
  //!
  //! For the 2.5D wavenumber sweep: PARDISO factorizes A(k0) once and reuses
  //! the L D L^T factors as a preconditioner (phase-33 solve) for neighbouring
  //! wavenumbers solved with COCR. Refresh is adaptive: if the previous COCR
  //! solve needed more than 'maxApplyIter' applications (~ iterations), the
  //! factors are recomputed at the next Setup. M^-1 r is delegated to PARDISO,
  //! so the fill-reducing permutation and scaling are handled internally.
  template<typename T>
  class LDLPrecond : public BasePrecond {
  public:
    LDLPrecond( PtrParamNode precondNode, PtrParamNode olasInfo );
    ~LDLPrecond() override;

    //! Per wavenumber: (re)factorize or keep the frozen factors
    void Setup( BaseMatrix& sysmat ) override;

    //! z = M^-1 r via a frozen PARDISO solve
    void Apply( const BaseMatrix& sysmat, const BaseVector& r,
                BaseVector& z ) override;

    BasePrecond::PrecondType GetPrecondType() const override {
      return BasePrecond::LDL;
    }

  private:
    //! Owns the frozen PARDISO factorization
    shared_ptr<PardisoSolver<T>> pardiso_;

    //! A frozen factorization exists
    bool factored_ = false;

    //! Apply() calls since last Setup ~ last solve's iteration count
    UInt lastApplyCount_ = 0;

    //! Adaptive refresh threshold (XML, default 50)
    UInt maxApplyIter_ = 50;

    //! true -> recompute factors at this Setup
    bool RefreshNeeded() const;

    LDLPrecond() = delete;  // needs XML + info nodes
  };

}

#endif // OLAS_LDLPRECOND_HH