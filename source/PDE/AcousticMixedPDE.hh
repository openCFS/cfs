// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AcousticMixedPDE.hh
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

#ifndef ACOUSTICMIXEDPDE_HH_
#define ACOUSTICMIXEDPDE_HH_

#include "SinglePDE.hh"

namespace CoupledField{

class AcousticMixedPDE : public SinglePDE{

  public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    AcousticMixedPDE( Grid* aGrid, PtrParamNode paramNode );

    virtual ~AcousticMixedPDE(){};

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> >
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Initialize all the nodes by this value
    void SetInitialCondition();

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! Initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling) {
      EXCEPTION("Coupling not supported for acousticPDE.");
    };

    //! Calculate coupling terms
    void CalcOutputCoupling() {
      EXCEPTION("Coupling not supported for acousticPDE.");
    };

    //! Returns if PDE can compute the quantity
    bool HasOutput(SolutionType output) { return false;};

    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! perform postprocessing on data
    void CalcResults( shared_ptr<BaseResult> result );

  protected:

    //!  Define available postprocessing results
    void DefinePrimaryResults();

    //! Init the time stepping
    void InitTimeStepping();
  private:

    bool usePiola_;
};

}


#endif /* ACOUSTICMIXEDPDE_HH_ */
