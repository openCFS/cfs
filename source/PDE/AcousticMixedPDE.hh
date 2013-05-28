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

//! forward class declaration
class CoefFunctionMulti;

class AcousticMixedPDE : public SinglePDE{

  public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    AcousticMixedPDE( Grid* aGrid, PtrParamNode paramNode, 
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain );

    virtual ~AcousticMixedPDE(){};

  protected:

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> >
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    virtual void DefineNcIntegrators() {
      EXCEPTION("ncInterfaces are not implemented for AcousticMixedPDE");
    };

    //! define all (linearform) integrators needed for this pde
    void DefineRhsLoadIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Initialize all the nodes by this value
    void SetInitialCondition();

    //! define all (bilinearform) integrators needed for this pdewith template
    //! for the space dimension
    template<class DATA_TYPE, UInt DIM>
    void DefineIntegratorsTempl();

    template<class DATA_TYPE, UInt DIM>
    void DefineRhsLoadIntegratorsTempl();

    //!  Define available postprocessing results
    void DefinePrimaryResults();

    //! Init the time stepping
    void InitTimeStepping();

    //! create feFunction for meanFluidMech velocity
    void CreateMeanFlowFunction(StdVector<std::string> dofNames);

  private:

    bool usePiola_;

    //! Coefficient function for the flow field

    //! This coefficient function describes the flow field. As this 
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> meanFlowCoef_;
        

    bool penalized_;

    bool doFluxTerm_;
};

}


#endif /* ACOUSTICMIXEDPDE_HH_ */
