// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PIEZOCOUPLING_HH
#define FILE_PIEZOCOUPLING_HH

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField {
class BaseResult;
class EntityIterator;
template <class TYPE> class Matrix;
template <class TYPE> class Vector;

// Forward declarations
class BaseMaterial;
class SinglePDE;
class BiLinearForm;

  //! Implements the definition of pairwise piezo-coupling
  
  //! This class implements the piezoelectric coupling based on the 
  //! volumetric material law.
  //! Within this object, pde1_ refers to the mechanical PDE,
  //! whereas  pde2_ refers to the electric PDE.
  class PiezoCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
    PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode,
                   PtrParamNode infoNode,
                   shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~PiezoCoupling();

  protected:
    // In this special case we need to reset the time stepping schemes
    // defined by the single PDEs to ensure consistency (This is the case for almost every coupledPDE)
    // therefore it is also crucial, that this method is called after(!) the single
    // PDEs defined their schemes. We check for that and issue an error if its true.
    // In future, it would be more appropriate to call the InitTimeStepping from the driver rather than from the
    // SinglePDEs such there is exaclty one place to initialize timestepping and a clean overload can be ensured...
    //! specify the time stepping
    virtual void InitTimeStepping();

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    // for hysteresis we need rhs loads that contain information about both pdes
    void DefineRhsLoadIntegrators();
    
    void DefinePBCIntegrators(shared_ptr<BaseFeFunction>& fe1, shared_ptr<BaseFeFunction>& fe2);

    template<typename DATA_TYPE>
    BiLinearForm* GetNormalPiezoFluxIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec, Double factor,
                                               BiLinearForm::CouplingDirection cplDir, bool fluxOpA);

    template<typename DATA_TYPE>
    BiLinearForm* GetNormalPiezoStrainIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec, Double factor,
                                                 BiLinearForm::CouplingDirection cplDir, bool fluxOpA);

    //! Define available results
    void DefineAvailResults();
    
    //! \copydoc BasePairCoupling::DefinePostProcResults
    void DefinePostProcResults();
    
    //! Returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
    BaseBDBInt * GetStiffIntegrator( BaseMaterial* actSDMat,
                                     RegionIdType regionId,
                                     bool isComplex );

    //! Returns a stiffness integrator appropriate to the actual problem (e.g. 3D) with the material tensor scaled by a given factor
    BaseBDBInt* GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex, PtrCoefFct scalingFactor);

    // returns true if hysteretic behavior was found and false otherwise
    // case true:
    //  elecToMechInt = integrator that couples electric field to mechanics
    //  mechToElecInt = integrator that couples mechancis to electric field
    // case false:
    //  elecToMechInt = output of GetStiffIntegrator
    //  mechToElecInt = NULL and will not be used; instead isCounterpart will be set for assemblys
    bool GetStiffIntegratorHyst( BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex, BaseBDBInt** elecToMechInt, BaseBDBInt** mechToElecInt);

    //! Subtype of related mechanical PDE
    std::string subType_;
    
  private:

  };


} // end of namespace

#endif
