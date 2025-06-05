// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAG_BASE_PDE_HH
#define FILE_MAG_BASE_PDE_HH

/*
 * MagBasePDE
 * > Idea taken from NACS
 * > Collect all functions that are common to MagneticPDE and MagEdgePDE (eventually also interesting for MagEdgeMxiedAVPDE)
 *
 * > Currently only used for hysteresis related functions (other parts commented out)
 */

#include "SinglePDE.hh"

namespace CoupledField
{

  class MagBasePDE: public SinglePDE {

  public:

    //! Constructor
    MagBasePDE( Grid *aGrid, PtrParamNode paramNode,
      PtrParamNode infoNode,
      shared_ptr<SimState> simState, Domain* domain );

    //!  Destructor
    virtual ~MagBasePDE();

    /*
     * Common functions for NODAL and EDGE formulation
     */
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    virtual void InitNonLin();

    //! Init the time stepping
    void InitTimeStepping();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();

    void DefineMagFluxDensity();

    // =======================================================================
    //   COILS
    // =======================================================================
    //@{ \name Attributes related to coils
    //! Map CoilID to coil definition
    std::map<Coil::IdType, shared_ptr<Coil> > coils_;

    //! Map regionIds to coil definitions
    typedef std::map<RegionIdType, shared_ptr<Coil> > CoilRegionMap;
    CoilRegionMap coilRegions_;

    //! Coefficients holding the current density for each coil
    std::map<RegionIdType, PtrCoefFct> coilCurrentDens_;

    //! Storage for CoefFunctions of external current density as source
    std::map<shared_ptr<Coil::Part>, PtrCoefFct> coilPartsExtJ_;

    //! Tells if there are coils excited by voltage
    bool hasVoltCoils_;
    //@}

    /*
     * Functions related to coils
     * Mostly copied from NACS (thanks and credits to Simetris!)
     */
    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
    //! specifying coil results.
    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

    //! Query parameter object for information on coils
    void ReadCoils();

    //! define special integrators for coils
    //! addition compared to NACS-version: Scaling factor
    //! (needed for magstrict coupling where we have to multiply with -1.0)
    void DefineCoilIntegrators(Double scaling);

    //! This method returns the RHS-integrator for a given current density,
    //! which is basically a BU-Integrator with a given identity-operator
    //! depending on the formulation.
    virtual LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "") = 0;

    // =======================================================================
    //   HYSTERESIS and Fixed Magnetization
    // =======================================================================
    /*
     * Objects and variables related to or required for hysteresis
     */
    //! CoefFunctions for postprocessing results
    shared_ptr<CoefFunctionMulti> polarization_;
    shared_ptr<CoefFunctionMulti> magnetization_;
    shared_ptr<CoefFunctionMulti> fieldIntensity_;
    
    //! same as bRHSRegions_ but for prescribed magnetization (via mat file)
    //! required to subtract fixed magnetization for computation of H in non-hysteretic regions
    std::map<RegionIdType,PtrCoefFct> mRHSRegions_;
    bool magnetizationSet_;

    //! In case of hysteresis, we set the region approximation in InitHystCoefs which has to be called
    //! prior to DefineIntegrators (where the regionApproximation usually is set)
    bool regionApproxSet_;
    std::map<RegionIdType,shared_ptr<ElemList> > SDLists_;

    /*
     * Functions related to hysteresis and fixed magnetization
     */
    //! Check all regions for hysteretic material behavior;
    //!for each found hysteretic region create a CoefFunctionHyst
    //!and add it to the map hysteresisCoefs_
    void InitMagnetization();

    PtrCoefFct GenerateHystereticCoefFunctions(RegionIdType actRegion){

      PtrCoefFct hystPol = hysteresisCoefs_->GetRegionCoef(actRegion);
      // shared_ptr<CoefFunction >
      PtrCoefFct curCoef = hystPol->GenerateMatCoefFnc("Reluctivity");

//      PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("MagPolarization");
//      polarization_->AddRegion( actRegion, hystOutput);
//
//      PtrCoefFct hystOutput2 = hystPol->GenerateOutputCoefFnc("MagMagnetization");
//      magnetization_->AddRegion( actRegion, hystOutput2);
//
//      PtrCoefFct hystOutput3 = hystPol->GenerateOutputCoefFnc("MagFieldIntensityHyst");
//      fieldIntensity_->AddRegion( actRegion, hystOutput3);

      return curCoef;
    }

    virtual LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation ) = 0;

    virtual BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ) = 0;

    // =======================================================================
    //   MISC
    // =======================================================================
    /*
     * Common variables for NODAL and EDGE formulation
     */
    //! Formulation type of magnetic PDE
    typedef enum { BASE, NODAL, EDGE } MagFormulationType;

    //!true, if 3d computation
    bool is3d_;

    //! formulation
    MagFormulationType formulation_;

    bool fluxDensityDefined_;

    //! In case of a transient / harmonic 3D simulation using NODAL formulation, we need
    //! an additional scalar potential to ensure HCurl compatibility (A-V-formulation)
    bool isMixed_;

    //! in mixed case, we have to check if there are any regions with conductivity at all
    bool anyRegionHasConductivity_;

    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> relucTensor_;

    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> reluc_;

    //! Coefficient function, containing the overall conductivity
    shared_ptr<CoefFunctionMulti> conduc_;

    bool coilOptimization_;

    inline bool DoCoilOptimization() { return coilOptimization_; }

    void GenerateLorentzForceResults(CoupledField::StdVector<std::string> &vecComponents, boost::shared_ptr<CoupledField::CoefFunctionMulti> &tcdCoef,
      CoupledField::PtrCoefFct &bFunc, CoupledField::Global::ComplexPart &part, boost::shared_ptr<CoupledField::BaseFeFunction> &feFct) ;

    template<typename FE>
    void GenerateVwpForcesResults(CoupledField::StdVector<std::string> &vecComponents, CoupledField::PtrCoefFct &bFunc, 
      boost::shared_ptr<CoupledField::BaseFeFunction> &feFct) ;

    //! This coefficient function describes the velocity field.
    shared_ptr<CoefFunctionMulti> VelocityCoef_;
    
    //! Function to get the velocity field of each region defined in the xml-Script for all childs of the base PDE
    void ReadRegionVelocityField(std::string velocityId, shared_ptr<ElemList> actSDList, RegionIdType actRegion, bool& coefUpdateGeo);

  };
} // end of namespace
#endif
