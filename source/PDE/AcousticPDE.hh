#ifndef ACOUSTICPDE_HH
#define ACOUSTICPDE_HH

#include "SinglePDE.hh"

namespace CoupledField{

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;
  class CoefFunctionMulti;

  class AcousticPDE : public SinglePDE{

  public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    AcousticPDE( Grid* aGrid, PtrParamNode paramNode,
                 PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain );

    virtual ~AcousticPDE(){};

    //! Indicate that acoustic PDE takes part in coupling to mechanics.
    void SetMechanicCoupling() {
      isMechCoupled_ = true;
    }

    //! Return acoustic formulation. Can either be pressure or potential.
    SolutionType GetFormulation() const { return formulation_; }

    /** @see virtual SinglePDE::GetNativeSolutionType()
      *  @return depending on the formulation: ACOU_POTENTIAL or ACOU_PRESSURE. */
     SolutionType GetNativeSolutionType() const { return formulation_; }

     //! check, if PDE has complex material parameters;
     bool IsMaterialComplex()  {
     	return complexFluidFormulation_;
     }

  protected:

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    // Struct to collect the coefficients of the rational function approximation in case of the TDEF formulation
    struct rationalCoeffsTDEF {
      StdVector<PtrCoefFct> AC, BC, CC, AlphaC, BetaC, GammaC;
      StdVector<PtrCoefFct> AV, BV, CV, AlphaV, BetaV, GammaV;
    };
    // Struct to collect the number of real and complex poles of the rational function approximation for the TDEF formulation
    struct NumPolesTDEF {
      Vector<unsigned int> realC;
      Vector<unsigned int> complC;
      Vector<unsigned int> realV;
      Vector<unsigned int> complV;
    };

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators();
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //!  Define available primary results
    void DefinePrimaryResults();

    //! Check if there is only one material 
    void CheckIfIsOnlyOneMaterial();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();
    
    //! Init the time stepping
    void InitTimeStepping();

    //! create feFunction for meanFluidMech velocity
    void CreateMeanFlowFunction(StdVector<std::string> dofNames);

    //! create transient PML integrators
    template<UInt DIM>
    void DefineTransientPMLInts(shared_ptr<ElemList> eList,std::string id,
    		                    RegionIdType actRegion, std::string tempId);

    //! create TDEF integrators (The time domain equivalent fluid formulation requires terms additional to the standard wave equation, see Diss Maurerlehner 2023)
    void DefineTDEFIntegrators(UInt iRegion, string polyId, string integId, shared_ptr<ElemList> actSDList, rationalCoeffsTDEF &TDEFcoeffs);

    //! create Surface integrators for ABC in case the BC is adjacent to a TDEF region (TDEF formulation applied)
    void DefineTDEFABCSurfaceIntegrators(PtrCoefFct &coeffStiffTDEF, RegionIdType &aRegion, std::set<RegionIdType> &adjVolRegion, shared_ptr<EntityList> &actSDList, rationalCoeffsTDEF &TDEFcoeffs);

    //! calculate the ABC coefficients for TDEF region (TDEF formulation applied)
    void CalcCoefsTDEFABC(UInt iRegion, PtrCoefFct &factor, PtrCoefFct &c0, PtrCoefFct &dens, PtrCoefFct &omegaTrg, PtrCoefFct &coeffDamp, PtrCoefFct &coeffStiffTDEF);

    //! read TDEF coefficients (coefficients of rational function approx. of the complex, frequency-dependent inverse equivalent fluid parameters) for the given region
    void ReadTDEFCoefficients(UInt iRegion, rationalCoeffsTDEF &TDEFcoeffs);

    //! evaluate the rational function approximation of the equivalent fluid parameters (inverse bulk modulus and inverse density) at given frequency for a given region
    void EvalTDEFRationalFncs(UInt iRegion, PtrCoefFct actFreq, rationalCoeffsTDEF &TDEFcoeffs);

    //    //! Set special RHS values
    //    virtual void SetRhsValues();

  private:
    //! Defines integrators for Nitsche coupling of an unknown on one specific
    //! interface.
    //! \param iface Interface for which the coupling is defined.
    template<UInt DIM, bool IS_COMPLEX>
    void DefineNitscheCoupling(NcInterfaceInfo &iface);

    // //! Compute a global factor for coupling to other PDEs. Currently only used for
    // //! the mechAcou coupling to get a symmetric matrix when in acouPotential formulation.
    // PtrCoefFct SetGlobalCouplingFactor();

    //! stores if the Acoustic PDE is in potential or pressure form
    SolutionType formulation_;

    //! if true, speed of sound squared ist at Laplace term!
    bool sosAtLaplace_;

    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;

    //! acoustic source density
    shared_ptr<CoefFunctionMulti> acousticSourceDensityCoef_;

    //! acoustic rhs source density
    shared_ptr<CoefFunctionMulti> acousticRhsDensityCoef_;

    //! Coefficient function for the flow field

    //! This coefficient function describes the flow field. As this 
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> meanFlowCoef_;
    
    //! This coefficient function describes the divergence of the flow field. As this
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> divMeanFlowCoef_;

    //! This coefficient function describes the temperatur field. As this
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> meanTemperatureCoef_;

    //! This coefficient function describes the lamb vector field. As this
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> lambCoef_;

    //! store convective bilinear forms
    std::map<RegionIdType, BaseBDBInt*> convectiveInts_;

    //! convective coefficient function for martial derivatve
    shared_ptr<CoefFunctionFormBased> convectiveCoef_;

    //! compute speed of sound, which may depend on temperature
    void ComputeSOS(PtrCoefFct& sos, PtrCoefFct dens, PtrCoefFct blk,
    		        PtrCoefFct regionTemp, std::string tempId);

    //! compute special PML coefficient: c0^2
    void ComputeSOS_SQR(PtrCoefFct& sosPML, PtrCoefFct dens, PtrCoefFct blk,
    		            PtrCoefFct regionTemp, std::string tempId);

    //! compute mech acou factor
    void CalcMechAcouFac(PtrCoefFct& coef, PtrCoefFct dens);
    void CalcMechAcouFacWithCoef(PtrCoefFct& coef, PtrCoefFct dens, PtrCoefFct surfDens, Double& scalFactor, Global::ComplexPart type) ;

    //! override from Single PDE due to convective operators
    virtual void FinalizePostProcResults();

    //! stores if the Acoustic PDE is coupled to mechanics
    bool isMechCoupled_;

    //! flag for transient PML
    bool isTimeDomPML_;

    //! flag indicating if we have almost PML (better stability in 3D)
    bool isAPML_;

    //! need wave-PDE for changing density
    bool complexFluidFormulation_;

    //! flag for checking if only one material is defined in the whole computational domain 
    bool isOnlyOneMaterial_;

    // Variables for the Time-Domain Equvivalent Fluid formulation
    bool timeDomainEqFluidFormulation_; //! flag indicating if the TDEF formulation is applied
    StdVector<bool> isTDEFReg_;         // flag indicating if a region is a TDEF region (i.e., has TDEF material definition)
    rationalCoeffsTDEF TDEFcoeffs_;     // struct with coeffs of the rational function approximation of the equivalent fluid parameters
    NumPolesTDEF TDEFpoleNumber_;       // // struct no. poles in the rational function approximation of the equivalent fluid parameters
    PtrCoefFct invTDEFBlk_;
    PtrCoefFct invTDEFDens_;

    //! Definition of convective integrators (Pierce Operator)
    //! \param actRegion  region id
    //! \param curRegNode current region node
    //! \param actSDList  list of elements
    //! \param coeffM     material coefficient of mass integrators
    template <UInt DIM, bool IS_COMPLEX>
    void DefineConvectiveIntegrators(RegionIdType actRegion, PtrParamNode curRegNode, shared_ptr<ElemList> actSDList, PtrCoefFct coeffM);

    //! This function defines Perfectly Matched Layer (PML) integrators for a given region.
    //! It outsources the definition of the stiffness and mass integrator for this region.
    //! \tparam DIM The dimension of the problem.
    //! \param actRegion The active region identifier.
    //! \param actSDList Shared pointer to the active element list.
    //! \param curRegNode Pointer to the current region node.
    //! \param c0 Shared pointer to the coefficient function c0.
    //! \param coeffK Shared pointer to the coefficient function for stiffness.
    //! \param coeffM Shared pointer to the coefficient function for mass.
    //! \param tempId Temporary identifier string.
    //! \param stiffInt Pointer to the base stiffness integrator (output parameter).
    //! \param massInt Pointer to the base mass integrator (output parameter).
    template <UInt DIM>
    void DefinePMLIntegrators(RegionIdType actRegion, shared_ptr<ElemList> &actSDList, PtrParamNode &curRegNode, PtrCoefFct &c0, PtrCoefFct &coeffK, PtrCoefFct &coeffM, std::string &tempId, BaseBDBInt *&stiffInt, BaseBDBInt *&massInt);

    //! This function assigns the integrator context for all defined integrators of the region.
    //! \param stiffInt Reference to a pointer to the stiffness integrator.
    //! \param massInt Reference to a pointer to the mass integrator.
    //! \param actRegion The active region identifier.
    //! \param actSDList Shared pointer to the active element list.
    //! \param coeffK Pointer to the coefficient function for stiffness.
    //! \param coeffM Pointer to the coefficient function for mass.
    void SetIntegratorContext(BaseBDBInt *&stiffInt, BaseBDBInt *&massInt, RegionIdType actRegion, shared_ptr<ElemList> &actSDList, PtrCoefFct &coeffK, PtrCoefFct &coeffM);
  };
}

#endif