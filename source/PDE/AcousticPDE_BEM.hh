#ifndef ACOUSTICPDE_BEM_HH
#define ACOUSTICPDE_BEM_HH

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

  class AcousticPDE_BEM : public SinglePDE
  {
    public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    AcousticPDE_BEM( Grid* aGrid, PtrParamNode paramNode,
                 PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain );

    virtual ~AcousticPDE_BEM(){};

    int dbg = 1;

    //! Indicate that acoustic PDE takes part in coupling to mechanics.
    void SetMechanicCoupling() {
      isMechCoupled_ = false;
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

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //!  Define available primary results
    void DefinePrimaryResults();
    
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

//    //! Set special RHS values
//    virtual void SetRhsValues();

  private:

    //! stores if the Acoustic PDE is in potential or pressure form
    SolutionType formulation_;

    //! if true, speed of sound squared ist at Laplace term!
    bool sosAtLaplace_;

    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;
    

    //! acoustic source density
    shared_ptr<CoefFunctionMulti> acousticSourceDensityCoef_;

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
  };
}

#endif
