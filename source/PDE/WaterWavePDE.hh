#ifndef WATERWAVEPDE_HH
#define WATERWAVEPDE_HH

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

  class WaterWavePDE : public SinglePDE{

  public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    WaterWavePDE( Grid* aGrid, PtrParamNode paramNode,
                 PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain );

    virtual ~WaterWavePDE(){};

    //! Indicate that acoustic PDE takes part in coupling to mechanics.
    void SetMechanicCoupling() {
      isMechCoupled_ = true;
    }

  protected:

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> >
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators() {;};

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

    //! create transient PML integrators
    template<UInt DIM>
    void DefineTransientPMLInts(shared_ptr<ElemList> eList,std::string id);

  private:

    //! stores if the Acoustic PDE is in potential or pressure form
    SolutionType formulation_;

    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;
    
    //! override from Single PDE due to convective operators
    virtual void FinalizePostProcResults();

    //! stores if the Acoustic PDE is coupled to mechanics
    bool isMechCoupled_;

    //! flag for transient PML
    bool isTimeDomPML_;

    //! flag indicating if we have almost PML (better stability in 3D)
    bool isAPML_;

    //! true for modeling surface gravity waves
    bool isSurfaceGravityWave_;
  };

}

#endif
