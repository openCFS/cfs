#ifndef ACOUSTICSPLITPDE_HH /// ADAPT
#define ACOUSTICSPLITPDE_HH /// ADAPT

#include "SinglePDE.hh"

namespace CoupledField{

  // forward class declaration	/// ADAPT ??
  class BaseResult;
  class ResultHandler;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;
  class CoefFunctionMulti;

  class AcousticSplitPDE : public SinglePDE{ /// ADAPT

  public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    AcousticSplitPDE( Grid* aGrid, PtrParamNode paramNode,	/// ADAPT
                 PtrParamNode infoNode,
                 shared_ptr<SimState> simState, Domain* domain );

    virtual ~AcousticSplitPDE(){}; /// ADAPT


    //! Return acoustic formulation. Can either be pressure or potential.
    SolutionType GetFormulation() const { return formulation_; }   /// ADAPT# scalarPotential


  protected:

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> >
    CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );	/// ADAPT
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //!  Define available primary results
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    
    //! Init the time stepping
    void InitTimeStepping();


    //! create transient Mapping integrators
    template<UInt DIM>	/// ADAPT# mapping
    void DefineTransientMapping(shared_ptr<ElemList> eList,std::string id);

  private:

    //! stores if the Split PDE is scalar or vector potential
    SolutionType formulation_;

    //! override from Single PDE due to convective operators
    virtual void FinalizePostProcResults();

  };

}

#endif
