/** DOCUMENTATION HEADER
 * This PDE is used to calculate the magnetostatic field via
 * the magnetic scalar potential using the following weak form
 * 
 * (gradPsi,gradW)_Omega = 0 with Psi = Psi_0 at Gamma_e .
 * 
 */
#ifndef FILE_MAGNETICSCALARPOTENTIALADJPDE_NEW
#define FILE_MAGNETICSCALARPOTENTIALADJPDE_NEW

#include "SinglePDE.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>

namespace CoupledField
{
  
  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class linElecInt;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;
  
  //! Class for magnetostatic scalar potential PDE
  class MagneticScalarPotentialAdjPDE : public SinglePDE {
    
  public:

    //! Constructor
    /*!
     \param grid pointer to grid
     \param paramNode pointer to the corresponding parameter node
    */
    MagneticScalarPotentialAdjPDE( Grid* grid, PtrParamNode paramNode,
                                PtrParamNode infoNode,
                                shared_ptr<SimState> simState, Domain* domain );
    
    //! Destructor
    virtual ~MagneticScalarPotentialAdjPDE(){};
       
    //! Initialize NonLinearities
    virtual void InitNonLin();

  protected:

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
    PtrParamNode infoNode );
    
    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );
    
    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators(){}; // curly braces are here, beacuse in the .cc file is no implementation of the method!
    
    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators(){}; // curly braces are here, beacuse in the .cc file is no implementation of the method!
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators(); 
    
    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    void FinalizePostProcResults();
    
    //! Init the time stepping:
    void InitTimeStepping();

    //! map containing the magnetic source field, e.g., from a previous sequence step
    std::map<RegionIdType, PtrCoefFct> hPostproc_;

    //! Coefficient function, containing the overall permeability
    shared_ptr<CoefFunctionMulti> perm_;
   
    //stores the flux for hystersis and nonlinear models
    std::map<RegionIdType, shared_ptr<CoefFunctionMulti> >  nlFluxCoefm_;

    //! map containing the magnetic field intensity of forward simulation, in case of design parameters
    std::map<RegionIdType, PtrCoefFct> hPostprocParam_;

    //! map containing the derivative of the anhysteresis curve
    std::map<RegionIdType, PtrCoefFct> magAnhystDerivA_;
    std::map<RegionIdType, PtrCoefFct> magAnhystDerivPs_;

    //yes, when permeability depends on parameters
    bool muDerivParam_;
  };
  
#ifdef DOXYGEN_DETAILED_DOC
  
  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================
  
  //! \class MagneticScalarPotentialADJPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! electrostatic equation in 3D. 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 
  
#endif
  
} // end of namespace

#endif
