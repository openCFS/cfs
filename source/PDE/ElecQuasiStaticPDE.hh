// ================================================================================================
/*!
 *       \file     ElecQuasiStaticPDE.hh
 *       \brief    This PDE models electric current conduction and electric displacement
 *                 current with the scalar electric potential
 *
 *       \date     January, 2021
 *       \author   Manfred Kaltenbacher, TU Graz
 */
//================================================================================================

#ifndef FILE_ELEQUASISTATIC_NEW
#define FILE_ELEQUASISTATIC_NEW

#include "SinglePDE.hh" 


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
  
  //! Class for electric current conduction PDE
  class ElecQuasistaticPDE : public SinglePDE {

  public:

    //! Constructor
    /*!
      \param grid pointer to grid
      \param paramNode pointer to the corresponding parameter node
    */
    ElecQuasistaticPDE( Grid* grid, PtrParamNode paramNode,
             PtrParamNode infoNode,
             shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~ElecQuasistaticPDE(){};

  protected:

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );
    
    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    virtual void DefineSurfaceIntegrators(){};

    //! Define all RHS linearforms for load / excitation
    void DefineRhsLoadIntegrators();

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping();

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    // ======================================================
    // COUPLING SECTION
    // ======================================================

  public:
    
  protected:
    
    //! SubType of electrostatic section
    std::string subType_;
    
    // *****************
    //  POSTPROCESSING
    // *****************
    
    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    // //! \copydoc SinglePDE::FinalizePostProcResults
    // void FinalizePostProcResults();

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

  private:

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class ElecCurrentPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! electric current conduction equation in 3D. 
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
