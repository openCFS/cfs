// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LINFLUIDMECHPDE
#define FILE_LINFLUIDMECHPDE

#include "SinglePDE.hh" 

namespace CoupledField
{

  // forward class declaration
  class BaseResult;
  class LinearFormContext;
  class CoefFunctionMulti;
  
  //! Class for linearized perturbed formulation of Navier-Stoke's equations.
  class LinFlowPDE : public SinglePDE {

  public:

    //! Constructor
    /*!
      \param grid pointer to grid
      \param paramNode pointer to the corresponding parameter node
    */
    LinFlowPDE( Grid* grid, PtrParamNode paramNode,
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~LinFlowPDE(){};

  protected:

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this PDE
    void DefineIntegrators( );

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators() {
      EXCEPTION("ncInterfaces are not implemented for LinFlowPDE");
    }

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators(){};

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping() {;};

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! Calculate field variables at arbitrary points
    void CalcField( SolutionType solType, StdVector<const Elem*>& elems,
                    StdVector<LocPoint>& points, SingleVector& values );

    //! SubType of electrostatic section
    std::string subType_;

    // *****************
    //  POSTPROCESSING
    // *****************
    
    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
     
    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );


  private:

    //! create feFunction for meanFluidMech velocity
    void CreateMeanFlowFunction(StdVector<std::string> dofNames);
    
    BaseBDBInt *
    GetStiffIntegrator( BaseMaterial* actSDMat,
                        RegionIdType regionId,
                        bool isComplex );

    //! Coefficient function for the flow field
    
    //! This coefficient function describes the flow field. As this 
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> meanFlowCoef_;
    shared_ptr<CoefFunction> meanFlowCoefScattered_;
    shared_ptr<BaseFeFunction> meanFlowFeFct_;
    shared_ptr<FeSpace> meanFlowFeSpace_;

    //! Surface regions on which pressure surface integral has to be defined
    StdVector<RegionIdType> presSurfaces_;

    //! compressible formulation
    bool isCompressible_;

    //! considers a varying background flow
    bool enableC2_;

    //! amplitude factor for first convective term in case of varying background flow
    Double factorC1_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class LinFlowPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! a linearized perturbed formulation of Navier-Stoke's equations in 3D.
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

} // end of namespace CoupledField

#endif // FILE_FLUIDMECHPERTURBEDPDE
