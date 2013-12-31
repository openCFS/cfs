// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FLUIDMECHPDE
#define FILE_FLUIDMECHPDE

#include "SinglePDE.hh" 

namespace CoupledField
{

  // forward class declaration
  class BaseResult;
  class LinearFormContext;
  class CoefFunctionMulti;
  
  //! Class for linearized perturbed formulation of Navier-Stoke's equations.
  class FlowPDE : public SinglePDE {

  public:

    //! Constructor
    /*!
      \param grid pointer to grid
      \param paramNode pointer to the corresponding parameter node
    */
    FlowPDE( Grid* grid, PtrParamNode paramNode,
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~FlowPDE(){};

  protected:

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this PDE
    void DefineIntegrators( );

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators() {
      EXCEPTION("ncInterfaces are not implemented for FlowPDE");
    }

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators(){};

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping();

    //! Nothing to do
//    void SetTimeStep(const Double dt) {;};

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

    BaseBDBInt *
    GetStiffIntegrator( BaseMaterial* actSDMat,
                        RegionIdType regionId,
                        bool isComplex );

    //! Surface regions on which pressure surface integral has to be defined
    StdVector<RegionIdType> presSurfaces_;

    //! type of stabilization
    bool stabilizedBochev_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class FlowPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! the incompressibel Navier-Stoke's equations in 3D.
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

#endif // FILE_FLUIDMECHPDE
