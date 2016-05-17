// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HEATCONDPDE_2001
#define FILE_HEATCONDPDE_2001

#include "SinglePDE.hh"
 

namespace CoupledField {

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class linElecInt;
  class LinearFormContext;
  class BaseBDBInt;
  class ResultFunctor;

  //! Class for heat conduction equation
  class HeatPDE: public SinglePDE
  {

  public:

    //!  Constructor. Here we read integration parameters.
    /*!
      \param aGrid pointer to grid
    */
    HeatPDE( Grid* aGrid, PtrParamNode paramNode,
             PtrParamNode infoNode,
             shared_ptr<SimState> simState, Domain* domain  );

    //! Destructor
    virtual ~HeatPDE(){};

    //! Is heat source (RHS) definition driven by interface between solid and void?
    inline bool HasInterfaceDrivenRHS() { return interfaceDrivenHeatSource_; }

  protected:
    
    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! Initialize NonLinearities
    virtual void InitNonLin();

   //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){};

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();
    
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //!  Define available postprocessing results
    void DefinePrimaryResults();

    //! Define available postprocessing results
    void DefinePostProcResults();

    //! Init the time stepping
    void InitTimeStepping();
    
    SolutionType GetNativeSolutionType() const { return HEAT_TEMPERATURE; }

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the heat conduction PDE 
    //! are given in the following table. Here 'Keyword' and 'Result Type'
    //! refer to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the HeatPDE class that is set, if the
    //! keyword is specified.\n
    //! \todo Specification of ReadStoreResults for HeatPDE!!!
    void ReadStoreResults();

    //! Extend public class of inhom. Neumann boundary condition in order
    //! account for special temperature coefficients
    struct RobinBc : public HomDirichletBc {

      //! Heat transfer coefficient
      std::string HTC;

      //! Initial temperature of the surrounding fluid
      std::string bulkTemp;
    };

    typedef  StdVector<shared_ptr<RobinBc> > RobinBcList;

    //! special neumann boundary conditions
    RobinBcList robinBcs_;
    
    bool interfaceDrivenHeatSource_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class HeatPDE
  //! 
  //! \purpose 
  //! This class is derived from class SinglePDE.
  //! It is used for solving heat conduction equation on one time step.  
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
