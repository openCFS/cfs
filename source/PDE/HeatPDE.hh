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

    //! set coupling to LinFlow PDE and whteher to use a symmetric form
    void SetLinFlowPDECouplingFlags(bool useSymmetricForm);

    /** constants for test temperature gradients (rhs), used for homogenization. We depend on the int values! */
    typedef enum {X=0, Y=1, Z=2} TestStrain;

    //! Is heat source (RHS) definition driven by interface between solid and void?
    inline bool HasInterfaceDrivenRHS() { return interfaceDrivenHeatSource_; }

    /** Add the integrators for the test strains for homogenization to the linear forms, similar as in multiple load case;
     * called from Excitation::ReadLoads or Excitation::SetHomogenizationTestStrains() (optimization)
     * @param test is an enum
     * @param linForms set to append linear Forms to, if NULL use assemble_
     * very similar to MechPDE::DefineTestStrainIntegrator*/
    void DefineTestStrainIntegrator(const TestStrain test, StdVector<LinearFormContext*>* linForms = NULL);

    /** Stores test temperature gradients for asymptotic homogenization*/
    static Enum<TestStrain> testStrain;

  protected:
    
    //! Read special boundary conditions
    void ReadSpecialBCs(){};

    /** Check what stabilisation type is selected for the given region
     * and compute the stabilisation factor CoefFunction.
     * @param ent region to check
     * @param factor pointer to write the CoefFunction of the stabilisation factor to
     * @return the stabilisation type
     */
    CoupledField::BasePDE::StabilisationType GetStabilisation(RegionIdType ent, PtrCoefFct &factor);

    /** Check if volume neighbour has SUPG enabled. */
    bool VolNeighbourHasSUPG(EntityIterator entit);

    void ReadDampingInformation();

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
    void DefineSurfaceIntegrators(){};

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
    
    //! Incorporate heat transfer boundary conditions
    void HeatTransferBC();

    //! Incorporate thermal radiation boundary conditions
    void ThermalRadiationBC();


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
      Double HTC;

      //! Initial temperature of the surrounding fluid
      Double bulkTemp;
    };

    typedef  StdVector<shared_ptr<RobinBc> > RobinBcList;

    //! special neumann boundary conditions
    RobinBcList robinBcs_;
    
    //! flag for heat source on design boundary, need this for tracking optimization
    bool interfaceDrivenHeatSource_;

    //! Tensor type
    SubTensorType tensorType_;

    //! Coefficient function for the convective velocity

    //! This coefficient function describes the velocity field. As this
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> convecVelCoef_;

    //! store convective bilinear forms
    std::map<RegionIdType, BaseBDBInt*> convectiveInts_;

    //! true, if coupled to LinFlow PDE
    bool isLinFlowPDECoupled_;
    //! Whether to use a symmetric formulation or not in coupling (to LinFlowPDE)
    bool isCouplingFormulationSymmetric_;

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
