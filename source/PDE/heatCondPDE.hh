// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HEATCONDPDE_2001
#define FILE_HEATCONDPDE_2001

#include "SinglePDE.hh"
 
namespace CoupledField {

  class LinearFormContext;

  //! Class for heat conduction equation
  class HeatCondPDE: public SinglePDE
  {

  public:

    //!  Constructor. Here we read integration parameters.
    /*!
      \param aGrid pointer to grid
    */
    HeatCondPDE( Grid* aGrid, PtrParamNode paramNode );

    //! Destructor
    virtual ~HeatCondPDE(){};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    /** region load 
     * @param regionLoads as returned from ReadRegionLoadsFromXML
     * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefineRegionLoadIntegrators(std::map<RegionIdType, RegionLoad>& regionLoads, StdVector<LinearFormContext*>* linForms = NULL);

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! perform postprocessing on data
    void CalcResults( shared_ptr<BaseResult> result );

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    void CalcOutputCoupling();
    
    //! Set the flag for a thermo-electrostatic coupling, this
    //! means to add the thermo integrator block into the damping 
    //! matrix and set off the thermo integrator block from the 
    //! mass matrix
    void SetElectroCoupling();

    
    //! Set the flag for a thermo-mechanical coupling, this
    //! means to add the thermo integrator block into the damping 
    //! matrix and set off the thermo integrator block from the 
    //! mass matrix
    void SetMechCoupling();
    
    
    //! Initialize all the nodes by this value
    void SetInitialCondition();
    
    
  protected:
	  
    //! Define availabe result types
    void DefineAvailResults();

    //! Init the time stepping
    void InitTimeStepping();
    

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the heat conduction PDE 
    //! are given in the following table. Here 'Keyword' and 'Result Type'
    //! refer to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the HeatCondPDE class that is set, if the
    //! keyword is specified.\n
    //! \todo Specification of ReadStoreResults for HeatCondPDE!!!
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
    
    //! flag for thermoelectro-coupling
    bool isElectroCoupled_;
    
    //! flag for thermoelastic-coupling
    bool isMechCoupled_;
    
    NonLinMethodType nonLinMethod_;

    //Double InitialTemp_;

    //! vector containing regionIds of non-conforming interfaces
    StdVector<RegionIdType> ncIFaces_;
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class HeatCondPDE
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
