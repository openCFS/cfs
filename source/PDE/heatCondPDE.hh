#ifndef FILE_HEATCONDPDE_2001
#define FILE_HEATCONDPDE_2001

#include "SinglePDE.hh"
 
namespace CoupledField {

  //! Class for heat conduction equation
  class HeatCondPDE: public SinglePDE
  {

  public:

    //!  Constructor. Here we read integration parameters.
    /*!
      \param aGrid pointer to grid
    */
    HeatCondPDE( Grid* aGrid, ParamNode* paramNode );

    //! Destructor
    virtual ~HeatCondPDE(){};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

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
    struct InhomHeatNeumannBc : public InhomNeumannBc {

      //! Heat transfer coefficient
      std::string htc;

      //! Initial temperature of body
      std::string tSolid;

      //! Ambient temperature
      std::string tFluid;
    };

    typedef  StdVector<shared_ptr<InhomHeatNeumannBc> > HeatInBcList;

    //! special neumann boundary conditions
    HeatInBcList heatInBcs_;
    
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
