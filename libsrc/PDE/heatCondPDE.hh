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
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    HeatCondPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile );

    //! Destructor
    virtual ~HeatCondPDE(){};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! perform postprocessing on data
    void PostProcess(){;};

    //! write results in file
    //! \param kstep actual time step number
    //! \param asteptime time corresponding to kstep
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteResultsInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);

    //! write history results in file
    //! \param kstep actual time step number
    //! \param asteptime time corresponding to kstep
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteHistoryInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    void CalcOutputCoupling();


  protected:
    //! Init the time stepping
    void InitTimeStepping();

    // variation of inhomogeneous Neumann boundary condition    
    StdVector<Double> htc_; //!< heat transfer coefficient
    StdVector<Double> tSolid_; //!< initial temperature of body
    StdVector<Double> tFluid_; //!< ambient temperature 



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
