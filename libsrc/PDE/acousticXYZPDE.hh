#ifndef FILE_ACOUSTICXYZPDE_2005
#define FILE_ACOUSTICXYZPDE_2005

#include "SinglePDE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "Domain/GridAdaption/GridAdaption.hh"
 
namespace CoupledField {

  //! Class for acoustic equation (no adaptivity)


  class AcousticXYZPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    AcousticXYZPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile );

    //! Destructor
    virtual ~AcousticXYZPDE(){};

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();
  
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! perform postprocessing on data
    void PostProcess();

    //! calculate pressure from acoustic potential
    void CalcElemPressure();

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

#ifdef ADAPTGRID
    //! test error of computation
    virtual Boolean TestError(const UInt level);
#endif


    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling)
    {Error("InitCoupling not implemented",__FILE__,__LINE__);};
  
    //! calculate coupling terms
    void CalcOutputCoupling()
    {Error("CalcOutputCoupling not implemented",__FILE__,__LINE__);};

    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output)
    {Error("HasOutput not implemented",__FILE__,__LINE__);};


  protected:

    //! Init the time stepping
    void InitTimeStepping();

    // ========================
    // time stepping
    // ========================
    // solving of nonlinear acoustics
    NodeStoreSol<Double> sol_der1Array_, sol_der2Array_;

    // ========================
    // Postprocessing results
    // ========================

    BaseNodeStoreSol * press_;    //!< pressure

    NodeStoreSol<Double> solDeriv1_; //!< contains 1st derivative of solution
    NodeStoreSol<Double> solDeriv2_; //!< contains 2nd derivative of solution

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the acoustics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the AcousticXYZPDE class that is set, if the
    //! keyword is specified.\n
    //! \todo Specification of ReadStoreResults for AcousticXYZPDE!!!
    void ReadStoreResults();

    // reads in the PML data
    void ReadDataPML(std::string& typePML, Matrix<Double>& inner, 
		     Double& dampPML, RegionIdType actRegion);

    //! computes the PML layer dimensions
    void GetPMLLayerData(Matrix<Double>& inner, Matrix<Double>& outer,
			 UInt actSD);

    //!
    Boolean savePress_;

    //!
    Boolean savePressHist_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcousticXYZPDE
  //! 
  //! \purpose 
  //! This class is derived from class SinglePDE.
  //! It is used for solving acoustic equation on one time step.  
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
