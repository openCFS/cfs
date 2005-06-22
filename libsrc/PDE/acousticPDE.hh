#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "SinglePDE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"
 
namespace CoupledField {

  //! Class for acoustic equation (no adaptivity)


  class AcousticPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    AcousticPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile );

    //! Destructor
    virtual ~AcousticPDE(){};

	//! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation( Grid *aptgrid );

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators();
  
    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! perform postprocessing on data
    void PostProcess(){};

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteResultsInFile(const UInt kstep,
                                    const Double asteptime,
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0);

    //! return size of solution
    UInt getSize() const 
    { return size_;}

#ifdef ADAPTGRID
    //! test error of computation
    virtual Boolean TestError(const UInt level);
#endif


    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output);
  
    //! calculate the vector of coupling forces to the mechanical PDE
    void CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                              StdVector<UInt> & couplingNodes,
                              Vector<Double>& elemCouplingSols,
                              UInt couplingdof );
  
    //! 
    void SetMechanicCoupling() {
      isMechCoupled_ = TRUE;
    }
  


  protected:

    //! Init the time stepping
    void InitTimeStepping();

    //! indicator for mechanic coupling
    Boolean isMechCoupled_;

	//! total number of unknowns (equations)
    UInt size_;
	//! variable in which PDE is formulated
    SolutionType formulation_;

	//! list of boundaries( for absorbing BCs)
    StdVector<RegionIdType> absBCs_;
	//! switch for absorbing boundary conditions
    Boolean absorbingBCs_;                

    // solving of nonlinear acoustics
    NodeStoreSol<Double> sol_der1Array_, sol_der2Array_;
    Vector<Double> RhsLinVal_;

    // Postprocessing results
    NodeStoreSol<Double> solDeriv1_; //!< contains 1st derivative of solution
    NodeStoreSol<Double> solDeriv2_; //!< contains 2nd derivative of solution
    NodeStoreSol<Double> rhs_;

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;

    //! bubbledensity
    Double bubbleDensity_;

    Boolean plotRHS_; //!< Flag for saving of rhs for output


  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the acoustics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the AcousticPDE class that is set, if the keyword
    //! is specified.\n
    //! \todo Specification of ReadStoreResults for AcousticPDE!!!
    void ReadStoreResults();

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class AcousticPDE
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
