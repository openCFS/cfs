#ifndef FILE_NRBCPDE_2001
#define FILE_NRBCPDE_2001

#include "SinglePDE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "Domain/GridAdaption/GridAdaption.hh"
 
namespace CoupledField {

  //! Class for nrbc equation


  class nrbcPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read the pde parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    nrbcPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile,
            std::string pdeNameWithIndex, StdVector<SolutionType> localsolType);

    //! Destructor
    virtual ~nrbcPDE(){};

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
    //! write history results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteHistoryInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);


#ifdef ADAPTGRID
    //! test error of computation
    virtual bool TestError(const UInt level);
#endif


    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);
  
    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);
  
    //! calculate the vector of coupling RHS term to the acoustic PDE
    void CalcAcouCouplingRHS( StdVector<Elem*> * couplingElems, 
                              StdVector<UInt> & couplingNodes,
                              Vector<Double>& elemCouplingSols,
                              UInt couplingdof );
  
    //! 
    void SetAcouCoupling() {
      isAcouCoupled_ = true;
    }
  


  protected:

    //! Init the time stepping
    void InitTimeStepping();

    //! order of absorbtion
    UInt order_;

    //! indicator for mechanic coupling
    bool isAcouCoupled_;
    
    //! variable in which PDE is formulated
    SolutionType formulation_;
    
    //! list of boundaries( for absorbing BCs)
    StdVector<RegionIdType> absBCs_;
    //! switch for absorbing boundary conditions
    bool absorbingBCs_;                

    bool m_bWriteSpecialBCs;            //!< switch for special bcs in combination with slicing technique

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

    bool plotRHS_; //!< Flag for saving of rhs for output

    // DODO
    // the grid adaption object
    GridAdaption *m_pGridAdaption;

    UInt indexofPDE_;
    char * auxpdeName_;

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the acoustics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the nrbcPDE class that is set, if the keyword
    //! is specified.\n
    //! \todo Specification of ReadStoreResults for nrbcPDE!!!
    void ReadStoreResults();

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class nrbcPDE
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
