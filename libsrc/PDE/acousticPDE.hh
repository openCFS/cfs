#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "SinglePDE.hh"
#include "ODEDescr/KellerMiksis.hh"
#include "ODEDescr/Gilmore.hh"
#include "ODESolve/ODESolver_RKF45.hh"
#include "Domain/GridAdaption/GridAdaption.hh"
 
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
    void ReadDampingInformation();

    //! Read special boundary conditions (here: bubble information)
    void ReadSpecialBCs();

    //! Initialize NonLinearities
    void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();
  
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! perform postprocessing on data
    void PostProcess();

    //! calculate Force acting on specified surface elements
    void CalcForce( StdVector<Elem*> & saveElems );

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

    //! calculate the heat source term for heatConduction PDE
    void CalcHeatCouplingRHS( Vector<Double> & energy, 
                              StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                              UInt actCoupling, UInt numCouplingNodes );
  
    //! 
    void SetMechanicCoupling() {
      isMechCoupled_ = TRUE;
    }
  


  protected:

    //! Init the time stepping
    void InitTimeStepping();

    // ========================
    // set solution information
    // ========================    
    Boolean isMechCoupled_; //!< indicator for mechanic coupling
    
    UInt size_; //!< total number of unknowns (equations)

    SolutionType formulation_; //!< variable in which PDE is formulated

    StdVector<RegionIdType> absBCs_; //!< subdomains, which form absorbing BCs

    Boolean absorbingBCs_; //!< switch for absorbing BCs     
    
    Boolean fracDamping_; //!< switch indicating use of fractional damping
    
    //! switch for special bcs in combination with slicing technique
    Boolean m_bWriteSpecialBCs;


    // ========================
    // time stepping
    // ========================
    // solving of nonlinear acoustics
    NodeStoreSol<Double> sol_der1Array_, sol_der2Array_;
    Vector<Double> RhsLinVal_;


    // ========================
    // coupling
    // ========================
    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 


    // ========================
    // Postprocessing results
    // ========================
    NodeStoreSol<Double> solDeriv1_; //!< contains 1st derivative of solution
    NodeStoreSol<Double> solDeriv2_; //!< contains 2nd derivative of solution
    NodeStoreSol<Double> rhs_; //!< right hand side vector

    // force calculation on surface elements
    ElemStoreSol<Double> acouForce_; //!< contains force on surface elements
    Double sumAcouForce_; //!< contains force acting on all surface elements
    StdVector<std::string> saveElemForceHist_;//!< name of elements to be saved

    // calculate pressure from acoustic potential
    StdVector<RegionIdType> calcElemPressure_; //!< contains the regions
    ElemStoreSol<Double>  acouPressure_; //!< conatins acoustic pressure
    StdVector<std::string> saveElemPressureHist_;//!< name of elements

    //! Attribute describing model for bubble dynamics
    BubbleDynType bubbleDynType_;
    //! bubbledensity
    Double bubbleDensity_;

    Boolean plotRHS_; //!< Flag for saving of rhs for output

    // DODO
    // the grid adaption object
    GridAdaption *m_pGridAdaption;

  private:

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the acoustics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the AcousticPDE class that is set, if the
    //! keyword is specified.\n
    //! \todo Specification of ReadStoreResults for AcousticPDE!!!
    void ReadStoreResults();

    // reads in the PML data
    void ReadDataPML(std::string& typePML, Matrix<Double>& inner, 
		     Double& dampPML, RegionIdType actRegion);

    //! computes the PML layer dimensions
    void GetPMLLayerData(Matrix<Double>& inner, Matrix<Double>& outer,
			 UInt actSD);

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
