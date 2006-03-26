#ifndef FILE_SINGLEPDE
#define FILE_SINGLEPDE


#include "PDE/StdPDE.hh"

#include <list>

#include "Driver/assemble.hh"
#include "Driver/MHassemble.hh"
#include "Domain/GridAdaption/GridAdaption.hh"

#include "DataInOut/Scripting/scriptable.hh"
#include "Utils/mathParser.hh"


namespace CoupledField
{
  // forward class declaration
  class SpaceErrorEstimator;
  class BasePairCoupling;
  
  //! Base class for all kinds of single field problems.

  class SinglePDE : public StdPDE, public Scriptable
  {
  
  public:

    // friend declaration
    friend class BasePairCoupling;

    Boolean BooleanComplexMaterialData_;
 
    void Init(UInt sequenceStep = 0,
              std::string  bcSequenceTag = "anyTag");
  
    // ---------------------- ***** --------------------------------

    //! destructor
    virtual ~SinglePDE();
  

    //! returns the load names
    StdVector<std::string>& GetLoadDom()
    {return assemble_->loadDom_;};
  
    //! returns the load dofs
    StdVector<std::string>& GetLoadDof()
    {return assemble_->loadDof_;};
  
    //! returns the load values
    StdVector<std::string>& GetLoadVals()
    {return assemble_->loadVals_;};
  
    //!returns the load functions
    StdVector<std::string>& GetLoadFncs()
    {return assemble_->fncname_loads_;};
  
    //! MpCCI gets the geometry
    virtual void PreparePDE4Computation() {;};
  
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! define algebraic system 
    void DefineAlgSys();
  
  
    // ======================================================
    // ADATPTIVITY SECTION
    // ======================================================
#ifdef ADAPTGRID
    //@{
    //! \name Methods used for performing adaptivity
  
    //! test error of calculation. return TRUE, if it is more then tolerance
    virtual Boolean TestError(const Integer level);
  
    //! refine mesh
    virtual void RefineMesh(const Integer level=0);

    // Does this method belong to postproc section?
    //! write information about relative error of calculation
    void WriteErrorInfo(WriteResults * ptmeshes);
  
    //@}
#endif

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling (only done once)
    virtual void InitCoupling(PDECoupling * Coupling) = 0;

    //! reset coupling counters and data (done after each timestep)
    virtual void ResetCoupling();

    //! Fill in input coupling terms
    virtual void CalcInputCoupling();
  
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling() = 0;

  
    // ======================================================
    // GET /SET  METHODS
    // ======================================================

    //! Activate the direct coupling
    void SetDirectCoupling ();
    
    //! get algsys identification tag of PDE
    PdeIdType GetPDEId()
    { return pdeId_; }

    //! set algebraic system object
    void SetAlgebraicSystem( BaseSystem *algSys);

    //! set solveStep object
    void SetSolveStep ( StdSolveStep *solveStep);

    //! set pointer to time stepping
    void SetTimeStepping(TimeStepping *timeStepping);

    //! return subtype
    virtual std::string GetSubType() {return subType_;}

    //! return number of restraints
    UInt GetNumRestraints( );
 
    //! set boundary condition
    //! \param atimestep         time step of claculation
    void SetBCs(const Double atimestep);

    //! Method for modifying an inhomogeneous boundary condition
    void SetIDBC( const std::string &name,
                  const std::string &dofString, 
                  const std::string &value, 
                  const std::string &phase,
                  const std::string &dynamics );
  
    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();
  
    // ======================================================
    // METHODS FOR ASSEMBLING
    // ======================================================
  
    //! specify type of system matrix for AlgebraicSystem
    void AssembleMatrices( );
  
    //! setup source term
    void AssembleSrcRHS( const Double time = 0.0 );
  
    //!  assemble a nonlinear RHS part
    void AssembleNLRHS( const Double time = 0.0 );
  
    //!  assemble a spring into the system matrix
    void AssembleSprings( const Double time = 0.0 );
  
    //! Initialize all matrices with nonlinear behavior
    void InitNonLinMatrices();
  
    //! constructes the matrix graph by providing to the algebraic system the element connectivities
    void SetupMatrixGraph();

    //! sets the actual frequency (just needed for harmonic analysis)
    void SetFrequency(Double actFreq);
    
    //! trigger the reassembling of the matrices
    void SetReassemble();

    // ======================================================
    // METHODS & MEMBERS FOR POST PROCESSING
    // ======================================================


    StdVector<RegionIdType> calcCharge_;
    ElemStoreSol<Double> charges_;
    ElemStoreSol<Complex> chargesComplex_;
    StdVector<RegionIdType> chargeNeighborRegion_;


  protected:

  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    SinglePDE(Grid *aptgrid, 
              WriteResults *aOutFile, 
              TimeFunc *aTimeFunc);

    //! private copy constructor
    SinglePDE & operator= (const StdPDE & myPDE) {
      Error( "Not implemented", __FILE__, __LINE__ );
      
      // For compiler
      return *this;
      ;}

   
    // ======================================================
    // INITIALIZATION METHODS
    // ======================================================
  
    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    virtual void ReadStoreResults() = 0;


    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( )=0;

    //! read damping information
    virtual void ReadDampingInformation( ){
      ENTER_FCN( "SinglePDE::ReadDampingInformation" );
    };
    
    //! read material data
    virtual void ReadMaterialData();

    //! read from config-file info about BCs
    void ReadBCs();

    // overloaded version of ReadBCs for special
    // boundary conditions in derived classes
    virtual void ReadSpecialBCs(){}

    //! Initialize NonLinearities
    virtual void InitNonLin(){};

#ifdef ADAPTGRID  
    //! ----------------- functions for adaptivity
    //! in this fnc we delete old pointer to Error-object & create new one
    void ConstructorError();
#endif
    
    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    void SaveSolution( const Double * ptSol, UInt size );
    void SaveSolution( const Complex * ptSol, UInt size );
    //@}


    // ======================================================
    // SCRIPTING SECTION
    // ======================================================
    //@{
    //! \name Scripting Methods
    
    //! Register scriptable functions
    virtual void RegisterFunctions();

    //! Wrapper for setting an inhomogeneous boundary condition
    void Wrap_IDBC();

    //! Wrapper for getting the nodal result
    void Wrap_GetValue();


    //@}
    
    // ======================================================
    // DATA SECTION
    // ======================================================

  
    // -----------------------------------------------------------------------
    // Storing information
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to storing information
    
    //! vector containing solutiontypes of PDE
    StdVector<SolutionType> solTypes_;

    //! vector containgin dofs of solutiontypes
    StdVector<UInt> solDofs_;

    //! TRUE, if at least one Result is calculated and written
    Boolean hasOutput_;

    //! TRUE, if solution should be written to result file
    Boolean saveSol_;

    //! TRUE, if RHS values should be written to result file
    Boolean saveRHSval_;

    //! TRUE, if first derivative of solution should be written to result file
    Boolean saveDeriv1_;

    //! TRUE, if second derivative of solution should be written to result file
    Boolean saveDeriv2_;

    //! TRUE, if solution should be written to history file
    Boolean saveSolHist_;

    //! TRUE, if special nodes should be written to be used for a 2nd run
    Boolean m_bReadSpecialBCs;
    
    // the grid adaption object
    GridAdaption *m_pGridAdaption;

    //! linear damping for slicing technique
    Double m_dStartDamping_;
    Double m_dPulseTime_;
    Double m_dPulseOffset_;

    //! TRUE, if RHS values should be written to history file
    Boolean saveRHSvalHist_;

    //! TRUE, if first derivative of solution should be written to history file
    Boolean saveDeriv1Hist_;

    //! TRUE, if second derivative of solution should be written to history file
    Boolean saveDeriv2Hist_;

    //! outputFormat for complex numbers
    ComplexFormat complexFormat_;  

    //@}

    // -----------------------------------------------------------------------
    // Adaptivity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to adaptivity

    //! object with methods for error estimation
    SpaceErrorEstimator * ptError_;

    //! array where  we store the number of refinement for the element
    StdVector<Double> markingElems_;

    Vector<Double> errorMap_;  //!< array with error map
    Double tolSpaceErr_;       //!< tolerance
    //@}

    // -----------------------------------------------------------------------
    // Solver parameters
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to parameters for solver
    UInt maxnumiter_;       //!< maximum of iterations (for iterative solver)
    Double eps_;               //!< accuracy
    Double dampiter_;          //!< damping parameter within iterative solution
    Double coarsealpha_;       //!< coarsening factor (just for AMG)

    //@}

    // -----------------------------------------------------------------------
    // Miscellaneous paramters
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellanous paramters
  
    //! Identifier of the PDE which is used in the algebraic system
    PdeIdType  pdeId_;
  
    //! flag for direct coupling
    Boolean isDirectCoupled_;

    //! Handler for MathParser object
    MathParser::HandlerType mHandler_;

    //@}
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SinglePDE
  //! 
  //! \purpose 
  //! This class serves as base class for all single field problems, 
  //! like electrostatic,  acoustic, mechanic and others.
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
