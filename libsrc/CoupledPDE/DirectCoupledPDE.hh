#ifndef FILE_DIRECT_COUPLED_PDE
#define FILE_DIRECT_COUPLED_PDE

#include "PDE/StdPDE.hh"

namespace CoupledField {

  // forward class declaration
  class SinglePDE;
  class BasePairCoupling;
  class Assemble;
  class BaseSystem;
  
  //! This class implements the direct coupling of StdPDEs.

  //! This class handles the direct coupling of two or more
  //! SinglePDEs. Therefore it needs references to SinglePDEs
  //! and the according pair coupling objects.

  class DirectCoupledPDE : public StdPDE
  {
  public:
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    DirectCoupledPDE(Grid *aptgrid, 
                     WriteResults *aOutFile, 
                     TimeFunc *aTimeFunc);
    
    //! Destructor
    virtual ~DirectCoupledPDE();
  
    //! Set SinglePDEs
    void SetSinglePDEs( const StdVector<SinglePDE*> &pdes);

    //! Set Coupling objects
    void SetCouplings( const StdVector<BasePairCoupling*> &couplings);

    //! Initialization routine
    void Init(UInt sequenceStep = 0,
              std::string  bcSequenceTag = "anyTag");

    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();

  
    //! set boundary condition
    //! \param atimestep         time step of claculation
    void SetBCs(const Double atimestep);

    //! define algebraic system 
    void DefineAlgSys();

    //! return number of restraints
    UInt GetNumRestraints();
  
    //! Get types of needed matrices (sysmtem, stiffness,..)
    void GetMatrixTypes( std::set<FEMatrixType> &matTypes);

    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    void SaveSolution( const Double * ptSol, UInt size );
    void SaveSolution( const Complex * ptSol, UInt size );
    //@}

    //!
    virtual void InitTimeStepping();
    

    // ======================================================
    // POSTPROC SECTION
    // ======================================================
  
    //@{
    //! \name Methods performing post-processing
  
    //! Do Postprocessing as descriped in conf file
    void PostProcess();
  
    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteResultsInFile(const UInt kstep = 0,
                            const Double asteptime = 0.0,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);
    //@}


    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling (only done once)
    void InitCoupling(PDECoupling * Coupling);

    //! reset coupling counters and data (done after each timestep)
    void ResetCoupling();
  
    //! Fill in input coupling terms
    void CalcInputCoupling();
  
  
    //! calculate coupling terms
    void CalcOutputCoupling();

  
    // ======================================================
    // METHODS FOR ASSEMBLING
    // ======================================================
  
    //! specify type of system matrix for AlgebraicSystem
    void AssembleMatrices();
  
    //! setup source term
    void AssembleSrcRHS(const Double time=0);
  
    //!  assemble a nonlinear RHS part
    void AssembleNLRHS(const Double time=0);
  
    //!  assemble a spring into the system matrix
    void AssembleSprings(const Double time=0);
  
    //! Initialize all matrices with nonlinear behavior
    void InitNonLinMatrices();
  
    //! constructes the matrix graph by providing to the algebraic system the element connectivities
    void SetupMatrixGraph();

  private:
  
    //! define the SolutionStep-Driver
    virtual void DefineSolveStep();

    //! References to SinglePDEs
    StdVector<SinglePDE*> singlePDEs_;

    //! References to pairwise coupling objects
    StdVector<BasePairCoupling*> couplings_;

    //! Number of total unknowns 

    //! Number of total unknowns, which is the sum of each PDE's
    //! number of equations times its equations dof
    UInt totalUnknowns_;

  };

} // end of namespace

#endif
