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

    //! returns all single PDEs
    StdVector<SinglePDE*>& GetSinglePDEs()
    { return singlePDEs_;};

    // get total number of unknowns
    UInt GetTotalUnknowns()
    { return totalUnknowns_;};

    //! Set Coupling objects
    void SetCouplings( const StdVector<BasePairCoupling*> &couplings);

    //! Get couplings object
    StdVector<BasePairCoupling*>& GetCouplingsObject() 
    { return couplings_;};

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
  
    //! write a restart file "simname_pdename.restart"
    void WriteRestart(const UInt nstep, UInt totalUnknowns=0);

    //! read a restart file "simname_pdename.restart"
    void ReadRestart(UInt &startStep, UInt totalUnknowns=0);

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteResultsInFile(const UInt kstep = 0,
                            const Double asteptime = 0.0,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);

    //! write history of nodes/elements to file
    void WriteHistoryInFile(const UInt kstep = 0,
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

  private:

    //! Auxilliary method for building re-ordering into BaseEQN objects

    //! This auxilliary method is used to obtain the permutation vector
    //! for the equation numbers from the algebraic system and pass it on
    //! to the BaseEQN data objects associated with the individual directly
    //! coupled PDEs. The algebraic system passes the pointer and also the
    //! responsibility for handling the memory associated with the
    //! permutation vector to CFS. After the BaseEQN data object has
    //! incorporated the new ordering, we should, thus, delete the
    //! permutation vector.
    //! \note The correct place to call this method depends on the type
    //! of algebraic system employed. In the case of a StandardSystem this
    //! method must be called only after assemble of the graph has completed.
    //! In the case of an SBM_System the method must be called <b>after</b>
    //! all PDEs have assembled their sub-graphs and <b>before</b> the
    //! coupling objects start the assembly of their sub-graphs.
    void IncorporateReordering();

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

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class DirectCoupledPDE
  //! 
  //! \purpose This class is implementes the concept of directly coupling two
  //! or more SinglePDE's together by the use of BasePairCoupling 
  //! objects. It has the same interface as a StdPDE, so the solveStep-classes
  //! can work transparently with it.
  //! 
  //! \collab This class is instantiated in the class Domain. It gets one or 
  //! more instances of SinglePDE (which get assembled on a main diagonal 
  //! position in the resulting algebraic system) and the according instances
  //! of BasePairCoupling, which implement the pairwise coupling between each 
  //! two SinglePDEs.
  //! 
  //! \implement
  //! 
  //! \status In use
  //! 
  //! \unused
  //! 
  //! \improve 
  //! - At the moment only a Newmark timestepping scheme is applied. In some
  //! cases this will not be sufficient.
  //! - Harmonic coupling has not been tested yet, as well as nonlinear 
  //! analysis
  //!

#endif

} // end of namespace

#endif
