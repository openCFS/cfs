#ifndef FILE_BASEPDE
#define FILE_BASEPDE

#include <list>

#include "Utils/StdVector.hh"
#include "General/environment.hh"


namespace CoupledField
{


  // forward class declarations
  class BaseSolveStep;
  class Assemble;

  //! Base class for partial differential equations

  //! Class BasePDE is the base class from which different types of classes
  //! describing individual types of PDEs are derived. This class reflects the basic
  //! functionality of all solvable problems (single, direct-coupled, iterative coupled).

  class BasePDE
  {

  public:

    //! Constructor
    BasePDE();
    
    //! Destructor
    virtual ~BasePDE();
    
    //! Initializes PDE 
    
    //! Initializes the PDE. This function is only called one time.
    //! \param bcSequenceId (input) name of the tag for current set of 
    //! boundary condition
    // virtual void Init(Integer sequenceStep = 0,
    //                std::string  bcSequenceTag = "anyTag") 
    //     {
    //       Error( "Init not implemented here");
    //   }

    //! write general defines (BCs, loads, etc.) to info-file
    virtual void WriteGeneralPDEdefines() = 0;

    //!
    virtual Assemble * getPDE_assemble() = 0;

    //! Return pointer to the SolveStep object
    virtual BaseSolveStep * GetSolveStep() = 0;

    //! set time step
    //! \params dt Current time step
    virtual void SetTimeStep(const Double dt) = 0;

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //@{
    //! \name Methods performing post-processing

    //! Do Postprocessing as descriped in conf file
    virtual void PostProcess() = 0;

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const Integer kstep = 0,
                                    const Double asteptime = 0.0, 
                                    Integer stepOffset = 0,
                                    Double timeOffset = 0.0) = 0; 
    //@}

    
    // ======================================================
    // GETTER METHODS
    // ======================================================

    //@{
    //! \name Getter methods

    //! return name of pde
    virtual std::string GetName() {return pdename_;}
    //@}

   
  protected:
    
    //! define the SolutionStep-Driver
    virtual void DefineSolveStep() {
      Error("DefineSolveStep not implemented");
    }

    
    // ======================================================
    // DATA SECTION
    // ======================================================
    
    //@{
    //! \name Attributes connected to the handling of boundary conditions

    //! tag of current set of boundary conditions

    //! tag of current set of boundary conditions. For a multiSequence-simulation
    //! this id determines, which set of boundary conditions is applied.
    std::string bcSequenceTag_;

    //! index of current set of boundary conditions. For a multiSequence-simulation
    //! this index determines, which set of boundary conditions is applied.
    Integer bcSequenceIndex_;
    //@}

    //! name of the PDE
    std::string pdename_;



  };

  

} // end of namespace

#endif

