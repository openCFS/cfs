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

  class BasePDE
  {

  public:

    bool converged_; //!< needed for coupling with MpCCI

    //! Constructor
    BasePDE();
    
    //! Destructor
    virtual ~BasePDE();

    //! write general defines (BCs, loads, etc.) to info-file
    virtual void WriteGeneralPDEdefines() = 0;

    //!
    virtual Assemble * getPDE_assemble() = 0;

    //! Return pointer to the SolveStep object
    virtual BaseSolveStep * GetSolveStep() = 0;


    //! write the PDE state (pdememento) to a restart file "simname_pdename.restart"
    virtual void WriteRestart(const UInt nstep, UInt totalUnknowns = 0) = 0;

    //! read the PDE state (pdememento)from a restart file: "simname_pdename.restart"
    virtual void ReadRestart(UInt &startStep, UInt totalUnknowns = 0) = 0;

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //@{
    //! \name Methods performing post-processing

    //! Do Postprocessing as descriped in conf file
    virtual void PostProcess() = 0;

    //! write results in file
    //! \param stepOffset offset for starting (time/frequency)step
    //! \param timeOffset offset for starting time / frequency
    virtual void WriteResultsInFile(const UInt kstep, 
                                    const Double actTimeFreq, 
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0) = 0; 

    //! write history (selected nodes/elements) results in file
    //! \param stepOffset offset for starting (time/frequency)step
    //! \param timeOffset offset for starting time / frequency
    virtual void WriteHistoryInFile(const UInt kstep, 
                                    const Double actTimeFreq, 
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0) = 0; 
    //@}

    
    // ======================================================
    // SET /GET METHODS
    // ======================================================

    //@{
    //! \name Getter methods

    //! return name of pde
    virtual std::string GetName() {return pdename_;}
    //@}

   
  protected:
    
    //! define the SolutionStep-Driver
    virtual void DefineSolveStep() {
      Error( "DefineSolveStep not implemented", __FILE__, __LINE__ );
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
    UInt bcSequenceIndex_;
    //@}

    //! name of the PDE
    std::string pdename_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BasePDE
  //! 
  //! \purpose
  //! Class BasePDE is the base class from which different types of classes
  //! describing individual types of PDEs are derived. This class reflects the 
  //! basic functionality of all solvable problems (single, direct-coupled, 
  //! iterative coupled).
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

