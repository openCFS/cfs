// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEPDE
#define FILE_BASEPDE

#include <list>

#include "Utils/StdVector.hh"
#include "General/environment.hh"
#include "General/exception.hh"

namespace CoupledField
{


  // forward class declarations
  class BaseSolveStep;
  class Assemble;
  class ParamNode;

  //! Base class for partial differential equations

  class BasePDE
  {

  public:

    bool converged_; //!< needed for coupling with MpCCI

    //! Constructor
    BasePDE( ParamNode* paramNode );
    
    //! Destructor
    virtual ~BasePDE();

    //! write general defines (BCs, loads, etc.) to info-file
    virtual void WriteGeneralPDEdefines() = 0;

    //!
    virtual Assemble * getPDE_assemble() = 0;

    //! Return pointer to the SolveStep object
    virtual BaseSolveStep * GetSolveStep() = 0;


    //! write the PDE state (pdememento) to a restart file "simname_pdename.restart"
    virtual void WriteRestart( ) = 0;

    //! read the PDE state (pdememento)from a restart file: "simname_pdename.restart"
    virtual void ReadRestart(UInt &startStep ) = 0;

    //! perform cleanup and do last computations
    virtual void Finalize() {};

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //@{
    //! \name Methods performing post-processing

    //! write results in file
    //! \param stepOffset offset for starting (time/frequency)step
    //! \param timeOffset offset for starting time / frequency
    virtual void WriteResultsInFile(const UInt kstep, 
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
      EXCEPTION("DefineSolveStep not implemented");
    }

    
    // ======================================================
    // DATA SECTION
    // ======================================================
    
    //@{

    //! index of current set of boundary conditions. For a multiSequence-simulation
    //! this index determines, which set of boundary conditions is applied.
    UInt sequenceStep_;

    //! ParamNode of current pde
    ParamNode * myParam_;
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

