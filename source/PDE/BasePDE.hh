#ifndef FILE_BASEPDE
#define FILE_BASEPDE

#include <list>

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "General/Exception.hh"
#include "General/Enum.hh"
#include "Forms/IntScheme.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{


  // forward class declarations
  class BaseSolveStep;
  class Assemble;
  class ParamNode;
  class BaseFeFunction;


  //! Base class for partial differential equations

  class BasePDE
  {

  public:
    
    bool converged_; //!< needed for coupling with MpCCI

    //! Constructor
    BasePDE( PtrParamNode paramNode );
    
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

    //! Set solution step in case of mulitlevel solution
   virtual void SetSolutionStep(UInt solStep_ ) {};
   
    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //@{
    //! \name Methods performing post-processing

    //! write results in file
    virtual void WriteResultsInFile(const UInt kstep, 
                                    const Double actTimeFreq ) = 0; 

    //@}

    
    // ======================================================
    // SET /GET METHODS
    // ======================================================

    //@{
    //! \name Getter methods

    /** Name of the PDE, overwritten in the coupled PDE */
    virtual const std::string& GetName() const { return pdename_; }
    //@}

    /** do string/enum conversion via BasePDE::analysisType */
    typedef enum {STATIC, TRANSIENT, HARMONIC, EIGENFREQUENCY, 
                  MULTI_SEQUENCE } AnalysisType;
    
    /** Helper method which determines if an AnalyisType is complex. */
    static bool IsComplex(AnalysisType type);
                  
    /** This enum is for the string/enum conversion of AnalysisType.
     * Don't be confused with the actual analysis type in  analysis_.
     * BasePDE::SetEnums() has to be called once before! */              
    static Enum<AnalysisType> analysisType;

    /** Sets up the Enums */
    static void SetEnums();


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
    PtrParamNode myParam_;
    //@}

    //! name of the PDE
    std::string pdename_;
    
    //! Solution strategy for problem
    SolStrategyType solStrategy_;
    
    //! In case of several multilevel solution, this variable carries the step
    
    //! This variable indicates the step in case of a two- or mulit-level
    //! solution strategy. It is 1 based.
    UInt solStep_;


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

