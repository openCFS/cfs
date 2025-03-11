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
  
  class BaseFeFunction;
  class SimState;
  class Domain;
  class MathParser;

  //! Base class for partial differential equations
  class BasePDE
  {

  public:
    
    //! Constructor
    BasePDE( PtrParamNode paramNode, PtrParamNode infoNode, 
             shared_ptr<SimState> simState, Domain* domain );
    
    //! Destructor
    virtual ~BasePDE();

    //! write general defines (BCs, loads, etc.) to info-file
    virtual void WriteGeneralPDEdefines() = 0;

    //! Return object for assembling the element matrices and vectors
    virtual Assemble * GetAssemble() = 0;

    //! Return pointer to the SolveStep object
    virtual BaseSolveStep * GetSolveStep() = 0;
    
    //! Return pointer to domain
    Domain * GetDomain()  {
      return domain_;
    }

    //! Return all Rhs FeFunctions
    virtual
    std::map<SolutionType, shared_ptr<BaseFeFunction> >& GetRhsFeFunctions() {
    	EXCEPTION("GetRhsFeFunctions not implemented in Base PDE");
    }

    //! Update PDE due to updated step in multistep solution strategy
    virtual void UpdateToSolStrategy() = 0;

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
    typedef enum {NO_ANALYSIS, STATIC, TRANSIENT, HARMONIC, MULTIHARMONIC,
    			  EIGENFREQUENCY, INVERSESOURCE, MULTI_SEQUENCE, BUCKLING, EIGENVALUE } AnalysisType;

    //! Enums for advection-diffusion stabilisation 
    typedef enum {NO_STABILISATION, ARTIFICIAL_DIFFUSION, SUPG } StabilisationType;
    
    /** Helper method which determines if an AnalyisType is complex. */
    bool IsComplex();
                  
    /** This enum is for the string/enum conversion of AnalysisType.
     * Don't be confused with the actual analysis type in  analysis_.
     * BasePDE::SetEnums() has to be called once before! */              
    static Enum<AnalysisType> analysisType;

    static Enum<StabilisationType> stabilisationType;

    /** Sets up the Enums */
    static void SetEnums();

    /** Obtain the matrix derivative map i.e. which fematrix is associated to
     *  which derivative. in general we would have
     *  Stiffness -> zero order
     *  Damping -> First Order
     *  Mass -> Second Order
     *  Auxiliary -> PDE dependent but mostly stiffness
     *  So we hardcode this here but it can be overwritten by the PDE
     */
    virtual std::map<FEMatrixType,Integer> GetMatrixDerivativeMap();

    virtual PtrParamNode GetMyParam() {
    		return myParam_;
    }

    virtual void SetSourceApproxType( bool type ) {
    	approxSourceWithDeltaFnc_ = type;
    }

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
    
    //! Info ParameterNode of current pde
    PtrParamNode myInfo_;
    //@}

    //! Name of the PDE
    std::string pdename_;
    
    //! Pointer for saving the internal simulation state
    shared_ptr<SimState> simState_;
    
    //! Pointer to domain
    Domain * domain_;
    
    //! Pointer to MathParser
    MathParser * mp_;
    
    //!
    bool approxSourceWithDeltaFnc_;

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

