#ifndef FILE_SINGLEPDE
#define FILE_SINGLEPDE


#include "PDE/StdPDE.hh"

#include <list>

#include "Domain/GridAdaption/GridAdaption.hh"
#include "DataInOut/Scripting/scriptable.hh"
#include "Utils/mathParser.hh"
#include "Domain/resultDof.hh"
#include "Utils/elemstoresol.hh"
#include "DataInOut/timefunc.hh"
#include "Domain/bcs.hh"

namespace CoupledField
{
  // forward class declaration
  class SpaceErrorEstimator;
  class BasePairCoupling;
  class DirectCoupledPDE;
  class Assemble;
  
  //! Base class for all kinds of single field problems.

  class SinglePDE : public StdPDE, public Scriptable
  {
  
  public:

    // friend declaration
    friend class BasePairCoupling;
    friend class DirectCoupledPDE;

    bool boolComplexMaterialData_;
 
    virtual void Init(UInt sequenceStep = 0,
                      std::string  bcSequenceTag = "anyTag");
  
    // ---------------------- ***** --------------------------------

    //! destructor
    virtual ~SinglePDE();
  
    //! MpCCI gets the geometry
    virtual void PreparePDE4Computation() {;};
  
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! define algebraic system 
    virtual void DefineAlgSys();
 
  
    // ======================================================
    // ADATPTIVITY SECTION
    // ======================================================
#ifdef ADAPTGRID
    //@{
    //! \name Methods used for performing adaptivity
  
    //! test error of calculation. return true, if it is more then tolerance
    virtual bool TestError(const Integer level);
  
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

    //! get algsys identification tag of PDE
    PdeIdType GetPDEId()
    { return pdeId_; }

    //! return subtype
    virtual std::string GetSubType() {return subType_;}

    //! Set Direct coupling information
    virtual void SetDirectCoupling();

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
    // METHODS & MEMBERS FOR POST PROCESSING
    // ======================================================

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
    
    //! true, if at least one Result is calculated and written
    bool hasOutput_;

    //! true, if solution should be written to result file
    bool saveSol_;

    //! true, if RHS values should be written to result file
    bool saveRHSval_;

    //! true, if first derivative of solution should be written to result file
    bool saveDeriv1_;

    //! true, if second derivative of solution should be written to result file
    bool saveDeriv2_;

    //! true, if solution should be written to history file
    bool saveSolHist_;

    //! true, if RHS values should be written to history file
    bool saveRHSvalHist_;

    //! true, if first derivative of solution should be written to history file
    bool saveDeriv1Hist_;

    //! true, if second derivative of solution should be written to history file
    bool saveDeriv2Hist_;

    //! outputFormat for complex numbers
    ComplexFormat complexFormat_;  

    //@}
    
    // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------    

    //@{
    //! \name boundary conitions

    //! Homogeneous Dirichlet boundary coniditions
    StdVector<HomDirichletBc> hdBc_;

    //! Inhomogeneous Dirichlet boundary conditions
    StdVector<InhomDirichletBc> idBC_;
    
    //! Inhomogeneous Neumann boundary conditions
    StdVector<InhomNeumannBc> inBC_;

    //! Right hand side load definitions
    LoadList loads_;

    //! Number of additional in. dirichlet boundary equations due to coupling
    UInt numCouplingBcs_;

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
    // Miscellaneous paramters
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellanous paramters
  
    //! Identifier of the PDE which is used in the algebraic system
    PdeIdType  pdeId_;
  
    //! flag for direct coupling
    bool isDirectCoupled_;

    //! Handle for MathParser object
    MathParser::HandleType mHandle_;

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
