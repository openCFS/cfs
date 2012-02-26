#ifndef FILE_STDPDE
#define FILE_STDPDE

#include <fstream>
#include "PDE/BasePDE.hh"

#include <set>

#include "Materials/Composite.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField {


  // forward class declarations
  class PDECoupling;
  class BasePairCoupling;
  class WriteResults;
  class BaseNodeStoreSol;
  class StdSolveStep;
  class PDECoupling;
  class ParamNode;
  class BiotSavart;
  class BaseFeFunction;
  
  //! Base class for all single-field and direct-coupled problems

  class StdPDE : public BasePDE
  {
  
  public:

    // friend cStdlass declarations
    friend class PDECoupling;

    // public typedefs
    typedef StdVector<shared_ptr<ResultInfo> > ResultInfoList;
    //! typedefs for result handling
    typedef std::set<shared_ptr<ResultInfo> > ResultSet;


    //! Virtual destructor
    virtual ~StdPDE();
    
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! Define the algebraic system 
    virtual void DefineAlgSys();

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();
    
    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling (only done once)
    virtual void InitCoupling(PDECoupling * Coupling) = 0;

    //! reset coupling counters and data (done after each timestep)
    virtual void ResetCoupling() = 0;
  
    //! Fill in input coupling terms
    virtual void CalcInputCoupling() = 0;
  
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling() = 0;


    // ======================================================
    // GET/SET METHODS
    // ======================================================

    //! Returns the feFunction which holds a result related to the specified solutionType
    virtual shared_ptr<BaseFeFunction> GetFeFunction( SolutionType solType);
    
    //! Return all solution FeFunctions
    virtual  
    std::map<SolutionType, shared_ptr<BaseFeFunction> > GetFeFunctions( ) {
      return feFunctions_;
    }
    
    //! Return all Rhs FeFunctions
    virtual 
    std::map<SolutionType, shared_ptr<BaseFeFunction> > GetRhsFeFunctions() {
      return rhsFeFunctions_;
    }
    
    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();
    
    //! Return vector with resultInfo types
    ResultInfoList& GetResultInfos() { return results_;}

    //! Return result for given solutionType

    //! Returns the resultInfo related to the specified solutionType
    virtual shared_ptr<ResultInfo> GetResultInfo( SolutionType solType );


    virtual AnalysisType GetAnalysisType() {
      return analysistype_;
    }
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output)
    {
      EXCEPTION("not implemented");
      return false;
    }
  
    //! Also for fractional damping model do obtain
    virtual UInt GetFracMemory() {
      return fracMemory_;
    }
    
    virtual bool GetFracDamping() {
      return fracDamping_;
    }

    virtual bool GetIsaxi() {
      return isaxi_;
    }
    
    //! Return Biot-Savart object
    shared_ptr<BiotSavart> GetBiotSavart() {
      return biotSavart_;
    };

    //! check, if Biot-Savart is set
    bool IsBiotSavart() {
      return isBiotSavart_;
    };

    //! Return pointer to paramNode of current pde
    PtrParamNode GetParamNode() { return myParam_; }

    //!
    //! \for computing and adding RHS to PDE in case of special sources 
    virtual void ComputeRHS(const Double atime) {;};

    //!
    //! \for computing vortex source both analytically and with complex 
    //! \potential function 
    virtual  void VortexAnalytical(Double & press, Vector<Double>& dTij_di,
                                   const Double x,
                                   const Double y, const Double t, 
                                   const UInt outType){
      EXCEPTION("VortexAnalytical is only implemented in acouFlowNoisePDE");};
  
    //! set boundary condition OBSOLETE in the future
    virtual  void SetBCs() = 0;

    virtual  void InitStabParams( ) {
      EXCEPTION("Not Implemented, only needed for fluidMech and smooth");
    };

    virtual  void PrintStabParams( ) {
      EXCEPTION("Not Implemented, only needed for fluidMech");
    };


    
    //! Init the time stepping
    virtual void InitTimeStepping()
    {EXCEPTION("InitTimeStepping not implemented");};

    //! Get couplings object
    virtual StdVector<BasePairCoupling*>* GetCouplingsObject() 
    {
      return NULL;
    };

    virtual void AcouSourceCalc(){EXCEPTION("AcouSourceCalc not implemented");};

    // ======================================================
    // COMMUNICATION ROUTINES FOR PARAMETER IDENTIFICATION
    // ======================================================

    //@{


    std::map<RegionIdType, BaseMaterial*>  getPDEMaterialData()
    {return materials_;};
    
    Assemble * getPDE_assemble(){return assemble_;}

    StdVector<RegionIdType> getPDE_subdoms(){return subdoms_;}

    AlgebraicSys * getPDE_algsys(){return algsys_;}
    
    Grid * getPDE_grid() {return ptgrid_;}

    // set if PDE is nonlinear
    virtual void SetNonLinearity(bool nonLin){
      nonLin_=nonLin;};

    // set if PDE is nonlinear (material dependency)
    virtual void SetMaterialNonLinearity(bool nonLin){
      nonLinMaterial_=nonLin;};

    //@}

    //@{
    //!  Get functions concerning nonlinearity

    bool IsNonLin() 
    { return nonLin_;};

    bool IsNonLinMaterial() 
    { return nonLinMaterial_;};

    bool IsHysteresis() 
    { return isHysteresis_;};

    bool IsIterCoupled() 
    { return isIterCoupled_;};

    bool& IsFirstTimeStepStatic()
    { return firstTimeStepStatic_;};

    std::map<RegionIdType, StdVector<NonLinType> >& GetNonLinRegionTypes() 
    { return regionNonLinTypes_;};

    UInt& GetIterCoupledCounter() 
    { return iterCoupledCounter_;};

    PDECoupling* GetCoupling()
    {return ptCoupling_;};

    //! List of inhomogeneous Dirichlet boundary conditions
    IdBcList GetIDBCList(){
      return idBcs_;};

    //! List of inhomogeneous Dirichlet boundary conditions read from file
    IdFileBcList GetIdFiBCList(){
      return idFiBcs_;};
    
    //! List of inhomogeneous Neumann boundary conditions
    const InBcList& GetINBCList() { return inBcs_; }

    //! Return material class
    MaterialClass GetMaterialClass() const { return pdematerialclass_; }
    //@}
      
  protected:
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    StdPDE(Grid *aptgrid, PtrParamNode paramNode );
  
    //! private copy constructor
    StdPDE & operator= (const StdPDE & myPDE) {
      EXCEPTION("Not implemented");

      // The following line is only to satisfy the compiler
      return *this;
    };
  

    //! Get coefficient for damping matrix in fractional damping model
    //! \todo This function has to be removed when the fractional
    //! damping model gets implemented in a separate Forms-class
    virtual Double GetFracDampMatrixCoeff(RegionIdType region);

    //! retruns, if PDE needs to store previous solution or not
    virtual bool NeedsPrevSol() {
      return  needSolPrev_;
    }

    // ======================================================
    // DATA SECTION
    // ======================================================
  
    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
    //@{
    //! \name Attributes related to geometry and node numbering

    Grid * ptgrid_;        //!< pointer to grid object

    //! subdomain-levels belonging to PDE
    StdVector<RegionIdType> subdoms_;

    //@}

    /** This is our pde info node. To be set/overwritten in each PDE! */ 
    PtrParamNode infoNode_; 
    
    //! Parameter node for OLAS
    PtrParamNode olasNode_;

    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering 
    //! defines subtype of mechanic PDE: plainStrain, 3d, ...
    std::string subType_;
    //@}

      // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------    

    //@{
    //! \name boundary conitions

    //! Homogeneous Dirichlet boundary coniditions
    HdBcList hdBcs_;

    //! Inhomogeneous Dirichlet boundary conditions
    IdBcList idBcs_;
    
    //! Inhomogeneous Dirichlet boundary conditions
    IdFileBcList idFiBcs_;
    
    //! List of inhomogeneous Neumann boundary conditions
    InBcList inBcs_;

    //! List of constraints
    ConstraintList constraints_;

    //! Right hand side load definitions
    LoadList loads_;

    //! Number of additional in. dirichlet boundary equations due to coupling
    UInt numCouplingBcs_;

    //@}

    // -----------------------------------------------------------------------
    // Non-linearity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to nonlinearity
    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    bool isHysteresis_;     //!< flag for hysteresis

    //! map for each region the type of nonlinearity
    std::map<RegionIdType, StdVector<NonLinType> > regionNonLinTypes_;

    //! map for each nonlinearity the id
    std::map<std::string, NonLinType> nonLinTypes_;

    // type of nonlinear algorithm (e.g., Newton)
    NonLinMethodType nonLinMethod_;

    //@}


    // -----------------------------------------------------------------------
    // Material data
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes handling info on material data

    //! Maps regions and (simple) materials
    std::map<RegionIdType, BaseMaterial*> materials_;    
    
    //! Maps regions and composite materials
    std::map<RegionIdType, Composite> compositeMaterials_;

    //! material class
    MaterialClass pdematerialclass_;   

    //@}


    // -----------------------------------------------------------------------
    // PDE coupling
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to handling PDE coupling
    bool isIterCoupled_;        //!< PDE couples with others
    Vector<Double> matParam_;      //!< change to material parameter
    bool updateCouplingBCs_ ;  //!< flag if coupling BC were already set
    PDECoupling *ptCoupling_;     //!< pointer to coupling object
  
    //! nodes at which coupling terms are calculated
    std::list<UInt> couplingNodes;    
  
    //! elements at which coupling terms are calculated
    StdVector<Elem*> couplingElements;
  
    //! iteration counter for coupled PDE solution process
    UInt iterCoupledCounter_;
    //@}


    // -----------------------------------------------------------------------
    // Time stepping
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to time stepping
    bool diagMass_;           //!< use of diagonal mass matrix in explicit time stepping
    bool firstTimeStepStatic_; //!< needed for coupled, iterative methods
    //@}



    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    //! Vector containing the results calculated by this PDE
    //! OBSOLETE
    ResultInfoList results_;
    
    //! Set containing the types of possible results
    ResultSet availResults_;
    //!

    //! flag indicating if this PDE needs the algebraic system
    bool needsAlgsys_;

    AnalysisType analysistype_; //!< analysis type
    bool isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    UInt dim_;                  //!< space dimension of pde
    bool isaxi_;             //!< true: axisymmetric problem
    bool isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    bool needSolPrev_;          //! true, if solution at time step n has also to bve stored

    //! list of damping types for all regions
    std::map<RegionIdType,DampingType> dampingList_;
    
    //! use of complex material data per region
    std::map<RegionIdType,bool> complexMatData_;

    bool fracDamping_; //!< true: fractional damping model
    UInt fracMemory_;     //!< number of old time steps to be saved (for fractional damping)
    
    //! object, handling the computation by Biot-Savart fundamental field
    shared_ptr<BiotSavart> biotSavart_;
        
    //! excitation computed by Biot-Savart
    bool isBiotSavart_;

    //! type of interpolation (for fractional damping)
    InterpolType inType_;
    
    //! checks, if we have for the coupling a incremental solution
    bool isIncrFormulation_;    
    
    //! if yes, PDE is computed on deformed geometry
    bool updatedLagrangeForm_;

    //! flag for knowing if we have to call ComputeRHS() in the harmonic driver
    bool ComputeRHSforHarm_;    

    //! Pointer to object of analysis (Static, Trans, Harm or Eig)
    Assemble * assemble_;

    //! pointer to SolveStep classes
    StdSolveStep * solveStep_;
    
    AlgebraicSys * algsys_;      //!< pointer to algebraic system
    
    //! Pointer to solution strategy object
    shared_ptr<SolStrategy> solStrat_;
  
    /** This is the node for linear system responsible for this pde. */
    PtrParamNode olasInfo_;
    
    //@}

    //! map which associates a Postprocessing result to its corresponding Primary Result
    std::map<SolutionType,SolutionType> postProcResults_;

    //! Map Storing FeSpaces for each unknown of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions_;
    
    //! Map storing the feFunctions of the previous step
    std::map<SolutionType, shared_ptr<BaseFeFunction> > prevFeFunctions_;

    //! Map storing time derivatives of FeFunctions
    std::map<SolutionType, shared_ptr<BaseFeFunction> > timeDerivFeFunctions_;
    
    //! Map storing the feFunctions of the RHS
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions_;

  }; // class StdPDE

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class StdPDE
  //! 
  //! \purpose 
  //! This class serves as base class for all single-field and 
  //! direct coupled problems. The idea is, that an iterative coupled
  //! pde can have one or more objects of StdPDE, which are then
  //! solved iteratively.
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
