// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STDPDE
#define FILE_STDPDE
#include <fstream>
#include "PDE/basePDE.hh"

#include <set>

#include "PDE/timestepping.hh"
#include "Domain/Composite.hh"
#include "Elements/fefunction.hh"
#include "Elements/fespace.hh"

namespace CoupledField {


  // forward class declarations
  class PDECoupling;
  class WriteResults;
  class TimeStepping;
  class BaseNodeStoreSol;
  class StdSolveStep;
  class PDECoupling;
  class ParamNode;
  class InfoNode;
  class BaseFeFunction;
  
  //! Base class for all single-field and direct-coupled problems

  class StdPDE : public BasePDE {
  
  public:

    // friend cStdlass declarations
    friend class PDECoupling;

    // public typedefs
    typedef StdVector<shared_ptr<ResultInfo> > ResultInfoList;


    //! Virtual destructor
    virtual ~StdPDE();
    
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! Define the algebraic system 
    virtual void DefineAlgSys() = 0;

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();
    
    /** Finds the xml node of the linear system
     * @return might be NULL */
    ParamNode* FindLinearSystem(const std::string& sysName);
    
    //! Transfer parameters from CFS++ to OLAS parameter object

    //! This method reads the parameters specified for the linear system
    //! associated with the PDE from the parameter object of CFS++, adapts
    //! them using an expert module and passes them to OLAS. It relies on
    //! The SetParams() method of the CFSOLASParams class for doing this.
    //! \param sysName name of the linear system in the XML file (in the
    //!                case of a single PDE this coincides with pdename_)
    void ReadOlasParams( std::string sysName );


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
    //! and a entityList of the specified name 
    virtual shared_ptr<BaseFeFunction> GetFeFunction( SolutionType solType, std::string name);

    //! Returns the feFunction which holds a result related to the specified solutionType
    //! and a entityList of the specified name 
    virtual shared_ptr<BaseFeFunction> GetFeFunction( SolutionType solType, RegionIdType regID);

    //! Returns the feFunction which holds a result related to the specified solutionType
    //! and a entityList of the specified name 
    StdVector<FunctionDescription> GetFunctionDescriptors(SolutionType type);

    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();

    //! Return vector with resultInfo types
    ResultInfoList& GetResultInfos() { return results_;}

    //! Return result for given solutionType

    //! Returns the resultInfo related to the specified solutionType
    virtual shared_ptr<ResultInfo> GetResultInfo( SolutionType solType );

    //! return pointer to vector with subdomains, on which we calculate the PDE
    virtual StdVector<RegionIdType> * getSDsPDE()
    { return &subdoms_;}

    virtual AnalysisType GetAnalysisType() {
      return analysistype_;
    }
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output)
    {
      EXCEPTION("not implemented");
      return false;
    }
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS1() const;
  
    //! return pointer to vector with second derivative of solution
    virtual const Vector<Double> & getS2() const;

    //! return pointer to vector with last solution
    virtual const Vector<Double> & getOld1() const;
    

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

    //! Return pointer to paramNode of current pde
    ParamNode * GetParamNode() { return myParam_; }

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

    //! Transforms a given BoundaryCondition value according to Timestepping (i.e. TransientSim)
    virtual void TransformBC(Double& transVal, Double initValue, Integer eqnNumber) = 0;

    virtual  void InitStabParams( ) {
      EXCEPTION("Not Implemented, only needed for fluidMech and smooth");
    };

    virtual  void PrintStabParams( ) {
      EXCEPTION("Not Implemented, only needed for fluidMech");
    };

    
    //! Initialize all/some the nodes by this value
    virtual void SetInitialCondition(){;};
    bool IsSetInitialCondition() { return isSetInitialCondition_;};
    double getInitialCondition(){return InitialCondition_;}; 
    
    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    virtual void SaveSolution( const Double * ptSol, UInt size ) = 0;
    virtual void SaveSolution( const Complex * ptSol, UInt size ) = 0;
    virtual void SavePrevSolution( const Double * ptSol, UInt size ) = 0;
    //@}

    //@{
    //! Save load part of RHS to private variable
    virtual void SaveRHS( const Double * ptSol, UInt size ) = 0;
    virtual void SaveRHS( const Complex * ptSol, UInt size ) = 0;
    //@}

    /** Get the data vector of the current solution of the algebraic system */
    virtual SingleVector* GetSolutionVector() { return solVec_;} 

    /** Get the data vector of the current rhs of the algebraic system. */
    SingleVector* GetRHSVector() { return rhsVec_;}

    //! get the data vector of the previous solution of a PDE.
    virtual SingleVector * GetPrevSolutionVector();
    
    /// returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement( Vector<Double>& sol, const EntityIterator& it, 
                             shared_ptr<ResultInfo> res );
    void GetSolVecOfElement( Vector<Complex>& sol, const EntityIterator& it,
                             shared_ptr<ResultInfo> res );
   
    
    /// returns the vector of time derivative of the solution belonging 
    /// to all nodes of the actual element
    void GetDerivSolVecOfElement( Vector<Double>& sol, const EntityIterator& it,
                                  shared_ptr<ResultInfo> res );
    void GetDerivSolVecOfElement( Vector<Complex>& sol, const EntityIterator& it,
                                  shared_ptr<ResultInfo> res );
    
    /// returns the vector of 2nd time derivative of the solution belonging to all nodes 
    /// of the actual element
    void GetDeriv2SolVecOfElement( Vector<Double>& sol, const EntityIterator& it,
                                   shared_ptr<ResultInfo> res );
    void GetDeriv2SolVecOfElement( Vector<Complex>& sol, const EntityIterator& it, 
                                   shared_ptr<ResultInfo> res);

    //! Init the time stepping
    virtual void InitTimeStepping()
    {EXCEPTION("InitTimeStepping not implemented");};
    
    virtual void AcouSourceCalc(){EXCEPTION("AcouSourceCalc not implemented");};

    // ======================================================
    // COMMUNICATION ROUTINES FOR PARAMETER IDENTIFICATION
    // ======================================================

    //@{
    //!  The following methods are used only durig parameter
    //!  identification process! Maybe one day a more to CFS++ consistent 
    //!  nomenclature would be nice ...

    //shared_ptr<EqnMap> GetEqnMap() { return eqnMap_; }
    //
    

    std::map<RegionIdType, BaseMaterial*>  getPDEMaterialData()
    {return materials_;};
    
    BaseNodeStoreSol * getPDESolution() {return sol_;};

    BaseNodeStoreSol * getPDESolutionPrev() {return solPrev_;};
    
    BaseSystem * getPDE_algsys(){return algsys_;};
  
    UInt getPDE_numElems(){return numElems_;};

    UInt GetNumPdeEquations(){return numPdeEquations_;};

    UInt GetNumPdeUnknowns(){return numPdeUnknowns_;};

    UInt GetNumInHomBcs(){return numPdeInHomDirBc_;};
  
    UInt getPDE_spaceDim(){return dim_;};
  
    Grid * getPDE_grid(){return ptgrid_;};  
  
    Assemble * getPDE_assemble(){return assemble_;}
  
    StdVector<RegionIdType> getPDE_subdoms(){return subdoms_;}
   
    Vector<Complex> getPDE_complexValuedCharge()
    {return complexValuedCharge_;};
    
    void setPDE_MatDataType(Global::ComplexPart & pMatType){
      matDataType_ = pMatType;};
  
    Global::ComplexPart getPDE_MatDataType()
    {return matDataType_;}

    //! get pointer to time stepping object
    TimeStepping * getTimeStepping() { return TS_alg_;};

    void setPDE_complexValuedCharge(Vector<Complex> chargeVec)
    {complexValuedCharge_=chargeVec;};

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

    bool IsComputeRHS4HarmSet()
    { return ComputeRHSforHarm_;};

    bool IsInstationary()  
    { return isInstationary_;};


    Double GetRhsL2Norm(Vector<Double>& actRHS) 
    { return RhsL2Norm(actRHS);};

    std::map<RegionIdType, NonLinType>& GetNonLinRegionTypes() 
    { return regionNonLinType_;};

    UInt& GetIterCoupledCounter() 
    { return iterCoupledCounter_;};

    PDECoupling* GetCoupling()
    {return ptCoupling_;};

    //! List of inhomogeneous Dirichlet boundary conditions
    IdBcList GetIDBCList(){
      return idBcs_;};

    //@}
      
  protected:
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    StdPDE(Grid *aptgrid, ParamNode* paramNode );
  
    //! private copy constructor
    StdPDE & operator= (const StdPDE & myPDE) {
      EXCEPTION("Not implemented");

      // The following line is only to satisfy the compiler
      return *this;
    };
  


    //! stores an algsys_ vector into a StdVector
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    //! calculates L2-norm of RHS regarding entries due to penalty formulation
    virtual Double RhsL2Norm(Vector<Double>& stdVec);

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
    InfoNode* infoNode_; 

    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering 
    UInt numPdeEquations_;  //!< Overall Number of Equations in this PDE 
    UInt numPdeUnknowns_;  //!< Overall number of Free(!) equations for this PDE
    UInt numPdeInHomDirBc_;  //!< overall number of inhomogenious BCs for this PDE
    UInt numElems_;     //!< number of elements in subdomains
  
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
    std::map<RegionIdType, NonLinType> regionNonLinType_;

    //! map for each region the id of the nonlinearity
    std::map<RegionIdType, std::string> regionNonLinId_;

    //! map for each id the nonlinearity
    std::map<std::string, NonLinType> nonLinIdType_;
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
  
    //! Data Type which decides wheather material is real or complex
    Global::ComplexPart  matDataType_;
    //! contains element results of complex valued charge 
    Vector<Complex> complexValuedCharge_;
    Vector<Complex> complexValuedEfield_;
    //@}


    // -----------------------------------------------------------------------
    // PDE coupling
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to handling PDE coupling
    bool isIterCoupled_;        //!< PDE couples with others
    Vector<Double> matParam_;      //!< change to material parameter
    bool updateCouplingBCs_ ;  //!< flag if coupling BC were already set
    UInt couplingBCsCounter_;  //!< counter for number of coupling BCs
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
    TimeStepping * TS_alg_;       //!< handles the time stepping
    bool effectiveMass_;       //!< use effective mass formulation for transient analysis
    bool diagMass_;           //!< use of diagonal mass matrix in explicit time stepping
    bool firstTimeStepStatic_; //!< needed for coupled, iterative methods
    bool isInstationary_;    //!< flag for stationary/instationary PDEs (like fluidMech)

    //@}



    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    //! Vector containing the results calculated by this PDE
    //! OBSOLETE
    ResultInfoList results_;

    //!

    //! flag indicating if this PDE needs the algebraic system
    bool needsAlgsys_;

    AnalysisType analysistype_; //!< analysis type
    bool isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    UInt dim_;                  //!< space dimension of pde
    bool isaxi_;             //!< true: axisymmetric problem
    bool isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    //shared_ptr<EqnMap> eqnMap_; //!< new equation handling
    bool needSolPrev_;          //! true, if solution at time step n has also to bve stored

    BaseNodeStoreSol * sol_;    //!< solution
    SingleVector * solVec_;        //! needed in iterative coupled computation 
    SingleVector * rhsVec_;        //! needed when writing the RHS to file
    BaseNodeStoreSol * solPrev_;    //!< solution at time step n
    SingleVector * solVecPrev_;        //! needed in coupled computation 

    // ======================================================
    // FeSpace and FeFunction
    // ======================================================
    // @{
    // ! \name FeSpace and Function of the Unknown

    ////! FeSpace Vector associated with the PDE
    //StdVector<shared_ptr<FeSpace>> feSpaces_;

    ////! FeFunction of the unknown
    //shared_ptr<BaseFeFunction> feFunction_;
    //@}
    
    //! list of damping types for all regions
    std::map<RegionIdType,DampingType> dampingList_;

    bool fracDamping_; //!< true: fractional damping model
    UInt fracMemory_;     //!< number of old time steps to be saved (for fractional damping)
    


    //! type of interpolation (for fractional damping)
    InterpolType inType_;
    
    //! checks, if we have for the coupling a incremental solution
    bool isIncrFormulation_;    
    
    //! flag for knowing if we have to call ComputeRHS() in the harmonic driver
    bool ComputeRHSforHarm_;    

    //! Pointer to object of analysis (Static, Trans, Harm or Eig)
    Assemble * assemble_;

    //! pointer to SolveStep classes
    StdSolveStep * solveStep_;
    
    BaseSystem * algsys_;      //!< pointer to algebraic system
  
    OLAS_Params * olasParams_; //!< pointer to paramter object of OLAS
    OLAS_Report * olasReport_; //!< pointer to report object of OLAS
    
    /** This is the node for linear system responsible for this pde. */
    InfoNode* olasInfo_;
    
    //! flag to check if there are initial conditions in the set up
    bool isSetInitialCondition_;
    //! initial value.
    Double InitialCondition_;
    
    //@}

    //! Map Storing FeFunctions for each unknown of Function Definitions
    std::map<SolutionType, StdVector<FunctionDescription> > functions_;

    //! map which associates a Postprocessing result to its corresponding Primary Result
    std::map<SolutionType,SolutionType> postProcResults_;

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
