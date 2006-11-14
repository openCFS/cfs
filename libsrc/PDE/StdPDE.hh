#ifndef FILE_STDPDE
#define FILE_STDPDE
#include <fstream>
#include "PDE/basePDE.hh"

#include <set>

#include "PDE/timestepping.hh"
#include "Domain/Composite.hh"

namespace CoupledField {


  // forward class declarations
  class PDECoupling;
  class WriteResults;
  class TimeStepping;
  class TimeFunc;
  class BaseNodeStoreSol;
  class StdSolveStep;
  class PDECoupling;
  class EqnMap;
  
  //! Base class for all single-field and direct-coupled problems

  class StdPDE : public BasePDE {
  
  public:

    // friend cStdlass declarations
    friend class PDECoupling;

    // public typedefs
    typedef StdVector<shared_ptr<ResultDof> > ResultList;


    //! Virtual destructor
    virtual ~StdPDE();
    
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! Define the algebraic system 
    virtual void DefineAlgSys() = 0;

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();

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
    // POSTPROCESSING METHODS
    // ======================================================
    
    //! compute volume above a deformed surface
    virtual void ComputeVolDefSurf(StdVector<RegionIdType> &surfRegions, 
                                   StdVector<std::string> &strDir);
    
    //!
    virtual Double ComputeVolElem(BaseFE * ptSurfEl, Matrix<Double>& SurfCoord, 
                                  Vector<Double> disp);
    
    //!
    virtual Complex ComputeVolElem(BaseFE * ptSurfEl, Matrix<Double>& SurfCoord, 
                                   Vector<Complex> disp);

    // ======================================================
    // GET/SET METHODS
    // ======================================================

    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();

    //! Return vector with result types
    ResultList& GetResults() { return results_;}

    //! return pointer to vector with subdomains, on which we calculate the PDE
    virtual StdVector<RegionIdType> * getSDsPDE()
    { return &subdoms_;}

    virtual AnalysisType GetAnalysisType() {
      return analysistype_;
    }
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output)
    {
      Error("not implemented",__FILE__,__LINE__);
      return false;
    }
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS1() const;
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS2() const;

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

    //!
    //! \for computing and adding RHS to PDE in case of special sources 
    virtual void ComputeRHS(const Double atime) {;};

    //!
    //! \for computing vortex source both analytically and with complex 
    //! \potential function 
    virtual  void VortexAnalytical(Double & press, Vector<Double>& dTij_di, const Double x,
                                   const Double y, const Double t, 
                                   const UInt outType){
      Error("VortexAnalytical is only implemented in acouFlowNoisePDE",__FILE__,__LINE__);};
  
    //! set boundary condition
    //! \param atimestep         time step of claculation
    virtual  void SetBCs( const Double atimestep ) = 0;


    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    virtual void SaveSolution( const Double * ptSol, UInt size ) = 0;
    virtual void SaveSolution( const Complex * ptSol, UInt size ) = 0;
    //@}

    //! get the data vector of the current solution of a PDE.
    CFSVector * GetSolutionVector();
    
    /// returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement(Vector<Double>& sol, const EntityIterator& it);
    void GetSolVecOfElement(Vector<Complex>& sol, const EntityIterator& it);
    
    /// returns the vector of time derivative of the solution belonging 
    /// to all nodes of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, const EntityIterator& it);
    void GetDerivSolVecOfElement(Vector<Complex>& sol, const EntityIterator& it);
    
    /// returns the vector of 2nd time derivative of the solution belonging to all nodes 
    /// of the actual element
    void GetDeriv2SolVecOfElement(Vector<Double>& sol, const EntityIterator& it);
    void GetDeriv2SolVecOfElement(Vector<Complex>& sol, const EntityIterator& it);

    //! Init the time stepping
    virtual void InitTimeStepping()
    {Error("InitTimeStepping not implemented",__FILE__,__LINE__);};

    // ======================================================
    // COMMUNICATION ROUTINES FOR PARAMETER IDENTIFICATION
    // ======================================================

    //@{
    //!  The following methods are used only durig parameter
    //!  identification process! Maybe one day a more to CFS++ consistent 
    //!  nomenclature would be nice ...

    shared_ptr<EqnMap> GetEqnMap() { return eqnMap_; }

    std::map<RegionIdType, BaseMaterial*>  getPDEMaterialData()
    {return materials_;};
    
    BaseNodeStoreSol * getPDESolution() {return sol_;};
    
    BaseSystem * getPDE_algsys(){return algsys_;};
  
    UInt getPDE_numElems(){return numElems_;};

    UInt getPDE_numPDENodes(){return numPDENodes_;};
  
    UInt getPDE_spaceDim(){return dim_;};
  
    Grid * getPDE_grid(){return ptgrid_;};  
  
    Assemble * getPDE_assemble(){return assemble_;}
  
    StdVector<RegionIdType> getPDE_subdoms(){return subdoms_;}
   
    Vector<Complex> getPDE_complexValuedCharge()
    {return complexValuedCharge_;};
    
    void setPDE_MatDataType(DataType & pMatType){
      matDataType_ = pMatType;};
  
    DataType getPDE_MatDataType()
    {return matDataType_;}

    WriteResults * getPDE_outFile(){return outFile_;};

    //! get pointer to time stepping object
    TimeStepping * getTimeStepping() { return TS_alg_;};

    void setPDE_complexValuedCharge(Vector<Complex> chargeVec)
    {complexValuedCharge_=chargeVec;};

    //!
    void sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted);
    void sortStresses(Vector<Complex>& unsorted, Vector<Complex>& sorted);

    //@}

    //@{
    //!  Get functions concerning nonlinearity

    bool IsNonLin() 
    { return nonLin_;};

    bool GetNonlinLogging() 
    { return nonLinLogging_;};

    bool IsHysteresis() 
    { return isHysteresis_;};

    bool IsIterCoupled() 
    { return isIterCoupled_;};

    bool& IsFirstTimeStepStatic()
    { return firstTimeStepStatic_;};

    bool IsComputeRHS4HarmSet()
    { return ComputeRHSforHarm_;};

    Double  GetIncStopCrit() 
    { return incStopCrit_;};

    Double  GetResidualStopCrit() 
    { return residualStopCrit_;};

    Double GetRhsL2Norm(Vector<Double>& actRHS) 
    { return RhsL2Norm(actRHS);};

    UInt    GetNonlinMaxIter() 
    { return nonLinMaxIter_;};

    std::string GetLineSearch() 
    { return lineSearch_;};

    std::string GetNonlinMethod()  
    { return nonLinMethod_;};

    StdVector<NonLinPDE>& GetNonlinPDEName() 
    { return nonLinPDEName_;};

    UInt& GetIterCoupledCounter() 
    { return iterCoupledCounter_;};

    //@}

  protected:
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    StdPDE(Grid *aptgrid, WriteResults *aOutFile, TimeFunc *aTimeFunc);
  
    //! private copy constructor
    StdPDE & operator= (const StdPDE & myPDE) {
      Error ("Not implemented", __FILE__, __LINE__);

      // The following line is only to satisfy the compiler
      return *this;
    };
  


    //! stores an algsys_ vector into a StdVector
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    //! calculates L2-norm of RHS regarding entries due to penalty formulation
    Double RhsL2Norm(Vector<Double>& stdVec);

    //! Get coefficient for damping matrix in fractional damping model
    //! \todo This function has to be removed when the fractional
    //! damping model gets implemented in a separate Forms-class
    virtual Double GetFracDampMatrixCoeff(RegionIdType region);


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
  
    //! surface-domain-levels belongig to PDE
    StdVector<RegionIdType> surfdoms_;

    //@}

  
    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering
    UInt numPDENodes_;  //!< number of nodes in subdomains
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
    
    //! List of inhomogeneous Nuemann boundary conditions
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
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    bool nonLinLogging_;    //!< log progress of non-linear iterations
    bool isHysteresis_;     //!< flag for hysteresis

    std::string lineSearch_;   //!< switch for lineSearch
    StdVector<NonLinPDE> nonLinPDEName_;//!< some PDEs carry a name (->acoustics!)
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
    DataType  matDataType_;
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

    //@}


    WriteResults * outFile_;  //!< pointer to output file


    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    //! Vector containing the results calculated by this PDE
    ResultList results_;

    //! flag indicating if this PDE needs the algebraic system
    bool needsAlgsys_;

    AnalysisType analysistype_; //!< analysis type
    bool isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    UInt dim_;                  //!< space dimension of pde
    bool isaxi_;             //!< true: axisymmetric problem
    bool isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    shared_ptr<EqnMap> eqnMap_; //!< new equation handling

    BaseNodeStoreSol * sol_;    //!< solution
    CFSVector * solVec_;        //! needed in iterative coupled computation 

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
    TimeFunc * ptTimeFunc_;    //!< pointer to time functions

  
    OLAS_Params * olasParams_; //!< pointer to paramter object of OLAS
    OLAS_Report * olasReport_; //!< pointer to report object of OLAS
    //@}

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
