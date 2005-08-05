#ifndef FILE_STDPDE
#define FILE_STDPDE
#include <fstream>
#include "PDE/basePDE.hh"

#include <set>

#include "PDE/timestepping.hh"
#include "nodeEQN.hh"

namespace CoupledField {


  // forward class declarations
  class PDECoupling;
  class WriteResults;
  class TimeStepping;
  class TimeFunc;
  class BaseNodeStoreSol;
  class NodeEQN;
  class StdSolveStep;
  class PDECoupling;
  
  //! Base class for all single-field and direct-coupled problems

  class StdPDE : public BasePDE {
  
  public:

    // friend cStdlass declarations
    friend class StdSolveStep;
    friend class PDEMemento;
    friend class PDECoupling;

    //! Virtual destructor
    virtual ~StdPDE();
    
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! Define the algebraic system 
    virtual void DefineAlgSys() = 0;

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();
  
    //! Deletes the algebraic system
    void DeleteAlgSys() {
      ENTER_FCN( "StdPDE::DeleteAlgsys" );
      if ( algsys_ != NULL ) {
        delete algsys_;
        algsys_ = NULL;
      }
    };

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
    //! get the encapsulated state of the PDE
  
    //! returns the current state of the PDE (solution, derivative,
    //! coupling-objects) in an encapsulated object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the solution, the derivative and perhaps coupling 
    //! values like geometry update have to be passed. 
    //! The PDEMemento object encapsulates this information. 
    //! Later on the information can be given back to the PDE
    //! with the method SetMemento();
    //! \param memento (output) Object where the current state gets saved
    virtual void GetMemento(PDEMemento & memento);
  
    //! set the encapsulated state of the PDE
  
    //! set the current state of this PDE (solution, derivative,
    //! coupling-objects) from an encapsulated object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the solution, the derivative and perhaps coupling 
    //! values like geometry update have to be passed. 
    //! The PDEMemento object encapsulates this information. 
    //! With this method the previous stored information can be set
    //! to the current PDE.
    //! \param memento (input) Previously saved state of the PDE
    //! \param frequency   (input) : frequency of previous sequence
    virtual void SetMemento(PDEMemento & memento, std::string transFromTo,
                            Double frequency);
                   
    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();
    
    //! return number of restraints
    virtual UInt GetNumRestraints( ) = 0;

    virtual UInt GetTimeStepCounter();
   

    //! Get types of needed matrices (sysmtem, stiffness,..)
    virtual void GetMatrixTypes( std::set<FEMatrixType> &matTypes) = 0;

    //! return pointer to vector with subdomains, on which we calculate the PDE
    virtual StdVector<RegionIdType> * getSDsPDE()
    { return &subdoms_;}

    //! Get type of analysis
    AnalysisType GetAnalysisType() {
      return analysistype_;
    }
  
    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output)
    {
      Error("not implemented",__FILE__,__LINE__);
      return FALSE;
    }
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS1() const;
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS2() const;
  
    //! return size of solution
    virtual UInt getSize() const 
    {
      Error("Function getSize is not overloaded in this class",__FILE__,
            __LINE__);
      return 0;

    }

    //! Get coefficient for damping matrix in fractional damping model
    //! \todo This function has to be removed when the fractional
    //! damping model gets implemented in a separate Forms-class
    virtual Double GetFracDampMatrixCoeff(UInt actSD);

	//! Also for fractional damping model do obtain
	virtual UInt GetFracMemory() {
	  return fracMemory_;
	}
	virtual Boolean GetIsaxi() {
	  return isaxi_;
	}

    //! computes the coordinates of an element including the delta
    //! \param connect (input) global node numbers of element
    //! \param ptCoord (output) coordinates of the element nodes
    //!                (nrNodes \f$\times\f$ spaceDim);
    virtual void GetElemCoords(const StdVector< UInt > connect,
                               Matrix< Double > &coordMat );
  
    //!
    //! \todo What is the purpose of this function?
    virtual void ComputeRHS(const Double atime) {;};
  
    //! set boundary condition
    //! \param atimestep         time step of claculation
    virtual  void SetBCs( const Double atimestep ) = 0;
  
  
    //! transform solution and derivatives due to slicing technique
    //void TransformSol4Slice(UInt & nodeShift, UInt & shiftFactor, 
    //		    const UInt flag);
    void TransformSol4Slice(UInt & shiftFactor, UInt & nodeShift,
				    UInt & elemgrid, Double &  meshsize, const UInt flag);

    //! save solution of special nodes
    void SaveNodes(const UInt shiftFactor, const Double timeStep,
		   const UInt numShift, const Integer nodeShift, 
		   const UInt maxnumelemz_);


    //@{
    //! store the new solution returned by the algebraic system
    //! \param ptSol pointer to solution array
    //! \param size legnth of solution array
    virtual void SaveSolution( const Double * ptSol, UInt size ) = 0;
    virtual void SaveSolution( const Complex * ptSol, UInt size ) = 0;
    //@}

    //! get the data vector of the current solution of a PDE.
    CFSVector * GetSolutionVector();
    
    /// returns the time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolOfElement(Matrix<Double>& sol, StdVector<UInt>& connect_PDE);
    
    /// returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);
    
    /// returns the vector of time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);
    
    /// returns the vector of 2nd time derivative of the solution belonging to all nodes 
    /// of the actual element
    void GetDeriv2SolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);
    
    
    // ======================================================
    // METHODS FOR ASSEMBLING
    // ======================================================

    //! specify type of system matrix for AlgebraicSystem
    virtual void AssembleMatrices( ) = 0;
    
    //! setup source term
    virtual void AssembleSrcRHS( const Double time = 0.0 ) = 0;
    
    //!  assemble a nonlinear RHS part
    virtual void AssembleNLRHS( const Double time = 0.0) = 0;

    //!  assemble a spring into the system matrix
    virtual void AssembleSprings( const Double time = 0.0) = 0;

    //! Initialize all matrices with nonlinear behavior
    virtual void InitNonLinMatrices() = 0;

    //! constructes the matrix graph by providing to the algebraic 
    //! system the element connectivities
    virtual void SetupMatrixGraph() = 0;

    //! trigger the reassmbling of the matrices
    virtual void SetReassemble() = 0;

    //! sets the actual frequency (just needed for harmonic analysis)
    virtual void SetFrequency(Double actFreq) = 0;

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

    MaterialData * getPDEMaterialData()
    {return materialData_;};
    
    BaseNodeStoreSol * getPDESolution() {return sol_;};
    
    BaseSystem * getPDE_algsys(){return algsys_;};
  
    UInt getPDE_numElems(){return numElems_;};

    UInt getPDE_dofspernode(){return dofspernode_;};
  
    UInt getPDE_numPDENodes(){return numPDENodes_;};
  
    UInt getPDE_spaceDim(){return dim_;};
  
    NodeEQN * getPDE_eqnData(){return eqnData_;}; 
  
    Grid * getPDE_grid(){return ptgrid_;};  
  
    Assemble * getPDE_assemble(){return assemble_;}
  
    StdVector<RegionIdType> getPDE_subdoms(){return subdoms_;}
   
    Vector<Complex> getPDE_complexValuedCharge()
    {return complexValuedCharge_;};

    virtual void setBCs_id_phase_(UInt i, Double & phase);
    
    void setPDE_piezoMaterialType(piezoMaterialType & pMatType){
      piezoMaterialType_ = pMatType;};
  
    piezoMaterialType getPDE_piezoMaterialType()
    {return piezoMaterialType_;}
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

    //!
    void sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted);
    void sortStresses(Vector<Complex>& unsorted, Vector<Complex>& sorted);

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
  
    StdVector<RegionIdType> pressSurf_;  //!< surface of pressure loads
    StdVector<Double>      pressVals_;  //!< values of the pressure loads
    StdVector<Double>      pressPhase_;  //!< phase of the pressure loads 
    //!(in case of harmonic analysis)
    StdVector<std::string> pressFnc_;   //!< function names of pressure loads

    //@}

  
    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering
    EQNType eqnType_;      //!< type of equation numbering used
    UInt dofspernode_;  //!< number of unknowns per node
    UInt dofsperedge_;  //!< number of unknowns per edge
    UInt numPDENodes_;  //!< number of nodes in subdomains
    UInt numElems_;     //!< number of elements in subdomains
  
    //! defines subtype of mechanic PDE: plainStrain, 3d, ...
    std::string subType_;
    //@}

    // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to the handling of boundary conditions
  
    //! names of interfaces with homogeneous Dirichlet boundary conditions
    StdVector<std::string> bcs_hd_;
  
    //! names of interfaces with inhomogeneous Dirichlet boundary conditions
    StdVector<std::string> bcs_id_;
  
    //! names of surfaces with homogeneous von Neumann boundary conditions
    StdVector<std::string> bcs_nh_;
  
    //! names of surfaces with inhomogeneous von Neumann boundary conditions
    StdVector<std::string> bcs_ni_;
  
    //! names of surfaces with inhomogeneous Robin boundary conditions
    StdVector<std::string> bcs_rh_;
  
    //! names of surfaces with inhomogeneous Robin boundary conditions
    StdVector<std::string> bcs_ri_;
  
    //! values of solution/magnitude at inhomogeneous Dirichlet boundaries
    StdVector<Double> val_id_;
  
    //! values of phases at inhomogeneous Dirichlet boundaries
    StdVector<Double> bcs_id_phase_;
  
    //! names of the time functions for inhomogeneous Dirichlet BCs
    StdVector<std::string> fncnames_id_;
  
    //! degrees of freedom (e.g. ux) of homogenous Dirichlet BC
  
    //! This is a vector of strings, which describe the degrees of freedom
    //! (e.g. ux) that are fixed at a certain interface by the homogenous
    //! Dirichlet boundary conditions.
    StdVector<std::string> homDirichDof_;
  
    //! This is a vector of strings, which describe the degrees of freedom
    //! (e.g. ux) that are fixed at a certain interface by the inhomogenous
    //! Dirichlet boundary conditions.
    StdVector<std::string> inhomDirichDof_;

    //! number of dirichlet boundary conditions
    UInt numDirichletBCs_;

    //! number of built-in dirichlet conditions
    UInt numBuildInDirichletBCs_;

    //@}

    // -----------------------------------------------------------------------
    // Non-linearity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to nonlinearity
    Boolean nonLin_;           //!< flag for nonlinear calculations
    Boolean geoUpdate_;        //!< flag for geometric update
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    Boolean nonLinLogging_;    //!< log progress of non-linear iterations
    Boolean isHysteresis_;     //!< flag for hysteresis

    std::string lineSearch_;   //!< switch for lineSearch
    StdVector<NonLinPDE> nonLinPDEName_;//!< some PDEs carry a name (->acoustics!)
    //@}

    // -----------------------------------------------------------------------
    // Material data
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes handling info on material data
    // Char * charMaterialFileName_;    //!< name of material file
    // LoadMaterialData *loadMaterial_; //!< material reader
    MaterialData *materialData_;     //!< material data structure
    std::string pdematerialclass_;    //!< material class
  
    //! Data Type which decides wheather material is real or complex
    piezoMaterialType piezoMaterialType_;
    //! contains element results of complex valued charge 
    Vector<Complex> complexValuedCharge_;
    Vector<Complex> complexValuedEfield_;
    //@}


    // -----------------------------------------------------------------------
    // PDE coupling
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to handling PDE coupling
    Boolean isIterCoupled_;        //!< PDE couples with others
    Matrix<Double> deltCoords_;    //!< offset to grid coordinates
    Vector<Double> matParam_;      //!< change to material parameter
    Boolean updateCouplingBCs_ ;  //!< flag if coupling BC were already set
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
    Boolean effectiveMass_;       //!< use effective mass formulation for transient analysis
    Boolean firstTimeStepStatic_; //!< needed for coupled, iterative methods

    //@}


    WriteResults * outFile_;  //!< pointer to output file


    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    AnalysisType analysistype_; //!< analysis type
    Boolean isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    UInt dim_;                  //!< space dimension of pde
    Boolean isaxi_;             //!< TRUE: axisymmetric problem
    Boolean isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    NodeEQN * eqnData_;         //!< equation handling

    BaseNodeStoreSol * sol_;    //!< solution
    CFSVector * solVec_;        //! needed in iterative coupled computation 

    //! specifies the type of damping model (see environment.hh)
    DampingType dampingType_;
	//! list of damping types for all regions
    StdVector<DampingType> dampingList_;
	//! number of old time steps to be saved (for fractional damping)
    UInt fracMemory_;
	//! type of interpolation (for fractional damping)
    InterpolType inType_;
  
	//! flag indicating, if damping matrix is needed
	Boolean needsDampingMatrix_;
	
	//! checks, if we have for the coupling a incremental solution
    Boolean isIncrFormulation_;    

    //! set defining which type of matrices (stiffness, mass,...) is used
    std::set<FEMatrixType> matrixTypes_;

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
