#ifndef FILE_BASEPDE
#define FILE_BASEPDE

#include <list>

#include "Utils/StdVector.hh"
#include "Domain/elem.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "General/environment.hh"
#include "Domain/bcs.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/filetype.hh"
#include "DataInOut/writeresults.hh"
#include "Matrix/matrix.hh"

#ifdef USE_OLAS
#include <olas.hh>
#else
#include <multigrid.hh>
#endif 


#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/LoadMaterialDataFile.hh"
#include "DataInOut/MaterialData.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Driver/assemble.hh"
#include "timestepping.hh"
#include "nodeEQN.hh"
#include "pdememento.hh"

#ifdef USE_DATABASE
#include "DataInOut/LoadMaterialDataDatabase.hh"
#endif

namespace CoupledField
{


  class SpaceErrorEstimator;

  //! Base class for partial differential equations

  //! Class BasePDE is the base class from which different types of classes
  //! describing individual types of PDEs are derived. 

  class BasePDE
  {
  public:

    // friend class declarations
    friend class PDECoupling;
    friend class PDEMemento;

    //! Constructor
    /*!
      \param aptgrid pointer to grid
      \param aptBCs pointer to boundary condition object
      \param aInFile pointer to class FileType. input data.
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    BasePDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile,
	    WriteResults *aOutFile, TimeFunc *aTimeFunc);

    //! Destructor
    virtual ~BasePDE();

    //! Initializes PDE 

    //! Initializes the PDE. This function is only called one time.
    //! \param bcSequenceId (input) name of the tag for current set of 
    //! boundary condition
    virtual void Init(Integer sequenceStep = 0,
		      std::string  bcSequenceTag = "anyTag");
    

    //! read material data
    //    virtual void ReadMaterialData();
    
    //virtual void ReadSpecialBCs(){}
  
    // The following Methods are used only durig parameter Identification process!!

    MaterialData * getPDEMaterialData(){return materialData_;};
    BaseNodeStoreSol * getPDESolution(){return sol_;};
    BCs * getPDE_BCs(){return ptBCs_;};
    BaseSystem * getPDE_algsys(){return algsys_;};
    Integer getPDE_numElems(){return numElems_;};
    Integer getPDE_dofspernode(){return dofspernode_;};

    Integer getPDE_numPDENodes(){return numPDENodes_;};
    Integer getPDE_spaceDim(){return dim_;};
    NodeEQN * getPDE_eqnData(){return eqnData_;}; 
    Grid * getPDE_grid(){return ptgrid_;};  
    Assemble * getPDE_assemble(){return assemble_;}
    StdVector<std::string> getPDE_subdoms(){return subdoms_;}
    Boolean BooleanComplexMaterialData_;
    Boolean converged_; //!< needed for coupling with MpCCI
   
    void setBCs_id_phase_(Integer i, Double & phase){
      ENTER_FCN("basePDE::setBCs_id_phase");
      if (bcs_id_phase_.GetSize() <= i)
	Error("no such index in Vector bcs_id_phase_",__FILE__,__LINE__);
      bcs_id_phase_[i]=phase;
    };
    void setPDE_actFrequency(Double & freq){
      ENTER_FCN("basePDE::setPDE_actFrequency");
      actFrequency_ = freq;
    };
    void setPDE_actFreqStep(Integer & fstep){
      ENTER_FCN("basePDE::setPDE_actFreqStep");
      actFreqStep_ = fstep;
    };

    Vector<Complex> complexValuedCharge_;
    Vector<Complex> getPDE_complexValuedCharge(){return complexValuedCharge_;};

    piezoMaterialType piezoMaterialType_; 

    void setPDE_piezoMaterialType(piezoMaterialType & pMatType){
      piezoMaterialType_ = pMatType;};

    piezoMaterialType getPDE_piezoMaterialType(){return piezoMaterialType_;}
        
    
    //	StdVector< StdVector<BaseIntDescriptor *>* > * rhsSrcIntegrators_;
    
    //StdVector< StdVector<BaseIntDescriptor *>* > * getPDE_rhsSrcIntegrators(){return rhsSrcIntegrators_;};
    //NodeStoreSol * getPDESolution(){return sol_;};
    
    //! set boundary condition
    //! \param level             level of grid
    //! \param atimestep         time step of claculation
    virtual void SetBCs(const Integer level, 
			const Double atimestep);

    //! write general defines (BCs, loads, etc.) to info-file
    virtual void WriteGeneralPDEdefines();

    //! returns the load names
    StdVector<std::string>& GetLoadDom()
    {return assemble_->loadDom_;};

    //! returns the load dofs
    StdVector<std::string>& GetLoadDof()
    {return assemble_->loadDof_;};

    //! returns the load values
    StdVector<Double>& GetLoadVals()
    {return assemble_->loadVals_;};

    //!returns the load functions
    StdVector<std::string>& GetLoadFncs()
    {return assemble_->fncname_loads_;};

    //! MpCCI gets the geometry
    virtual void PreparePDE4Computation() {;};

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //@{
    //! \name Methods performing post-processing

    //! Do Postprocessing as descriped in conf file
    virtual void PostProcess(const Integer level) {;};

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(Integer stepOffset = 0,
				    Double timeOffset = 0.0) = 0; 
    //@}

    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================

    //! define algebraic system 
    virtual void SetAlgSys();

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();

    //! set time step
    //! \params dt Current time step
    virtual void SetTimeStep(const Double dt)
    {TS_alg_->Init(matrix_factor_, dt);}

    //! deletes the algebraic system
    void DeleteAlgSys(int as_id)
    {if (algsys_)
      assemble_->DeleteAlgSys();
    };
  
    //static analysis
  
  virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset) {;};

  virtual void SolveStepStatic(const Integer kstep, const Double asteptime,
			    const Integer level, const Boolean reset);

  virtual void StepStaticLin(const Integer kstep, const Double asteptime,
			    const Integer level, const Boolean reset);

  virtual void StepStaticNonLin(const Integer kstep, const Double asteptime,
			    const Integer level, const Boolean reset)
  {Error("StepStaticNonLin not implemented!",__FILE__,__LINE__);};

  virtual void PostStepStatic(const Integer kstep, const Double asteptime,
			      const Integer level) {;};

  //harmonic analysis
  virtual void PreStepHarmonic(const Integer freqStep, const Double frequency, 
			       Integer level, const Boolean reset);

  virtual void SolveStepHarmonic(const Integer freqStep, const Double frequency, 
				 Integer level, const Boolean reset);

 virtual void StepHarmonicLin(const Integer freqStep, const Double frequency, 
				 Integer level, const Boolean reset);

  virtual void StepHarmonicNonLin(const Integer freqStep, const Double frequency, 
				  Integer level, const Boolean reset)
  {Error("Harmonic step not implemented!",__FILE__,__LINE__);};

  virtual void PostStepHarmonic(const Integer freqStep, const Double frequency, 
				Integer level, const Boolean reset) {;};

  //transient analysis
  virtual void PreStepTrans(const Integer kstep, const Double asteptime,
			    const Integer level, const Boolean reset);

  virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean updatesysmat);

  virtual void StepTransLin(const Integer kstep, const Double asteptime,
			    const Integer level, const Boolean updatesysmat);

  virtual void StepTransNonLin(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean updatesysmat)
  {Error("Nonlinear Transient Step not implemented!",__FILE__,__LINE__);};


  virtual void PostStepTrans(const Integer kstep, const Double asteptime,
			     const Integer level);

  /// initialize PDEs before iteration (done for each time step)
  void InitStepTransCoupled(Double asteptime);

  // ======================================================
  // COUPLING SECTION
  // ======================================================

  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling)
  {Error("InitCoupling Not implemented",__FILE__,__LINE__);}


  //! Fill in input coupling terms
  virtual void CalcInputCoupling();


  //! calculate coupling terms
  virtual void CalcOutputCoupling()
  {Error("InitCoupling Not implemented",__FILE__,__LINE__);}

  
    // ======================================================
    // ADATPTIVITY SECTION
    // ======================================================
#ifdef ADAPTGRID
    //@{
    //! \name Methods used for performing adaptivity

    //! test error of calculation. return TRUE, if it is more then tolerance
    virtual Boolean TestError(const Integer level);
 
    //! refine mesh
    virtual void RefineMesh(const Integer level=0);

    //@}
#endif


    // ======================================================
    // GETTER METHODS
    // ======================================================

    //@{
    //! \name Getter methods

    //! return name of pde
    virtual std::string GetName() {return pdename_;}

    //! return subtype
    virtual std::string GetSubType() {return subType_;}

    //! return pointer to vector with subdomains, on which we calculate the PDE
    virtual StdVector<std::string> * getSDsPDE()
    { return &subdoms_;}

    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output)
    {
      Error("not implemented",__FILE__,__LINE__);
      return FALSE;
    }

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
    virtual void SetMemento(PDEMemento & memento);

    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS1() const;

    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS2() const;
    
    //! return size of solution
    virtual Integer getSize() const 
    {
      Error("Function getSize is not overloaded in this class",__FILE__,
	    __LINE__);
      return 0;
    }

    //! return number of restraints
    Integer GetNumRestraints(const Integer level=-1);


    //! computes the coordinates of an element including the delta
    //! \param connect (input) global node numbers of element
    //! \param ptCoord (output) coordinates of the element nodes
    //!                (nrNodes \f$\times\f$ spaceDim);
    //! \param level (input) index for multilevel hierarchy
    virtual void GetElemCoords(const StdVector< Integer > connect,
			       Matrix< Double > &coordMat,
			       const Integer level);

    //@}

    // Does this method belong to postproc section?
    //! write information about relative error of calculation
    void WriteErrorInfo(WriteResults * ptmeshes);

    //! Reset function (what does it do?)
    virtual void Reset()
    { Error("Fnc Reset is not implemented",__FILE__,__LINE__);}

	//! Get number of time step
	virtual Integer GetTimeStepCounter()
	{ return laststepcalc_; }

	//! Get coefficient for damping matrix in fractional damping model
	virtual Double GetFracDampMatrixCoeff(Integer actSD);

  protected:

    // ======================================================
    // INITIALIZATION METHODS
    // ======================================================
  
    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    virtual void ReadStoreResults() = 0;


    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators(const Integer level)=0;

    //! read material data
    virtual void ReadMaterialData();

    //! read from config-file info about BCs
    void ReadBCs();
    
    // overloaded version of ReadBCs for special
    // boundary conditions in derived classes
    virtual void ReadSpecialBCs(){}

    //! Initialize NonLinearities
    virtual void InitNonLin(){};

    //! Init the time stepping
    virtual void InitTimeStepping()
      {Error("InitTimeStepping not implemented",__FILE__,__LINE__);};

    //! return index of dof defined by keyword (e.g. 'ux')
    virtual Integer GetBCDof(const std::string keyword);

#ifdef ADAPTGRID  
    //! ----------------- functions for adaptivity
    //! in this fnc we delete old pointer to Error-object & create new one
    void ConstructorError();
#endif


    /// returns the time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolOfElement(Matrix<Double>& sol, StdVector<Integer>& connect_PDE);

    /// returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    /// returns the vector of time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    /// returns the vector of 2nd time derivative of the solution belonging to all nodes 
	/// of the actual element
    void GetDeriv2SolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    //!
    void sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted);
    void sortStresses(Vector<Complex>& unsorted, Vector<Complex>& sorted);

    // stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);
    

    // ======================================================
    // DATA SECTION
    // ======================================================

    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes related to geometry and node numbering
    EQNType eqnType_;      //!< type of equation numbering used
    Integer dofspernode_;  //!< number of unknowns per node
    Integer dofsperedge_;  //!< number of unknowns per edge
    Integer numPDENodes_;  //!< number of nodes in subdomains
    Integer numElems_;     //!< number of elements in subdomains
    Grid * ptgrid_;        //!< pointer to grid object

    //! subdomain-levels belonging to PDE
    StdVector<std::string> subdoms_;

    //! surface-domain-levels belongig to PDE
    StdVector<std::string> surfdoms_;

    StdVector<std::string> pressSurf_;  //!< surface of pressure loads
    StdVector<Double>      pressVals_;  //!< values of the pressure loads
    StdVector<std::string> pressFnc_;   //!< function names of pressure loads

    //! defines subtype of mechanic PDE: plainStrain, 3d, ...
    std::string subType_;
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

    //@}

    // -----------------------------------------------------------------------
    // PDE coupling
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to handling PDE coupling
    Boolean pdeIsCoupled_;        //!< PDE couples with others
    Matrix<Double> deltCoords_;    //!< offset to grid coordinates
    Vector<Double> matParam_;      //!< change to material parameter
    Boolean updateCouplingBCs_ ;  //!< flag if coupling BC were already set
    Integer couplingBCsCounter_;  //!< counter for number of coupling BCs
    Integer numDirichletBCs_;     //!< number of dirichlet boundary conditions
    PDECoupling *ptCoupling_;     //!< pointer to coupling object

    //! nodes at which coupling terms are calculated
    std::list<Integer> couplingNodes;    

    //! elements at which coupling terms are calculated
    StdVector<Elem*> couplingElements;

    //! iteration counter for coupled PDE solution process
    Integer iterCoupledCounter_;
    //@}

    // -----------------------------------------------------------------------
    // Time stepping
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to time stepping
    TimeStepping * TS_alg_;       //!< handles the time stepping
    Boolean effectiveMass_;       //!< use effective mass formulation for transient analysis
    Boolean firstTimeStepStatic_; //!< needed for coupled, iterative methods
    Double lasttimecalc_;    //!< Last time on which we have calculated solution
    
    //! Number of last timestep on which we have calculated our solution
    Integer laststepcalc_;
    
    Double  actFrequency_; //!< current frequency for harmonic analysis
    Integer actFreqStep_;  //!< current frequency step for harmonic analysis
  //@}

    // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to the handling of boundary conditions

    //! tag of current set of boundary conditions

    //! tag of current set of boundary conditions. For a multiSequence-simulation
    //! this id determines, which set of boundary conditions is applied.
    std::string bcSequenceTag_;

    //! index of current set of boundary conditions. For a multiSequence-simulation
    //! this index determines, which set of boundary conditions is applied.
    Integer bcSequenceIndex_;

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

    //! pointer to boundary condition object
    BCs *ptBCs_;
    //@}

    // -----------------------------------------------------------------------
    // Storing information
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to storing information
    
    //! vector containing solutiontypes of PDE
    StdVector<SolutionType> solTypes_;

    //! vector containgin dofs of solutiontypes
    StdVector<Integer> solDofs_;

	//! TRUE, if at least one Result is calculated and written
	Boolean hasOutput_;

    //! TRUE, if solution should be written to result file
    Boolean saveSol_;

    //! TRUE, if first derivative of solution should be written to result file
    Boolean saveDeriv1_;

    //! TRUE, if second derivative of solution should be written to result file
    Boolean saveDeriv2_;

    //! TRUE, if solution should be written to history file
    Boolean saveSolHist_;

    //! TRUE, if first derivative of solution should be written to history file
    Boolean saveDeriv1Hist_;

    //! TRUE, if second derivative of solution should be written to history file
    Boolean saveDeriv2Hist_;

    WriteResults * outFile_;  //!< pointer to output file
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
    Boolean recalc_;           //!< ???
    //@}

    // -----------------------------------------------------------------------
    // Solver parameters
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to parameters for solver
    Integer maxnumiter_;       //!< maximum of iterations (for iterative solver)
    OLAS_Params * olasParams_; //! pointer to paramter object of OLAS
    OLAS_Report * olasReport_; //! pointer to report object of OLAS
    Integer numeqcoarse_;      //!< number of unknowns on coarse level(just for AMG)
    Double  eps_;              //!< accuracy
    Double dampiter_;          //!< damping parameter within iterative solution
    Double coarsealpha_;       //!< coarsening factor (just for AMG)
   
    ComplexFormat complexFormat_;  //!< outputFormat for complex numbers

    //! specifies the type of damping model (see environment.hh)
    DampingType dampingType_;
    Boolean needsDampingMatrix_;
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
    Integer nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    Boolean nonLinLogging_;    //!< log progress of non-linear iterations

    std::string lineSearch_;   //!< switch for lineSearch
    

    //@}

    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    AnalysisType analysistype_; //!< analysis type
    Boolean isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    std::string pdename_;       //!< type of PDE (set in the derived classes)
    ShortInt dim_;              //!< space dimension of pde
    Boolean isaxi_;             //!< TRUE: axisymmetric problem
    Boolean isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    NodeEQN * eqnData_;         //!< equation handling

    BaseNodeStoreSol * sol_;        //!< solution

    Vector<Double> solIncr_;      //! needed in iterative coupled computation 
    Vector<Double> actSol_;       //! needed in iterative coupled computation 
    Boolean isIncrFormulation_;   //! checks, if we have for the coupling a incremental solution
    

    //! factors for constructing effective mass / stiffness matrix
    Double matrix_factor_[4];

    //! actual level (actual = current?)
    Integer actlevel_;

    //! ???
    Integer as_sysid_;

    //! Pointer to object of analysis (Static, Trans, Harm or Eig)
    Assemble * assemble_;
    BaseSystem * algsys_;     //!< pointer to algebraic system
    FileType * inFile_;       //!< pointer to input file
    TimeFunc * ptTimeFunc_;   //!< pointer to time functions
    //@}



  };

} // end of namespace

#endif

