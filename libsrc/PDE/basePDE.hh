#ifndef FILE_BASEPDE
#define FILE_BASEPDE

#include <Domain/elem.hh>
#include <Utils/storesol.hh>
#include <list>
#include <General/environment.hh>
#include <Domain/bcs.hh>
#include <DataInOut/timefunc.hh>
#include <DataInOut/filetype.hh>
#include <DataInOut/writeresults.hh>
#include <Matrix/matrix.hh>

#ifdef USE_OLAS
#include <olas.hh>
#else
#include <multigrid.hh>
#endif


#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/MaterialData.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Driver/assemble.hh"
#include "timestepping.hh"
#include "BaseEQN.hh"

namespace CoupledField
{
 
class SpaceErrorEstimator;

  //! Base class for partial differential equations

  //! Class BasePDE is the base class from which different types of classes
  //! describing individual types of PDEs are derived. 

  class BasePDE
  {
  public:

    friend class PDECoupling;

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

    //! read material data
    virtual void ReadMaterialData();

    //! set boundary condition
    //! \param level             level of grid
    //! \param update indicator: do we update boundary condition in algebraic
    //!                          system or set new ones?
    //! \param atimestep         time step of claculation
    virtual void SetBCs(const Integer level, const Integer update,
			const Double atimestep);

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators(const Integer level)=0;

    //! write general defines (BCs, loads, etc.) to info-file
    virtual void WriteGeneralPDEdefines();

    //! retruns the load names
    std::vector<std::string>& GetLoadDom()
    {return assemble_->loadDom_;};

    //! returns the load dofs
    std::vector<std::string>& GetLoadDof()
    {return assemble_->loadDof_;};

    //! returns the load values
    std::vector<Double>& GetLoadVals()
    {return assemble_->loadVals_;};

    //!returns the load functions
    std::vector<std::string>& GetLoadFncs()
    {return assemble_->fncname_loads_;};

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //@{
    //! \name Methods performing post-processing

    //! Do Postprocessing as descriped in conf file
    virtual void PostProcess(const Integer level) {;};

    //! write results in file
    virtual void WriteResultsInFile()=0;  
    //@}

    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================

    //! define algebraic system 
    /*! \param AS_sysid id of PDE in algebraic system  */
    virtual void SetAlgSys(int sysid);

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();

    //! define algebraic system solver Parameters
    virtual void SetSolverParameters();

    //! Init the time stepping
    //! \param dt time step
    virtual void InitTimeStepping(const Double dt)
    {Error("InitTimeStepping not implemented",__FILE__,__LINE__);};

    //! deletes the algebraic system
    void DeleteAlgSys(int as_id)
    {assemble_->DeleteAlgSys();};
  
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

  virtual void PreStepHarmonic(const Integer level) {;};

  virtual void SolveStepHarmonic(const Integer level)
  {Error("Harmonic step not implemented!",__FILE__,__LINE__);};
  virtual void StepHarmonicLin(const Integer level)
  {Error("Harmonic step not implemented!",__FILE__,__LINE__);};
  virtual void StepHarmonicNonLIn(const Integer level)
  {Error("Harmonic step not implemented!",__FILE__,__LINE__);};

  virtual void PostStepHarmonic(const Integer level) {;};

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

    //! return pointer to vector with subdomains, on which we calculate the PDE
    virtual std::vector<std::string> * getSDsPDE()
    { return &subdoms_;}

    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(std::string output)
    {
      Error("not implemented",__FILE__,__LINE__);
      return FALSE;
    }

    //! return solution
    virtual const BaseStoreSol& getS() {return *sol_;}

    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS1() const 
    { Error("Not implemented",__FILE__,__LINE__);}
    
    //! return size of solution
    virtual Integer getSize() const 
    {
      Error("Function getSize is not overloaded in this class",__FILE__,
	    __LINE__);
      return 0;
    }

    //! return number of restraints
    Integer GetNumRestraints(const Integer level=-1);

    //! return assignment array Mesh2PDENode
    inline std::vector<Integer> & GetMesh2PDENode() {return mesh2PDENode_;}

    //! return assignment array PDE2MeshNode
    inline std::vector<Integer> & GetPDE2MeshNode() {return pde2MeshNode_;}
    
    //! computes the coordinates of an element including the delta
    //! \param connect (input) global node numbers of element
    //! \param ptCoord (output) coordinates of the element nodes
    //!                (nrNodes \f$\times\f$ spaceDim);
    //! \param level (input) index for multilevel hierarchy
    virtual void GetElemCoords(const Vector< Integer > connect,
			       Matrix< Double > &coordMat,
			       const Integer level);

    //! returns the local PDE number of an array of nodes
    /*!
      \param PDENodes (output) Vector of PDE node numbers
      \param MeshNodes (input) Vector of mesh node numbers
      \param Mesh2PDENode (input) Vector assigning PDE to mesh node numbers
    */
    virtual void Mesh2PDENode(Vector<Integer> & PDENodes, 
			      const Vector<Integer> & MeshNodes,
			      const std::vector<Integer> & Mesh2PDENode);
  
    //! returns the local global Mesh node numbers of an array of nodes
    /*!
      \param MeshNodes (output) Vector of mesh node numbers    
      \param PDENodes (input) Vector of PDE node numbers
      \param PDE2MeshNode (input) Vector assigning mesh to PDE  node numbers
    */
    virtual void PDE2MeshNode(Vector<Integer> & MeshNodes, 
			      const Vector<Integer> & PDENodes,
			      const std::vector<Integer> & PDE2MeshNode);
    //@}

    // Does this method belong to postproc section?
    //! write information about relative error of calculation
    void WriteErrorInfo(WriteResults * ptmeshes);

    //! Reset function (what does it do?)
    virtual void Reset()
    { Error("Fnc is not implemented",__FILE__,__LINE__);}


  protected:
  
#ifndef XMLPARAMS
    //! check for saving parameters
    virtual void ReadSavings();
#else
    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! Currently we define a simple dummy implementation here in basePDE,
    //! but in the long-run this should become a pure virtual method.
    virtual void ReadStoreResults() {
      savesol_ = TRUE;
      savederiv1_ = TRUE;
      savederiv2_ = TRUE;
    };
#endif

    //! read from config-file info about BCs
    void ReadBCs(const std::string eq);
    
    //! return index of dof defined by keyword (e.g. 'ux')
    virtual Integer GetBCDof(const std::string keyword)
    {
      Error("GetBCDof not implemented",__FILE__,__LINE__);
      return 0;
    }


    //! assign local PDE node numbers to own subdomains
    /*!
      \param mesh2PDENode (output) Vector assigning mesh to PDE node numbers
      \param pde2MeshNode (output) Vector assigning PDE to mesh node numbers
      \param subdoms (input) Vector of subdomains which are to be mapped
    */
    virtual void AssignPDENodeNumbers(std::vector<Integer> & mesh2PDENode,
				      std::vector<Integer> & pde2MeshNode,
				      const std::vector<std::string> & subdoms);


    //! assign local PDE element numbers to own subdomains
    /*!
      \param mesh2PDENode (output) Vector assigning mesh to PDE element numbers
      \param pde2MeshNode (output) Vector assigning PDE to mesh element numbers
      \param subdoms (input) Vector of subdomains which are to be mapped
    */
    //! \todo Elena: Is there a Mapping of Elements in adapted grids?
    virtual void AssignPDEElemNumbers(std::vector<Integer> & mesh2PDEElem,
				      std::vector<Integer> & pde2MeshElem,
				      const std::vector<std::string> & subdoms);
    

#ifdef ADAPTGRID  
    //! ----------------- functions for adaptivity
    //! in this fnc we delete old pointer to Error-object & create new one
    void ConstructorError();
#endif


    /// returns the time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolOfElement(Matrix<Double>& sol, Vector<Integer>& connect_PDE);

    /// returns the vector of time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, Vector<Integer>& connect_PDE);

    /// calc the normal vector of a line element (for acoustic coupling)
    void CalcLineNormalVec(Vector<Double>& n, Matrix<Double>& ptCoord);
  
    /// calc the normal vector of an interface element (outwards of the pde subdobmain
    void CalcLineNormalVec(Vector<Double>& n, const Elem& interfaceElem, const Elem& neighbour);


    // ======================================================
    // DATA SECTION
    // ======================================================

    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes related to geometry and node numbering
    Integer dofspernode_;  //!< number of unknowns per node
    Integer dofsperedge_;  //!< number of unknowns per edge
    Integer numPDENodes_;  //!< number of nodes in subdomains
    Integer numElems_;     //!< number of elements in subdomains 
    Grid * ptgrid_;        //!< pointer to grid object

    //! subdomain-levels belonging to PDE
    std::vector<std::string> subdoms_;

    //! surface-domain-levels belongig to PDE
    std::vector<std::string> surfdoms_;

    // Assignment MeshNodeNumers <-> PDENodeNumbers
    // Note: PDE/Mesh-Node numbers start with 1, but the arrayindex always
    // starts with 0 so correct transformation reads as follows:
    // PDENode = Mesh2PDENode[PDENode - 1]

    //! vector containing PDE (=local) node numbers
    std::vector<Integer> mesh2PDENode_;

    //! vector containing Mesh (=global) node numbers
    std::vector<Integer> pde2MeshNode_;
   
    //! vector containing PDE (=local) element numbers
    std::vector<Integer> mesh2PDEElem_;

    //! vector containing Mesh (=global) element numbers
    std::vector<Integer> pde2MeshElem_;
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
    std::vector<Elem*> couplingElements;

    //! iteration counter for coupled PDE solution process
    Integer iterCoupledCounter_;
    //@}

    // -----------------------------------------------------------------------
    // Time stepping
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to time stepping
    Double stepTime_;             //!< time step;
    TimeStepping * TS_alg_;       //!< handles the time stepping
    Boolean effectiveMass_;       //!< use effective mass formulation for transient analysis
    Boolean firstTimeStepStatic_; //!< needed for coupled, iterative methods
    //@}

    // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to the handling of boundary conditions

    //! names of interfaces with homogeneous Dirichlet boundary conditions
    std::vector<std::string> bcs_hd_;

    //! names of interfaces with inhomogeneous Dirichlet boundary conditions
    std::vector<std::string> bcs_id_;

    //! names of surfaces with homogeneous von Neumann boundary conditions
    std::vector<std::string> bcs_nh_;

    //! names of surfaces with inhomogeneous von Neumann boundary conditions
    std::vector<std::string> bcs_ni_;

    //! names of surfaces with inhomogeneous Robin boundary conditions
    std::vector<std::string> bcs_rh_;

    //! names of surfaces with inhomogeneous Robin boundary conditions
    std::vector<std::string> bcs_ri_;

    //! values of solution at inhomogeneous Dirichlet boundaries
    std::vector<Double> val_id_;

    //! set, if BCs already set (shouldn't this better be a bool?)
    Integer updateBCs_;

    //! names of the time functions for inhomogeneous Dirichlet BCs
    std::vector<std::string> fncnames_id_;

    //! degrees of freedom (e.g. ux) of homogenous Dirichlet BC

    //! This is a vector of strings, which describe the degrees of freedom
    //! (e.g. ux) that are fixed at a certain interface by the homogenous
    //! Dirichlet boundary conditions.
    std::vector<std::string> homDirichDof_;

    //! This is a vector of strings, which describe the degrees of freedom
    //! (e.g. ux) that are fixed at a certain interface by the inhomogenous
    //! Dirichlet boundary conditions.
    std::vector<std::string> inhomDirichDof_;

    //! pointer to boundary condition object
    BCs *ptBCs_;
    //@}

    // -----------------------------------------------------------------------
    // Storing information
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to storing information
    
    //! TRUE, if solution should be written to result file
    Boolean savesol_;

    //! TRUE, if first derivative of solution should be written to result file
    Boolean savederiv1_;

    //! TRUE, if second derivative of solution should be written to result file
    Boolean savederiv2_;

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
    Vector<Double> markingElems_;

    Vector<Double> errorMap_;  //!< array with error map
    Double tolSpaceErr_;       //!< tolerance
    Boolean recalc_;           //!< ???
    //@}

    // -----------------------------------------------------------------------
    // Solver parameters
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to parameters for solver
    Integer maxnumiter_;     //!< maximum of iterations (for iterative solver)
#ifdef USE_OLAS
  SolverType solvertype_;    //!< type of solver (see environment.hh)
  PrecondType precondtype_;  //!< type of preconditioner (see environment.hh)
#else
    Integer solvertype_;     //!< type of solver (see las_environment.hh)
    Integer precondtype_;    //!< type of preconditioner (see las_environment.hh)
#endif
    Integer numeqcoarse_;    //!< number of unknowns on coarse level(just for AMG)
    Double  eps_;            //!< accuracy
    Double dampiter_;        //!< damping parameter within iterative solution
    Double coarsealpha_;     //!< coarsening factor (just for AMG)
    Double lasttimecalc_;    //!< Last time on which we have calculated solution

    //! Number of last timestep on which we have calculated our solution
    Integer laststepcalc_;

    //! specifies the type of damping model (see environment.hh)
    DampingType damping_type_;
    //@}

    // -----------------------------------------------------------------------
    // Non-linearity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to nonlinearity
    Boolean nonLin_;           //!<  flag for nonlinear calculations
    Boolean geoUpdate_;        //!< flag for geometric update
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    Boolean lineSearch_;       //!< switch for lineSearch

    //@}

    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    AnalysisType analysistype_; //!< analysis type
    std::string pdename_;       //!< type of PDE (set in the derived classes)
    ShortInt dim_;              //!< space dimension of pde
    Boolean isaxi_;             //!< TRUE: axisymmetric problem
    Boolean isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    BaseEQN * EqnData_;         //!< equation handling

      BaseStoreSol * sol_;      //!< solution

    Boolean initMatrices_; //!< true, if matrix is set up each iteration step

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

