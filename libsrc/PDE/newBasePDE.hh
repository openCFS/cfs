#ifdef NEWBASEPDE

#ifndef FILE_BASEPDE
#define FILE_BASEPDE

#include <list>
#include <General/environment.hh>
#include <AlgebraicSystem/abstractAlgSys.hh>
#include <multigrid.hh>
#include <Domain/bcs.hh>
#include <DataInOut/timefunc.hh>
#include <DataInOut/filetype.hh>
#include <DataInOut/writeresults.hh>

#include <DataInOut/conffile.hh>
#include <DataInOut/LoadMaterialData.hh>
#include <DataInOut/MaterialData.hh>
#include <CoupledPDE/pdecoupling.hh>
#include <Utils/array.hh>
#include <Driver/assemble.hh>
#include "timestepping.hh"
#include "BaseEQN.hh"

namespace CoupledField
{
 
class SpaceErrorEstimator;

 //! Base class for partial differential equation
  /*! Class BasePDE is base class, from which different type of PDEs are derived. 
   */

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
  BasePDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile, WriteResults *aOutFile, 
	  TimeFunc *aTimeFunc); 

  //! Deconstructor
  virtual ~BasePDE();


  //! read material data
  virtual void ReadMaterialData();

  //! set boundary condition
  /*!
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atimestep time step of claculation
  */
 virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);


  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level)=0;
  
  //! check for saving parameters
  virtual void ReadSavings();

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

  //! Do Postprocessing as descriped in conf file
  virtual void PostProcess(const Integer level) {;};


  //! write results in file
  virtual void WriteResultsInFile()=0;  



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
  /*!
    \param dt time step
  */
  virtual void InitTimeStepping(const Double dt)
  {Error("InitTimeStepping not implemented",__FILE__,__LINE__);};

  //! deletes the algebraic system
  void DeleteAlgSys(int as_id)
  {assemble_->DeleteAlgSys();};
  
  //static analysis

  virtual void PreStepStatic(const Integer level) {;};

  virtual void SolveStepStatic(const Integer level, const Double aTime=0);
  virtual void StepStaticLin(const Integer level, const Double aTime=0);
  virtual void StepStaticNonLin(const Integer level, const Double aTime=0)
  {Error("StepStaticNonLin not implemented!",__FILE__,__LINE__);};

  virtual void PostStepStatic(const Integer level) {;};

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

  virtual void StepTransLin(const Integer level, const Boolean updatesysmat);

  virtual void StepTransNonLin(const Integer level, const Boolean updatesysmat)
  {Error("Nonlinear Transient Step not implemented!",__FILE__,__LINE__);};


  virtual void PostStepTrans(const Integer level);

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
  //! test error of calculation. return TRUE, if it is more then tolerance
  virtual Boolean TestError(const Integer level);
 
   //! refine mesh
  virtual void RefineMesh(const Integer level=0);
#endif


  //------------------------ get functions----------------------------------------

  //! return name of pde
  virtual std::string GetName() {return pdename_;}

  //! return pointer to vector with subdomains, on which we calculate the PDE
  virtual std::vector<std::string> * getSDsPDE()
  { return &subdoms_;}

   //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output)
  {Error("not implemented",__FILE__,__LINE__);}

  //! return pointer to vector with solution
  virtual const Array<Double>& getS() {return sol_;}

  //! return pointer to vector with first derivative of solution
  virtual const Array<Double>& getS1() const 
  { Error("Not implemented",__FILE__,__LINE__);}



  //! return size of solution
  virtual Integer getSize() const 
  { 
    Error("Function getSize is not overloaded in this class",__FILE__,__LINE__);
  } 


  //! return number of restraints
  Integer GetNumRestraints(const Integer level=-1);

  //! return assignment array Mesh2PDENode
  std::vector<Integer> & GetMesh2PDENode() {return Mesh2PDENode_;}

  //! return assignment array PDE2MeshNode
  std::vector<Integer> & GetPDE2MeshNode() {return PDE2MeshNode_;}
    
  //! computes the coordinates of an element including the delta
  /*!
    \param connect (input) global node numbers of element
    \param ptCoord (output) coordinates of the element nodes (nrNodes \f$\times\f$ spaceDim);
    \param level (input) index for multilevel hierarchy
  */
  virtual void GetElemCoords(const Vector< Integer > connect, Matrix< Double > &coordMat, const Integer level);
  
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
  
  //! write information about relative error of calculation
  void WriteErrorInfo(WriteResults * ptmeshes);



  //!
  virtual void Reset()
  { Error("Fnc is not implemented",__FILE__,__LINE__);}



protected:
   //! read from .config-file info about BCs
   void ReadBCs(const std::string eq);

  /// return index of dof defined by keyword (e.g. 'ux')
  virtual Integer GetBCDof(const std::string keyword)
  { Error("GetBCDof not implemented",__FILE__,__LINE__);}

  //! maps the local node solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector referring to PDE node numbers
    \param PDE2MeshNode (input) Vector assigning PDE to mesh node numbers
  */
  virtual void TransformNodeSolution(Array<Double> & MeshSol, 
				     Array<Double> & PDESol,
				     const std::vector<Integer> & PDE2MeshNode);

  //! maps the local element solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector referring to PDE node numbers
    \param Elems (input) Vector of Elements to which the PDESol belongs to (same ordering!)
     */
  virtual void TransformElemSolution(Array<Double> & MeshSol, 
				     Array<Double> & PDESol, 
				     const std::vector<Elem*> & Elems); 

  //! maps the local element solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector
    \param Elems (input) Vector of subdomains to which PDESol belongs to
     */
  virtual void TransformElemSolution(Array<Double> & MeshSol, 
				     Array<Double> & PDESol, 
				     const std::vector<std::string> & SD);

  //! maps the local node solution to the coupling nodes
  virtual void NodeSolutionToCoupling(Array<Double>& CouplingSol,
				  const std::vector<Integer>& NodeNumbers);


  //! maps the local element solution to the coupling nodes
  void ElemSolutionToCoupling(Array<Double>& CouplingSol,
			      const std::vector<Elem*>& NodeNumbers,
			      Vector<Double>& forceOnElem);
  
				     

  //! assign local PDE node numbers to subdomains
  /*!
    \param Mesh2PDENode (output) Vector assigning mesh to PDE node numbers
    \param PDE2MeshNode (output) Vector assigning PDE to mesh node numbers
    \param subdoms (input) Vector of subdomains which are to be mapped
  */
  virtual void AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			    std::vector<Integer> & PDE2MeshNode,
			    const std::vector<std::string> & subdoms);

  //! assign local PDE node numbers to a list of elements (e.g. interface) 
  /*!
    \param Mesh2PDENode (output) Vector assigning mesh to PDE node numbers
    \param PDE2MeshNode (output) Vector assigning PDE to mesh node numbers
    \param subdoms (input) Vector of elements which are to be mapped
  */
//   virtual void AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
// 			    std::vector<Integer> & PDE2MeshNode,
// 			    const std::vector<Elem*> &Elements );
  

  /// Stores result vector in the multidimensional solution array sol_
  void StoreToSolArray(Double * ptSol);
  
  /// Stores result vector in the multidimensional solution array sol_q
  void StoreVecToSolArray(std::vector<Double>& sol);

#ifdef ADAPTGRID  
  //! ----------------- functions for adaptivity
  //! in this fnc we delete old pointer to Error-object & create new one
  void ConstructorError();
#endif


  /// returns the solution vector (sol1_x, sol1_y, sol2_x, sol2_y, ..) belonging to all nodes of the actual element
  void GetSolVecOfElement(Vector<Double>& sol, Vector<Integer>& connect_PDE);

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

  //! analysis type
  AnalysisType analysistype_;

  Boolean InitMatrices_;        //!< true, if matrix is set up each iteration step

  std::string pdename_; //!< type of PDE (set in the derived classes)
  ShortInt Dim_;         //!< space dimension of pde
  Boolean isaxi_;        //!< TRUE: axisymmetric problem
  Integer dofspernode_; //!< number of unknowns per node
  Integer dofsperedge_; //!< number of unknowns per edge
  std::vector<std::string> subdoms_;  //!< subdomain-levels belongig to PDE
  std::vector<std::string> surfdoms_; //!< surface-domain-levels belongig to PDE

  //Material data
  Char * charMaterialFileName_;    //!< name of material file
  LoadMaterialData *loadMaterial_; //!< material reader
  MaterialData *materialData_;     //!< material data structure
  std::string pdematerialclass_;    //!< material class

  Integer numPDENodes_;  //!< number of nodes in subdomains
  Integer numElems_;      //!< number of elements in subdomains 

  // pointers to objects
  Grid * ptgrid_;           //!< pointer to Grid
  BCs *ptBCs_;              //!< pointer to Boundary Condition  Object
  BaseSystem * algsys_;     //!< pointer to algebraic system
  FileType * InFile_;       //!< pointer tio input file
  WriteResults * OutFile_;  //!< pointer to output file
  TimeFunc * ptTimeFunc_;   //!< pointer to time functions
  PDECoupling * ptCoupling_;//!< pointer to Coupling Object

  //Equation handling
  BaseEQN * EqnData_;

  // Solution
  Array<Double> sol_;

  // Assignment MeshNodeNumers <-> PDENodeNumbers
  // Note: PDE/Mesh-Node numbers start with 1, but the arrayindex always starts with 0
  // so correct transformation reads as follows: PDENode = Mesh2PDENode[PDENode - 1]
  std::vector<Integer> Mesh2PDENode_;  //!< array containing PDE (=local) node numbers
  std::vector<Integer> PDE2MeshNode_;  //!< array containing Mesh (=global) node numbers
  
  // coupling parameters
  Boolean PDEisCoupled_;         //!< PDE couples with others
  Array<Double> deltCoords_;    //!< offset to grid coordinates
  Array<Double> matParam_;     //!< change to material parameter
  Boolean updateCouplingBCs_ ;         //!< flag if coupling BC were already set
  std::vector<Elem*> CouplingElements; //!< elements where coupling terms are calculated
  std::list<Integer> CouplingNodes;   //!< nodes where coupling terms are calculated
  Integer couplingBCsCounter_;        //!< counter for number of coupling BCs
  Integer numDirichletBCs_;           //!< number of dirichlet boundary conditions
  Integer iterCoupledCounter_;        //!< iteration counter for coupled PDE solution process

  //time stepping stuff

  Double StepTime_; //!< time step;
  Double matrix_factor_[4]; //!< factors for constructing effective mass / stiffness matrix

  TimeStepping * TS_alg_;  //<! handels the time stepping

  Integer actlevel_; //! actual level

  Integer as_sysid_;

  //! boundary conditions
  std::vector<std::string> bcs_hd_;  //!< homogeneous Dirichlet BC levels
  std::vector<std::string> bcs_id_;  //!< inhomogeneous Dirichlet BC levels
  std::vector<std::string> bcs_nh_;  //!< homogeneous Neumann BC levels
  std::vector<std::string> bcs_ni_;  //!< inhomogeneous Neumann BC levels
  std::vector<std::string> bcs_rh_;  //!< homogeneous Robin BC levels
  std::vector<std::string> bcs_ri_;  //!< inhomogeneous Robin BC levels

  std::vector<Double> val_id_;      //!< values of the inhomogeneous Dirichlet BC
  std::vector<std::string> fncnames_id_; //!< names of the time functions for inhomogeneous Dirichlet BC 

  std::vector<std::string> homDirichDof_;   //!< dof (e.g. ux) of homogenous Dirichlet BC
  std::vector<std::string> inhomDirichDof_; //!< dof (e.g. ux) of homogenous Dirichlet BC

  Integer updateBCs_;               //!< set, if BCs already set
  

  //saving information
  Boolean savesol_;     //!< TRUE, if solution should be written to result file
  Boolean savederiv1_;  //!< TRUE, if solution should be written to result file
  Boolean savederiv2_;  //!< TRUE, if solution should be written to result file

  //! data for adaptivity
  SpaceErrorEstimator * ptError_;   //!<  object with methods for error estimation
  Vector<Double> errorMap_;         //!< array with error map
  Vector<Double> markingElems_;     //!< array where  store number of refinement for the element
  Double tolSpaceErr_;              //!< tolerance
  Boolean Recalc_;

  //!solver parameters
  Integer maxnumiter_;    //!< maximum of iterations (for iterative solver)
  Integer solvertype_;    //!< type of solver (see las_environment.hh)
  Integer precondtype_;   //!< type of preconditioner (see las_environment.hh)
  Integer numeqcoarse_;   //!< numbver of unknowns on coarse level (just for AMG)
  Double  eps_;           //!< accuracy
  Double dampiter_;       //!< damping parameter within iterative solution
  Double coarsealpha_;    //!< coarsening factor (just for AMG)
  DampingType damping_type_; //!< specifies the type of damping model (see environment.hh)

  Double lasttimecalc_;  //!< Last time on which we have calculated solution
  Integer laststepcalc_; //!< Number of last timestep on which we have calculated our solution

  Assemble * assemble_;                //!< Pointer to object of analysis (Static, Trans, Harm or Eig)
  

  // nonlinear stuff
  Boolean nonLin_; //!<  flag for nonlinear calculations
  Boolean GeoUpdate_; //!< flag for geometric update
  
  Double incStopCrit_; //!< stopping criterion for incremental error
  Double residualStopCrit_;  //!< stopping criterion for residual error

  Boolean firstTimeStepStatic_; //!< needed for coupled, iterative methods

  /// switch for lineSearch
  Boolean lineSearch_;

  /// use effective mass formulation for transient analysis
  Boolean effectiveMass_;
};

} // end of namespace

#endif

#endif // #ifdef NEWBASEPDE
