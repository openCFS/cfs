#ifndef FILE_BASEPDE_2001
#define FILE_BASEPDE_2001

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
#include "timestepping.hh"

namespace CoupledField
{
 
class SpaceErrorEstimator;

 //! Base class for partial differential equation
  /*! Class BasePDE is base class, from which different type of PDEs are derived. 
   */

class BasePDE
{
public:

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

  //! define algebraic system 
  /*!
    \param AS_sysid id of PDE in algebraic system
   */
  virtual void SetAlgSys(const Integer AS_sysid);

  /// Initialize all necessary matrices 
  void InitMatrices();

  //! read material data
  virtual void ReadMaterialData();

  //! set matrix factors for construction of effective mass / stiffness matrix
  virtual void SetMatrixFactors()=0;

  //! define algebraic system identifictaion
  /*!
    \param AS_sysid id of PDE in algebraic system
   */
  virtual void SetAlgSys_id(const Integer AS_sysid);

  //! define the discrete PDE
  virtual void DiscreteParamsPDE()=0;

  //! define algebraic system solver Parameters
  virtual void SetSolverParameters();

  //! constructes the matrix graph by providing to the algebraic system the element connectivities
  virtual void SetupMatrixGraph(Integer numeq, Integer graphtype);

  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling)
  {Error("Not implemented");}

  //! Create the matrices and Solver as well as Preconditioner
  virtual void CreateMatrices_Solver();

  //! Fill in input coupling terms
  virtual void CalcInputCoupling();

   //! setup element matrices for AlgebraicSystem for assembling procedure
  /*!
    \param level level of grid
   */
  virtual void SetupMatrices(const Integer level)=0;

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param matrixclass out: class of matrix (real scalar, complex scalar, etc.)
    \param matrixtype out: defines matrix types (Stiffness, Mass, etc.)
    \param graphtype out: type of graph
    \param numdofpernode out: number of dof per node
    \param numdirichlets out: number of nodes for dirichlets conditions
    \param numconstraints out: number of nodes for constraints conditions
  */
  virtual void SpecifyMatrices(Integer &matrixclass, Integer *matrixtype, Integer &graphtype, 
                               Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints){;};

  //! set boundary condition
  /*!
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atimestep time step of claculation
  */
  virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);


  //! solve one step for static problem 
  /*!
    \param level level of grid
  */
  virtual void SolveStepStatic(const Integer level)=0;

  //! solve one step for nonlinear static problem 
  /*!
    \param level level of grid
  */
  virtual void SolveStepStaticNonLin(const Integer level)
  {Error("SolveStepStaticNonLin not implemented!",__FILE__,__LINE__);};
  
  /// solve one harmonic step
  /*!
    \param level level of grid
  */
  virtual void SolveStepHarmonic(const Integer level)
  {Error("Harmonic step not implemented!",__FILE__,__LINE__);};
  


  /*!
    \param kstep number of calculating step
    \param asteptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */

  virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, 
			      const Boolean updatesysmat)=0;


  //! Init the time stepping
  /*!
    \param dt time step
  */
  virtual void InitTimeStepping(const Double dt)
    {Error("Not implemented");}

  //! Do Postprocessing as descriped in conf file
  virtual void PostProcess(const Integer level) {;};

  //! calculate coupling terms
  virtual void CalcOutputCoupling()
  {Error("Not implemented");}
  
  //! write results in file
  virtual void WriteResultsInFile()=0;  


#ifdef ADAPTGRID
  //------------------------ functions for adaptivity process ---------------------

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
  { 
    return &subdoms_;
  }

   //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output)
  {
    Error("not implemented");
    return 0;
  }

  //! return pointer to vector with solution
  virtual const Array<Double>& getS() {return sol_;}

  //! return pointer to vector with first derivative of solution
  virtual const Array<Double>& getS1() const 
  { 
    Error("Not implemented");
    return DVec;
  }

  //! return pointer to vector with first derivative of solution, calculated on previous step
  virtual const Array<Double>& getS1old() const
  { 
    Error("Not implemented");
    return DVec;
  }

  //! return pointer to vector with second derivative of solution
  virtual const Array<Double>& getS2() const 
  { 
    Error("Function getS2 is not overloaded in this class");
    return DVec;
  }

 //! return pointer to vector with second derivative of solution, calculated on previous step
  virtual const Array<Double>& getS2old() const 
  { 
    Error("Function getS2old is not overloaded in this class");
    return DVec;
  }

  //! return parameter beta from Newmark method
  virtual Double getBeta() const 
  { 
    Error("Function getBeta is not overloaded in this class");
    return Dval;
  }

  //! return parameter gamma from Newmark method
  virtual Double getGamma() const 
  { 
    Error("Function getGamma is not overloaded in this class");
    return Dval;
  }

  //! return size of solution
  virtual Integer getSize() const 
  { 
    Error("Function getSize is not overloaded in this class");
    return Ival;
  } 

  //! return id of PDE in algebraic system
  Integer GetSysId() const { return as_sysid_;} 

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
  void DeleteAlgSys(const Integer algsys_id);

  //!
  virtual void Reset()
  { Error("Fnc is not implemented",__FILE__,__LINE__);}
  

protected:
  //! generates a multi-dof-matrix with similar entries for all dofs
  virtual void MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  
			    const Integer nrDofs);

   //! read from .config-file info about BCs
   void ReadBCs(const std::string eq);

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
  virtual void AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			    std::vector<Integer> & PDE2MeshNode,
			    const std::vector<Elem*> &Elements );
  

  /// Stores result vector in the multidimensional solution array sol_
  void StoreToSolArray(Double * ptSol);
  
  /// Stores result vector in the multidimensional solution array sol_q
  void StoreVecToSolArray(std::vector<Double>& sol);

#ifdef ADAPTGRID  
  //! ----------------- functions for adaptivity
  //! in this fnc we delete old pointer to Error-object & create new one
  void ConstructorError();
#endif

  //! analysis type
  AnalysisType analysistype_;

  //! paramters for discrete PDE
  Integer MatrixType_;          //!< type of matrix (real, complex, etc.)
  Integer GraphType_;           //!< type of graph (nodal, edge,..)
  Integer SystemMatrix_;        //!< need system matrix (TRUE/FALSE)
  Integer StiffnessMatrix_;     //!< need stiffness matrix (TRUE/FALSE)
  Integer DampingMatrix_;       //!< need damping matrix (TRUE/FALSE)
  Integer MassMatrix_;          //!< need mass matrix (TRUE/FALSE)
  Integer ConvectionMatrix_;    //!< need convective matrix (TRUE/FALSE)
  Boolean InitMatrices_;        //!< true, if matrix is set up each iteration step

  std::string pdename_; //!< type of PDE (set in the derived classes)
  ShortInt Dim_;         //!< space dimension of pde
  Integer dofspernode_; //!< number of unknowns per node
  Integer dofsperedge_; //!< number of unknowns per edge
  std::vector<std::string> subdoms_;  //!< subdomain-levels belongig to PDE

  //Material data
  Char * charMaterialFileName_;    //!< name of material file
  LoadMaterialData *loadMaterial_; //!< material reader
  MaterialData *materialData_;     //!< material data structure
  std::string pdematerialclass_;    //!< material class

  Integer NumPDENodes_;  //!< number of nodes in subdomains
  Integer NumElems_;      //!< number of elements in subdomains 

  // pointers to objects
  Grid * ptgrid_;           //!< pointer to Grid
  BCs *ptBCs_;              //!< pointer to Boundary Condition  Object
  BaseSystem * algsys_;     //!< pointer to algebraic system
  FileType * InFile_;       //!< pointer tio input file
  WriteResults * OutFile_;  //!< pointer to output file
  TimeFunc * ptTimeFunc_;   //!< pointer to time functions
  PDECoupling * ptCoupling_;//!< pointer to Coupling Object

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
  Integer numDirichletBCs_;                  //!< number of dirichlet boundary conditions

  //!solver parameters
  Integer maxnumiter_;    //!< maximum of iterations (for iterative solver)
  Integer solvertype_;    //!< type of solver (see las_environment.hh)
  Integer precondtype_;   //!< type of preconditioner (see las_environment.hh)
  Integer numeqcoarse_;   //!< numbver of unknowns on coarse level (just for AMG)
  Double  eps_;           //!< accuracy
  Double dampiter_;       //!< damping parameter within iterative solution
  Double coarsealpha_;    //!< coarsening factor (just for AMG)


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
  std::vector<std::string> bcs_loads_;  //!< load BC levels


  std::vector<Double> val_id_;   //<! values of the inhomogeneous Dirichlet BC
  std::vector<Double> val_loads_; //<! values of the load BC
  Integer updateBCs_;                 //!< set, if BCs already set
  

  //! data for adaptivity
  SpaceErrorEstimator * ptError_; //!<  object with methods for error estimation
  Vector<Double> errorMap_; //!< array with error map
  Vector<Double> markingElems_; //!< array where  store number of refinement for the element
  Double tolSpaceErr_; //!< tolerance

  //Dummies, just for SUN compiler
  Array<Double> DVec;
  Double Dval;
  Integer Ival;
  Boolean Dbool;

};

} // end of namespace

#endif
