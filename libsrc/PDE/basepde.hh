#ifndef FILE_BASEPDE_2001
#define FILE_BASEPDE_2001

#include <AlgebraicSystem/abstractAlgSys.hh>
#include <multigrid.hh>
#include <Domain/bcs.hh>
#include <DataInOut/timefunc.hh>
#include <DataInOut/filetype.hh>
#include <DataInOut/writeresults.hh>
#include <DataInOut/material.hh>
#include <DataInOut/conffile.hh>
#include <DataInOut/LoadMaterialData.hh>
#include <DataInOut/MaterialData.hh>

namespace CoupledField
{

 
class TimeErrorEstimator;

template<Integer dim>
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
  virtual void BasePDE::SetupMatrixGraph(Integer numeq, Integer graphtype);

  //! Create the matrices and Solver as well as Preconditioner
  virtual void BasePDE::CreateMatrices_Solver();

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


  /*!
    \param kstep number of calculating step
    \param asteptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */

  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat)=0;


  //! Do Postprocessing as descriped in conf file
  virtual void PostProcess(const Integer level) {;};
  
  //! write results in file
  virtual void WriteResultsInFile()=0;  

  //! Calculation parameters in Newmark method
  /*!
    \param adt size of timestep
  */
  virtual void CalcParameters(const Double adt)  
  { Error("Not implemented");}

  //! save received solution as solution on the previous step 
  virtual void SaveSolAsPrevStep()
  { Error("Function RestoreSol is not implemented in this class");}  


  //------------------------ get functions--------------------------------------------------------

  //! return pointer to vector with subdomains, on which we calculate the PDE
  virtual std::vector<std::string> * getSDsPDE()
  { 
    return &subdoms_;
  }

  //! return pointer to vector with solution
  virtual const Vector<Double> & getS() const =0;

  //! return pointer to vector with first derivative of solution
  virtual const Vector<Double> & getS1() const 
  { 
    Error("Not implemented");
    return DVec;
  }

  //! return pointer to vector with first derivative of solution, calculated on previous step
  virtual const Vector<Double> & getS1old() const
  { 
    Error("Not implemented");
    return DVec;
  }

  //! return pointer to vector with second derivative of solution
  virtual const Vector<Double> & getS2() const 
  { 
    Error("Function getS2 is not overloaded in this class");
    return DVec;
  }

 //! return pointer to vector with second derivative of solution, calculated on previous step
  virtual const Vector<Double> & getS2old() const 
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


  //---------------------- part for adaptivity---------------------------------------

  //! restore solution from previous step. it is used in time adaptive procedure
  virtual void RestoreSol()
  { Error("Function RestoreSol is not implemented in this class");}  

  //! calculate energy norm
  virtual  Double CalcEnergyNorm()
  { 
    Error("Function PDE::CalcEnergyNorm is not implemented in this class");
    return Dval;
  }
  
  //! refine mesh
  virtual void RefineMesh()
  { Error("Function BasePDE::RefineMesh is not implemented in this class");} 

  //! test space error
  virtual Boolean TestError()
  { 
    Error("Function BasePDE::TestError is not implemented in this class");
    return Dbool;
  }  

  //! write additional info (marked elements, relative error) to files with mesh
  virtual void PrintMeshesInfo(WriteResults * ptMehes)
  { Error("Function BasePDE::PrintMeshesInfo is not implemented in this class");}

  //! Create pointer to according class of time error estimation
  virtual TimeErrorEstimator * CreatePtTimeError()
  { 
    Error("Function CreatePtTimeError is not overloaded in this class");
    return DtimeErrorEst;
  }  

  //! solve one step for transient problem on new mesh. it is used in adaptive procedures for space
  /*!
    \param kstep number of calculating step
    \param asteptime time of calculation
    \param level level of grid
  */
  virtual void SolveStepTransNewMesh(const Integer kstep, const Double asteptime, const Integer level)
  { 
    Error("Not implemented",__FILE__,__LINE__);
  }

  //! returns the local PDE number of an array of nodes
  /*!
    \param PDENodes (output) Vector of PDE-Numbers
    \param MeshNodes (input) Vector of mesh (=global) node numbers
  */
  virtual void Mesh2PDENode(Vector<Integer> & PDENodes, Vector<Integer> & MeshNodes);
  
  //! returns the local global Mesh node numbers of an array of nodes
  /*!
    \param MeshNodes (output) Vector of mesh (=global) node numbers    
    \param PDENodes (input) Vector of PDE-Numbers
  */
  virtual void PDE2MeshNode(Vector<Integer> & MeshNodes, Vector<Integer> & PDENodes);


protected:
  /// generates a multi-dof-matrix with similar entries for all dofs
  virtual void MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  
			    const Integer nrDofs);

   //! read from .config-file info about BCs
   void ReadBCs(const std::string eq);

  //! maps the local node solution to the global solution
  /*!
    \param MeshSol (output) Solution vector referring to Mesh node numbers    
    \param PDESol (input) Solution vector referring to PDE node numbers
  */
  virtual void TransformNodeSolution(Vector<Double> & MeshSol,  Vector<Double> & PDESol);

  //! maps the local element solution to the global solution
  /*!
    \param MeshSol (output) Solution vector referring to Mesh node numbers    
    \param PDESol (input) Solution vector referring to PDE node numbers
  */
  virtual void TransformElemSolution(Vector<Double> & MeshSol,  Vector<Double> & PDESol);


  //! assign  arrays for Mesh and PDE node numbers
  void AssignPDENodeNumbers();

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


  std::string pdename_; //!< type of PDE (set in the derived classes)
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

  // Assignment MeshNodeNumers <-> PDENodeNumbers
  std::vector<Integer> Mesh2PDENode_;  //!< array containing PDE (=local) node numbers
  std::vector<Integer> PDE2MeshNode_;  //!< array containing Mesh (=global) node numbers 
 
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

  Integer actlevel_; //! actual level
  TimeErrorEstimator * ptTimeError_; //!< pointer to extimator

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

  //Dummies, just for SUN compiler
  Vector<Double> DVec;
  Double Dval;
  Integer Ival;
  Boolean Dbool;
  TimeErrorEstimator *DtimeErrorEst;
};

} // end of namespace

#endif
