#ifdef NEWBASEPDE

#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include "newBasePDE.hh"

 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
    It is used for solving mechanic equation on one time step.  
  */

class MechPDE: public BasePDE
{

public:

  //!  Constructor. here we read integration parameters
  /*!
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  MechPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~MechPDE();



  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);


    //! set boundary condition
  /*!
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system or set new
    \param atimestep time step of calculation
  */
  virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);


  /// return index of dof defined by keyword (e.g. 'ux')
  Integer GetBCDof(const std::string keyword);
  

  /// calculates L2-norm of RHS regarding entries due to penalty formulation
  Double RhsL2Norm(std::vector<Double>& stdVec);


  /// sets external forces and returns L2Norm of them
  Double SetExternalForces(const Integer level);


  /// reads the directions (e.g. for prestress) from the config-file
  void GetDirection(Directions& dir, const std::string keyword);


  /// returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
  BaseForm * GetStiffIntegrator(MaterialData& actSDMat, Boolean reducedInt=FALSE);
  

  // ======================================================
  // COUPLING SECTION
  // ======================================================
  
 //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);

  //! calculate coupling terms
  virtual void CalcOutputCoupling();
  
  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output);


  /// setup source term
  void SetupRHS(const Integer level);
  

// ======================================================
// SOLVING SECTION
// ======================================================

  
  //! solve one step for static problems
  /*! \param level level of grid  */
  virtual void SolveStepStatic(const Integer level);
  

  //! solve one step for transient problem 
  /*!
    \param kstep number of calculating step
    \param steptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */
  virtual void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
			      const Boolean updatesysmat);
  

  //! prepare for correct time stepping
  /*! \param dt time step  */
  virtual void InitTimeStepping(const Double dt);


  // ======================================================
  // POSTPROC SECTION
  // ======================================================

  //! write results in file
   virtual void WriteResultsInFile();



protected:

  
  Integer size_;        //!< total number of unknowns (equations)


private:
  /// calculates matrices D^_ and D^__ (see Hughes p. 217) for reduced integration
  void CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat, MaterialData& mat);


  // defines subtype of mechanic PDE: plainStrain, 3d, ...
  std::string subType_;

  Integer GetNrBCDof (const std::string & dofStartString);

  /// stores an algsys_ vector into a std::vector and returns that L2-norm
  void StoreAlgsysToVec(std::vector<Double>& stdVec, Double * pt);


  /// returns that L2-norm of an algsys vector
  Double AlgsysL2Norm(Double * pt);
  

  /// flag for nonlinear calculations
  Boolean nonLin_;

  /// flag for reduced Integration
  Boolean reducedInt_;
  


  //! solve one step for linear static problem 
  /*!
    \param level level of grid
  */
  virtual void SolveStepStaticLin(const Integer level);
  


  //! solve one step for nonlinear static problem 
  /*!
    \param level level of grid
  */
  virtual void SolveStepStaticNonLin(const Integer level);

  /// returns the solution vector belonging to all nodes of the actual element
  void GetSolOfElement( Matrix<Double>& elDisp, Vector<Integer>& connect_PDE);

  /// stopping criterion for incremental error
  Double incStopCrit_;

  /// stopping criterion for residual error
  Double residualStopCrit_;  

  /// value of prestress
  Double preStressVal_;

  /// direction of prestress
  Directions preStressDir_;
  
  /// dof (e.g. ux) of homogenous Dirichlet BC
  std::vector<std::string> homDirichDof_; 

  /// dof (e.g. ux) of homogenous Dirichlet BC
  std::vector<std::string> inhomDirichDof_; 

  /// dof (e.g. ux) of load condition
  std::vector<std::string> loadDof_; 


  /// material data for reduced integration
  MaterialData * lambdaMat;

  /// material data for reduced integration
  MaterialData * mueMat;
};

} // end of namespace
#endif

#endif //#ifdef NEWBASEPDE
