#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include "basePDE.hh"

 
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


  /// return index of dof defined by keyword (e.g. 'ux')
  virtual Integer GetBCDof(const std::string keyword);
  

  /// calculates L2-norm of RHS regarding entries due to penalty formulation
  Double RhsL2Norm(Vector<Double>& stdVec);


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

/// do one transient step
  void StepTransNonLin(const Integer kstep, const Double asteptime,
		       const Integer level, const Boolean reset);
  

  //! prepare for correct time stepping
  /*! \param dt time step  */
  virtual void InitTimeStepping(const Double dt);

  //!
  virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			     const Integer level, const Boolean reset);

  //! solve one step for nonlinear static problem 
  /*!
    \param level level of grid
    \param aTime sequence of different levels for RHS
  */
  virtual void StepStaticNonLin(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean reset);

  //!
  virtual void PostStepStatic(const Integer kstep, const Double asteptime,
			      const Integer level);
  

  // ======================================================
  // POSTPROC SECTION
  // ======================================================

  //! write results in file
   virtual void WriteResultsInFile();

  //!  return pointer to vector with first derivative of solution
  //virtual const Array<Double>& getS1() const { return TS_alg_->GetDeriv1();}
  virtual const Vector<Double>& getS1() const { return TS_alg_->GetDeriv1();}

  //! return pointer to vector with second derivative of solution
  //virtual const Array<Double>& getS2() const { return TS_alg_->GetDeriv2();}
  virtual const Vector<Double> & getS2() const { return TS_alg_->GetDeriv2();}
protected:

  
  Integer size_;        //!< total number of unknowns (equations)


private:

  /// calc rhs coupling to acoustic pde
  //void CalcAcousticCouplingRHS(std::vector<Elem*> * couplingElems, Vector<Double>& forceOnElem);
  
  /// calc rhs coupling to acoustic pde
  void CalcAcousticCouplingRHS(std::vector<Elem*> * couplingElems, 
			       std::vector<Integer>& couplingNodes,
			       std::vector<MaterialData*>* materials,
			       //Array<Double>& forceOnElem,
			       StoreSol<Double> & forceOnElem,
			       Integer couplingdof,
			       std::vector<Elem*> * neighbours);
  

  /// does a line search and returns the optimal residual norm
  Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
		    Double& etaLineSearch, Integer level, Boolean trans=FALSE);


  /// Write nonlin iteration norms to the cla-file
  void WriteClaNlNorms(const Integer iterationCounter, const Double residualL2Norm,
		       const Double extForcesL2Norm, const Double residualErr, 
		       const Double solIncrL2Norm, const Double actSolL2Norm, 
		       const Double incrementalErr);
  

  /// calculates matrices D^_ and D^__ (see Hughes p. 217) for reduced integration
  void CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat, MaterialData& mat);


  // defines subtype of mechanic PDE: plainStrain, 3d, ...
  std::string subType_;

  Integer GetNrBCDof (const std::string & dofStartString);

  /// stores an algsys_ vector into a std::vector and returns that L2-norm
  void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);


  /// returns that L2-norm of an algsys vector
  Double AlgsysL2Norm(Double * pt);
  

  /// flag for reduced Integration
  Boolean reducedInt_;

  /// returns the solution matrix belonging to all nodes of the actual element
  void GetSolOfElement( Matrix<Double>& elDisp, Vector<Integer>& connect_PDE);  
  

  /// value of prestress
  Double preStressVal_;

  /// direction of prestress
  Directions preStressDir_;
  
  /// material data for reduced integration
  MaterialData * lambdaMat;

  /// material data for reduced integration
  MaterialData * mueMat;

  /// external forces (for nonlin simulations)
  Vector<Double> extForces_;

};

} // end of namespace
#endif

