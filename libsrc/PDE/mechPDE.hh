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
  MechPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile,
	  WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~MechPDE();

  //! Initialize NonLinearities
  virtual void InitNonLin();


  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);


  /// calculates L2-norm of RHS regarding entries due to penalty formulation
  Double RhsL2Norm(Vector<Double>& stdVec);


  /// sets external forces and returns L2Norm of them
  Double SetExternalForces(const Integer level);


  /// reads the directions (e.g. for prestress) from the config-file
  void GetDirection(Directions& dir, const std::string keyword);


  /// returns a stiffness integrator appropriate to the actual problem (e.g.3D)
  BaseForm * GetStiffIntegrator(MaterialData& actSDMat,
				Boolean reducedInt=FALSE);
  

  // ======================================================
  // COUPLING SECTION
  // ======================================================
  
  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);

  //! calculate coupling terms
  virtual void CalcOutputCoupling();
  
  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(SolutionType output);


  /// setup source term
  void SetupRHS(const Integer level);
  

  // ======================================================
  // SOLVING SECTION
  // ======================================================

  /// do one transient step
  void StepTransNonLin(const Integer kstep, const Double asteptime,
		       const Integer level, const Boolean reset);
  
  //! Init the time stepping
  //! \param dt time step
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
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time  
  virtual void WriteResultsInFile(Integer stepOffset = 0,
				  Double timeOffset = 0.0);

  //! do PostProcessing step
  virtual void PostProcess(const Integer level);

  //!  return pointer to vector with first derivative of solution
  //virtual const Array<Double>& getS1() const { return TS_alg_->GetDeriv1();}
  virtual const Vector<Double>& getS1() const { return TS_alg_->GetDeriv1();}

  //! return pointer to vector with second derivative of solution
  //virtual const Array<Double>& getS2() const { return TS_alg_->GetDeriv2();}
  virtual const Vector<Double> & getS2() const { return TS_alg_->GetDeriv2();}

protected:

  
  Integer size_;        //!< total number of unknowns (equations)

#ifdef XMLPARAMS
    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the mechanics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the MechPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Keyword</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>displacement</td>
    //!     <td>nodeResults</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>velocity</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>acceleration</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();
#endif


private:

  // calc rhs coupling to acoustic pde
  // void CalcAcousticCouplingRHS(StdVector<Elem*> * couplingElems,
  // Vector<Double>& forceOnElem);
  
  /// calc rhs coupling to acoustic pde
  void CalcAcousticCouplingRHS(StdVector<Elem*> * couplingElems, 
			       StdVector<Integer>& couplingNodes,
			       StdVector<MaterialData*>* materials,
			       Vector<Double> & forceOnElem,
			       Integer couplingdof,
			       StdVector<Elem*> * neighbours);
  

  /// does a line search and returns the optimal residual norm
  Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
		    Double& etaLineSearch, Integer level, Boolean trans=FALSE);


  /// Write nonlin iteration norms to the cla-file
  void WriteClaNlNorms(const Integer iterationCounter,
		       const Double residualL2Norm,
		       const Double extForcesL2Norm, const Double residualErr, 
		       const Double solIncrL2Norm, const Double actSolL2Norm, 
		       const Double incrementalErr);
  

  /// calculates matrices D^_ and D^__ (see Hughes p. 217) for reduced integration
  void CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat,
		      MaterialData& mat);

  Integer GetNrBCDof (const std::string & dofStartString);

  /// stores an algsys_ vector into a StdVector and returns that L2-norm
  void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

  /// returns that L2-norm of an algsys vector
  Double AlgsysL2Norm(Double * pt);
  
  /// flag for reduced Integration
  Boolean reducedInt_;

  /// returns the solution matrix belonging to all nodes of the actual element
  void GetSolOfElement( Matrix<Double>& elDisp, StdVector<Integer>& connect_PDE);

  StdVector<std::string> pressSurf_;  //!< surface of pressure loads
  StdVector<Double>      pressVals_;  //!< values of the pressure loads
  StdVector<std::string> pressFnc_;   //!< function names of pressure loads

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

  //postprocessing
  ElemStoreSol<Double> Stress_;  //!< conatins magnetic field
  StdVector<std::string> calcStress_;  //!< contains the subdomains, on which the stress is computed

  //! contains mechanic velocity
  NodeStoreSol<Double> solDeriv1_;
  
  //! contains mechanic acceleration
  NodeStoreSol<Double> solDeriv2_;

};

} // end of namespace
#endif

