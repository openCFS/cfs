#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include "PDE/basePDE.hh"

namespace CoupledField
{

  // forward class declaration
  class IterSolveStep;


  //! Derived Class from BasePDE. It solves iteratively
  //! a list of PDES and calculates coupling terms in between

  class IterCoupledPDE : public BasePDE
{

  // friend declaration
  friend class IterSolveStep;

 public:

  //! Constructor
  IterCoupledPDE(StdVector<StdPDE*> & PDEs,
		 StdVector<PDECoupling*> & Couplings,
		 std::string sequenceTag); 

  //! Destructor
  ~IterCoupledPDE();

  //! calculates coupling interfaces
  void InitCoupling(Integer level);
  
  //! write general defines (BCs, loads, etc.) to info-file
  void WriteGeneralPDEdefines();
  
  //! Defines, which of the coupled PDEs are currently solved
  //! and which are neglected. This method is mainly needed
  //! for multiSequenceAnalysis, where in different steps
  //! only a subset of all coupled PDEs is solved.
  //! \param pdes (input) Unsorted list of PDEs which
  //! are currently solved 
  void DefineSolvingPDEs(StdVector<StdPDE*> & pdes);

  //! Solve static step
  //void SolveStepStatic(const Integer kstep, const Double asteptime, const Integer level, 
  //	       const Boolean updatesysmat);
  
  //! solve transient step
  //void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
  //		      const Boolean updatesysmat);
  
  //! write results in file
  void WriteResultsInFile(Integer stepOffset = 0,
			  Double timeOffset = 0.0);
  
  //! set time step
  //! \params dt Current time step
  void SetTimeStep(const Double dt);

  //! Do Postprocessing as descriped in conf file
  void PostProcess(const Integer level);
  
protected:

  //! Write information about coupling interfaces
  //! and coupling setup into ostream
  void WriteCouplingInfo(std::ostream &out);

  //! calculates the norm of a vector
  //Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval);

  Integer maxiter_;                        //!< maximum number of iterations per time step
  StdVector<Double> norms_;              //!< norm of coupling values

  //! Flag for nonlinear logging
  Boolean nonLinLogging_;
  
  // general PDE parameters
  AnalysisType analysistype_;         //!< type of analysis
  StdVector<StdPDE *> PDEs_;         //!< list of belonging PDEs
  StdVector<PDECoupling*> Couplings_; //!< vector of coupling objects
  Integer NumPDEs_;                   //!< number of PDEs 
  Integer actlevel_;                  //!< current level (for multigrid)
  std::string sequenceTag_;           //!< tag for current multisequence step
  
  //! vector of flags indicating if specified
  //! PDE gets solved. The ordering corresponds
  //! to that of PDEs_ in the base class.
  StdVector<Boolean> solvePDE_;

};

} // end of namespace

#endif
