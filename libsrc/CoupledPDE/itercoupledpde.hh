#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include "CoupledPDE/basecoupledpde.hh"
//#include "Utils/storesol.hh"
//#include "Utils/elemstoresol.hh"

namespace CoupledField
{

  //! Derived Class from BaseCoupledPDE. It solves iteratively
  //! a list of PDES and calculates coupling terms in between

class IterCoupledPDE : public BaseCoupledPDE
{

 public:

  //! Constructor
  IterCoupledPDE(StdVector<BasePDE*> & PDEs,
		 StdVector<PDECoupling*> & Couplings,
		 Grid *aptgrid, 
		 BCs *aptBCs, 
		 FileType *aInFile, 
		 WriteResults *aOutFile,
		 std::string sequenceTag); 

  //! Destructor
  ~IterCoupledPDE();

  //! calculates coupling interfaces
  virtual void InitCoupling(Integer level);
  
  //! Defines, which of the coupled PDEs are currently solved
  //! and which are neglected. This method is mainly needed
  //! for multiSequenceAnalysis, where in different steps
  //! only a subset of all coupled PDEs is solved.
  //! \param pdes (input) Unsorted list of PDEs which
  //! are currently solved 
  virtual void DefineSolvingPDEs(StdVector<BasePDE*> & pdes);

  //! Solve static step
  virtual void SolveStepStatic(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat);
  
  //! solve transient step
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat);
  
  //! write results in file
  virtual void WriteResultsInFile(Integer stepOffset = 0,
				  Double timeOffset = 0.0);


protected:

  //! Write information about coupling interfaces
  //! and coupling setup into ostream
  void WriteCouplingInfo(std::ostream &out);

  //! calculates the norm of a vector
  Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval);

  Integer maxiter_;                        //!< maximum number of iterations per time step
  StdVector<Double> norms_;              //!< norm of coupling values

  //! Flag for nonlinear logging
  Boolean nonLinLogging_;
  
};

} // end of namespace

#endif
