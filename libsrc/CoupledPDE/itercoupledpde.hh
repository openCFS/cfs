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
		 WriteResults *aOutFile); 

  //! Destructor
  ~IterCoupledPDE();

   //! calculates coupling interfaces
  virtual void InitCoupling(Integer level);
  
  //! Solve static step
  virtual void SolveStepStatic(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat);
  
  //! solve transient step
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat);
  
  //! write results in file
  virtual void WriteResultsInFile();


protected:
  //!
  void WriteCouplingInfo();

  //! calculates the norm of a vector
  Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval);

  Integer maxiter_;                        //!< maximum number of iterations per time step
  StdVector<Double> norms_;              //!< norm of coupling values
};

} // end of namespace

#endif
