#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include <CoupledPDE/basecoupledpde.hh>

namespace CoupledField
{

  //! Derived Class from BaseCoupledPDE. It solves iteratively
  //! a list of PDES and calculates coupling terms in between

class IterCoupledPDE : public BaseCoupledPDE
{
public:

  //! Constructor
  IterCoupledPDE(std::vector<BasePDE*> & PDEs,
		 std::vector<PDECoupling*> & Couplings,
		 Grid *aptgrid, 
		 BCs *aptBCs, 
		 FileType *aInFile, 
		 WriteResults *aOutFile); 

  //! Destructor
  ~IterCoupledPDE();
  
  //! calculates coupling interfaces
  virtual void InitCoupling(Integer level);
  
  //! Solve static step
  virtual void SolveStepStatic(const Integer level);
  
  //! solve transient step
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat);
  
  //! write results in file
  virtual void WriteResultsInFile();

};

} // end of namespace

#endif
