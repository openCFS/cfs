#ifndef FILE_DIRECT_COUPLED_PDE
#define FILE_DIRECT_COUPLED_PDE

#include "PDE/StdPDE.hh"

namespace CoupledField {

  // forward class declaration
  class SinglePDE;
  
  //! This class implements the direct coupling of StdPDEs.

  //! This class handles the direct coupling of two or more
  //! SinglePDEs. Therefore it needs references to SinglePDEs
  //! and the according pair coupling objects.

class DirectCoupledPDE : public StdPDE
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
  DirectCoupledPDE(Grid *aptgrid, BCs *aptBCs, 
		   FileType *aInFile, 
		   WriteResults *aOutFile, 
		   TimeFunc *aTimeFunc);
    
  //! Destructor
  virtual ~DirectCoupledPDE();

  //! write general defines (BCs, loads, etc.) to info-file
  void WriteGeneralPDEdefines();

  
  //! set boundary condition
  //! \param level             level of grid
  //! \param atimestep         time step of claculation
  void SetBCs(const Integer level, 
	      const Double atimestep);

  // ======================================================
  // POSTPROC SECTION
  // ======================================================
  
  //@{
  //! \name Methods performing post-processing
  
  //! Do Postprocessing as descriped in conf file
  void PostProcess(const Integer level);
  
  //! write results in file
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time  
  void WriteResultsInFile(const Integer kstep = 0,
			  const Double asteptime = 0.0,
			  Integer stepOffset = 0,
			  Double timeOffset = 0.0);
  //@}


  // ======================================================
  // COUPLING SECTION
  // ======================================================
  
  //! initalize PDE coupling (only done once)
  void InitCoupling(PDECoupling * Coupling);

  //! reset coupling counters and data (done after each timestep)
  void ResetCoupling();
  
  //! Fill in input coupling terms
  void CalcInputCoupling();
  
  
  //! calculate coupling terms
  void CalcOutputCoupling();

  
private:

  //! References to SinglePDEs
  StdVector<SinglePDE*> singlePDEs_;

};

} // end of namespace

#endif
