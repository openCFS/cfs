#ifndef FILE_BASECOUPLEDPDE_2003
#define FILE_BASECOUPLEDPDE_2003


#include <list>
#include "Domain/bcs.hh"
#include "DataInOut/filetype.hh"
#include "DataInOut/writeresults.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

#include <PDE/basePDE.hh>

namespace CoupledField
{
 
class PDECoupling;

 //! Base class for coupling between different PDEs
  /*! Class BaseCoupledPDE is base class for coupled PDE problems
   */
class BaseCoupledPDE
{
public:

  //! Constructor
  /*!    
    \param PDEs vector with pointers to pdes
    \param Couplings pointer to vector of coupling objects
    \param aptrgrid pointer to grid object
    \param aptBCs pointer to boundary condition object
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param tag for current multisequence step
  */
  BaseCoupledPDE(StdVector<BasePDE*> & PDEs,
		 StdVector<PDECoupling*> & Couplings,
		 Grid *aptgrid, 
		 BCs *aptBCs, 
		 FileType *aInFile, 
		 WriteResults *aOutFile,
		 std::string sequenceTag); 

  //! Deconstructor
  virtual ~BaseCoupledPDE();


  virtual std::string GetName()
  {return coupledpdename_;};
  
  //! calculates coupling interfaces
  virtual void InitCoupling(Integer level)=0;
  
  //! Defines, which of the coupled PDEs are currently solved
  //! and which are neglected. This method is mainly needed
  //! for multiSequenceAnalysis, where in different steps
  //! only a subset of all coupled PDEs is solved.
  //! \param pdes (input) Unsorted list of PDEs which
  //! are currently solved 
  virtual void DefineSolvingPDEs(StdVector<BasePDE*> & pdes) = 0;
  
  //! Solve static step
  virtual void SolveStepStatic(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat)=0;
  
  //! solve transient step
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat)=0;
  
  //! write results in file
  virtual void WriteResultsInFile(Integer stepOffset = 0,
				  Double timeOffset = 0.0)=0;


  //! Init the time stepping
  /*!
    \param dt time step
  */
  virtual void InitTimeStepping(const Double dt);

  /// time stepping params, provided by "TransientDriver"
  void SetTimeSteppingParams(Integer numstep, Double  firstdt, Integer isavebegin, 
			     Integer isaveend, Integer isaveincr);
  
  /// write general definitions (loads, bcs, ...) to every pde
  virtual void WriteGeneralPDEdefines();


  /// perform on every pde a pre step (before solving transient step
  virtual void PreStepTrans(const Integer kstep, const Double asteptime, 
			    const Integer level, const Boolean reset);

  /// perform on every pde a post step (after solving transient step
  virtual void PostStepTrans(const Integer kstep, const Double asteptime, 
			     const Integer level);

  /// initialize PDEs before iteration (done for each time step)
  void InitStepTransCoupled(Double asteptime);
  


  // ======================================================
  // POSTPROC SECTION
  // ======================================================

  //! Do Postprocessing as descriped in conf file
  virtual void PostProcess(const Integer level);


protected:
  
  // general PDE parameters
  std::string coupledpdename_;        //!< name of the coupled pde 
  AnalysisType analysistype_;         //!< type of analysis
  StdVector<BasePDE *> PDEs_;         //!< list of belonging PDEs
  StdVector<PDECoupling*> Couplings_; //!< vector of coupling objects
  Integer NumPDEs_;                   //!< number of PDEs 
  Integer actlevel_;                  //!< current level (for multigrid)
  std::string sequenceTag_;           //!< tag for current multisequence step
  
  //! vector of flags indicating if specified
  //! PDE gets solved. The ordering corresponds
  //! to that of PDEs_ in the base class.
  StdVector<Boolean> solvePDE_;

  // pointers to objects
  Grid * ptgrid_;           //!< pointer to Grid
  BCs *ptBCs_;              //!< pointer to Boundary Condition  Object
  FileType * InFile_;       //!< pointer tio input file
  WriteResults * OutFile_;  //!< pointer to output file


  Double dt_;                //!< time step size  
  Integer numstep_;          //!< nr. of calculated timesteps
  Double  firstdt_;          //!< time step size
  Integer isavebegin_;       //!< starting index of saving
  Integer isaveend_;         //!< end index of saving
  Integer isaveincr_;        //!< incremental step of saving
};

} // end of namespace

#endif
