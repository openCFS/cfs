#ifndef FILE_BASECOUPLEDPDE_2003
#define FILE_BASECOUPLEDPDE_2003


#include <list>
#include <Domain/bcs.hh>
#include <DataInOut/filetype.hh>
#include <DataInOut/writeresults.hh>

#ifndef NEWBASEPDE
#include <PDE/basepde.hh>
#else
#include <PDE/newBasePDE.hh>
#endif //#ifndef NEWBASEPDE

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
  */
  BaseCoupledPDE(std::vector<BasePDE*> & PDEs,
		 std::vector<PDECoupling*> & Couplings,
		 Grid *aptgrid, 
		 BCs *aptBCs, 
		 FileType *aInFile, 
		 WriteResults *aOutFile); 

  //! Deconstructor
  virtual ~BaseCoupledPDE();


  //! name
  virtual std::string GetName()
  {return coupledpdename_;};
  
  //! calculates coupling interfaces
  virtual void InitCoupling(Integer level)=0;
  
  //! Solve static step
  virtual void SolveStepStatic(const Integer level)=0;
  
  //! solve transient step
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat)=0;
  
  //! write results in file
  virtual void WriteResultsInFile()=0;


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
  virtual void PreStepTrans(const Integer level, const Boolean reset);

  /// perform on every pde a post step (after solving transient step
  virtual void PostStepTrans(const Integer level);

  /// initialize PDEs before iteration (done for each time step)
  void InitStepTransCoupled(Double asteptime);
  


  // ======================================================
  // POSTPROC SECTION
  // ======================================================

  //! Do Postprocessing as descriped in conf file
  virtual void PostProcess(const Integer level);


protected:
  
  // general PDE parameters
  std::string coupledpdename_;          //!< name of coupled PDE
  AnalysisType analysistype_;           //!< type of analysis
  std::vector<BasePDE *> PDEs_;         //!< list of belonging PDEs
  std::vector<PDECoupling*> Couplings_; //!< vector of coupling objects
  Integer NumPDEs_;                     //!< number of PDEs 
  Integer actlevel_;                    //!< current level (for multigrid)

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
