#ifndef FILE_SINGLEPDE
#define FILE_SINGLEPDE


#include "PDE/StdPDE.hh"


#include <list>

#include "Utils/StdVector.hh"
#include "Domain/elem.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "General/environment.hh"
#include "Domain/bcs.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/filetype.hh"
#include "DataInOut/writeresults.hh"
#include "Matrix/matrix.hh"

#include "olas.hh"


#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/LoadMaterialDataFile.hh"
#include "DataInOut/MaterialData.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Driver/assemble.hh"
#include "timestepping.hh"
#include "nodeEQN.hh"
#include "pdememento.hh"
#include "Driver/baseSolveStep.hh"

#ifdef USE_DATABASE
#include "DataInOut/LoadMaterialDataDatabase.hh"
#endif


namespace CoupledField
{
  // forward class declaration
  class SpaceErrorEstimator;
  
  //! Base class for all kinds of single field problems.

  //! This class serves as base class for all single field problems, 
  //! like electrostatic,  acoustic, mechanic and others.

class SinglePDE : public StdPDE
{
  
public:

  Boolean BooleanComplexMaterialData_;

  Boolean converged_; //!< needed for coupling with MpCCI
 

 
  void Init(Integer sequenceStep = 0,
	    std::string  bcSequenceTag = "anyTag");
  
  // ---------------------- ***** --------------------------------
  

  //! returns the load names
  StdVector<std::string>& GetLoadDom()
  {return assemble_->loadDom_;};
  
  //! returns the load dofs
  StdVector<std::string>& GetLoadDof()
  {return assemble_->loadDof_;};
  
  //! returns the load values
  StdVector<Double>& GetLoadVals()
  {return assemble_->loadVals_;};
  
  //!returns the load functions
  StdVector<std::string>& GetLoadFncs()
  {return assemble_->fncname_loads_;};
  
  //! MpCCI gets the geometry
  virtual void PreparePDE4Computation() {;};
  
  // ======================================================
  // ADATPTIVITY SECTION
  // ======================================================
#ifdef ADAPTGRID
  //@{
  //! \name Methods used for performing adaptivity
  
  //! test error of calculation. return TRUE, if it is more then tolerance
  virtual Boolean TestError(const Integer level);
  
  //! refine mesh
  virtual void RefineMesh(const Integer level=0);

  // Does this method belong to postproc section?
  //! write information about relative error of calculation
  void WriteErrorInfo(WriteResults * ptmeshes);
  
  //@}
#endif

  // ======================================================
  // COUPLING SECTION
  // ======================================================
  
  //! initalize PDE coupling (only done once)
  virtual void InitCoupling(PDECoupling * Coupling) = 0;

  //! reset coupling counters and data (done after each timestep)
  virtual void ResetCoupling();

  //! Fill in input coupling terms
  virtual void CalcInputCoupling();
  
  
  //! calculate coupling terms
  virtual void CalcOutputCoupling() = 0;

  
  // ======================================================
  // GETTER METHODS
  // ======================================================

  //! return subtype
  virtual std::string GetSubType() {return subType_;}

  //! return number of restraints
  Integer GetNumRestraints(const Integer level=-1);
 
  //! set boundary condition
  //! \param level             level of grid
  //! \param atimestep         time step of claculation
  void SetBCs(const Integer level, 
	      const Double atimestep);
  
  //! write general defines (BCs, loads, etc.) to info-file
  void WriteGeneralPDEdefines();
  
protected:

  
  //! Constructor
  /*!
    \param aptgrid pointer to grid
    \param aptBCs pointer to boundary condition object
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  SinglePDE(Grid *aptgrid, 
	    BCs *aptBCs, 
	    FileType *aInFile,
	    WriteResults *aOutFile, 
	    TimeFunc *aTimeFunc);


  //! destructor
  virtual ~SinglePDE();

  //! private copy constructor
  SinglePDE & operator= (const StdPDE & myPDE) {
    ;}

  // ======================================================
  // INITIALIZATION METHODS
  // ======================================================
  
  //! Obtain information on desired output quantities from parameter file
  //! This method is used to query the parameter handling object for the
  //! desired output quantities and translate their literal description into
  //! the internal format by setting the corresponding class attributes.
  virtual void ReadStoreResults() = 0;


  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level)=0;

  //! read material data
  virtual void ReadMaterialData();

  //! read from config-file info about BCs
  void ReadBCs();

  // overloaded version of ReadBCs for special
  // boundary conditions in derived classes
  virtual void ReadSpecialBCs(){}

  //! Initialize NonLinearities
  virtual void InitNonLin(){};

  //! Init the time stepping
  virtual void InitTimeStepping()
  {Error("InitTimeStepping not implemented",__FILE__,__LINE__);};

#ifdef ADAPTGRID  
  //! ----------------- functions for adaptivity
  //! in this fnc we delete old pointer to Error-object & create new one
  void ConstructorError();
#endif
    

  // ======================================================
  // DATA SECTION
  // ======================================================

  
  // -----------------------------------------------------------------------
  // Storing information
  // -----------------------------------------------------------------------
  
  //@{
  //! \name Attributes connected to storing information
    
  //! vector containing solutiontypes of PDE
  StdVector<SolutionType> solTypes_;

  //! vector containgin dofs of solutiontypes
  StdVector<Integer> solDofs_;

  //! TRUE, if at least one Result is calculated and written
  Boolean hasOutput_;

  //! TRUE, if solution should be written to result file
  Boolean saveSol_;

  //! TRUE, if first derivative of solution should be written to result file
  Boolean saveDeriv1_;

  //! TRUE, if second derivative of solution should be written to result file
  Boolean saveDeriv2_;

  //! TRUE, if solution should be written to history file
  Boolean saveSolHist_;

  //! TRUE, if first derivative of solution should be written to history file
  Boolean saveDeriv1Hist_;

  //! TRUE, if second derivative of solution should be written to history file
  Boolean saveDeriv2Hist_;

  //@}

  // -----------------------------------------------------------------------
  // Adaptivity
  // -----------------------------------------------------------------------

  //@{
  //! \name Attributes connected to adaptivity

  //! object with methods for error estimation
  SpaceErrorEstimator * ptError_;

  //! array where  we store the number of refinement for the element
  StdVector<Double> markingElems_;

  Vector<Double> errorMap_;  //!< array with error map
  Double tolSpaceErr_;       //!< tolerance
  //@}

  // -----------------------------------------------------------------------
  // Solver parameters
  // -----------------------------------------------------------------------

  //@{
  //! \name Attributes connected to parameters for solver
  Integer maxnumiter_;       //!< maximum of iterations (for iterative solver)
  Integer numeqcoarse_;      //!< number of unknowns on coarse level(just for AMG)
  Double  eps_;              //!< accuracy
  Double dampiter_;          //!< damping parameter within iterative solution
  Double coarsealpha_;       //!< coarsening factor (just for AMG)
   
  ComplexFormat complexFormat_;  //!< outputFormat for complex numbers

  //@}

};

} // end of namespace
#endif
