#ifndef FILE_BASEPDE_2001
#define FILE_BASEPDE_2001

#include "abstractAlgSys.hh"
#include "bcs.hh"
#include "timefunc.hh"
#include "filetype.hh"
#include "writeresults.hh"
#include "material.hh"
#include "conffile.hh"


namespace CoupledField
{

 
class TimeErrorEstimator;

template<Integer dim>
class SpaceErrorEstimator;

 //! Base class for partial differential equation
  /*! Class BasePDE is base class, from which different type of PDEs are derived. 
   */
class BasePDE
{
public:

  //! Constructor. here we read integration parameters
  /*!
    \param  aptalgsys pointer to class Algebraic system
    \param aMatFile pointer to class Material. material data.
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  BasePDE(AbstractAlgebraicSys * aptalgsys, Material * aMatFile, FileType * aInFile,  WriteResults * aOutFile, TimeFunc * aTimeFunc); 

  //! Deconstructor
  virtual ~BasePDE();

  //! set information for algebraic system about PDE. set matrix factors.
  virtual void SetMatrixFactors()=0;

  //! define algebraic system identifictaion
  /*!
    \param AS_sysid id of PDE in algebraic system
   */
  virtual void SetAlgSys_id(const Integer AS_sysid);

  //! setup element matrices for AlgebraicSystem for assembling procedure
  /*!
    \param level level of grid
   */
  virtual void SetupMatrices(const Integer level, BCs * ptBCs=NULL)=0;

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param matrixtype out: 0..NOCLASS, 1..RSPARSE, 2..CSPARSE, 3..RBLOCK, 4.. CBLOCK = 0,
 5..RFULL, 6..CFULL, 7..MIXED
    \param matrixsystype out:define need we memory for different types of element-matrix or not                     
    \param graphtype out: type of graph
    \param numdofpernode out: number of dof per node
    \param numdirichlets out:number of nodes for dirichlets conditions
    \param numconstraints out:number of nodes for constraints conditions
  */
  virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, 
                               Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)=0;

  //! specify type of solver for algebraic system. it is read from config-file
  /*!
    \param asolvertype  Richardson or CG
    \param aprecondtype ID or MG
    \param aeps relative accuracy in the precond. energy
    \param adampiter damping parameter for Jacobi, SSOR
    \param amaxnumit max number of iterations
    \param numeqcoarse number of equation for coarsing
    \param coarsealpha coarsing parameter for AMG
  */
  virtual void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter, Integer &amaxnumit, Integer &numeqcoarse, Double &coarsealpha)=0;

  //! set boundary condition
  /*!
    \param ptBCs pointer to boundary condition
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atimestep time step of claculation
  */
  virtual void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atimestep)=0;

  //! solve one step for static problem 
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param level level of grid
  */
  virtual void SolveStepStatic(BCs * ptBCs, const Integer level)=0;

  //! solve one step for transient problem 
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param kstep number of calculating step
    \param asteptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */
  virtual void SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean updatesysmat)
   { Error("Not implemented",__FILE__,__LINE__);}

  //! solve one step for transient problem on new mesh. it is used in adaptive procedures for space
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param kstep number of calculating step
    \param asteptime time of calculation
    \param level level of grid
  */
   virtual void SolveStepTransNewMesh(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level)
   { Error("Not implemented",__FILE__,__LINE__);}

  //! write results in file
  virtual void WriteResultsInFile()=0;  

  //! Calculation parameters in Newmark method
  /*!
    \param adt size of timestep
  */
  virtual void CalcParameters(const Double adt)  
  { Error("Not implemented");}

  //! Create pointer to according class of time error estimation
  virtual TimeErrorEstimator * CreatePtTimeError(){ Error("Function CreatePtTimeError is not overloaded in this class");}  

  //! return pointer to vector with subdomains, on which we calculate the PDE
  virtual std::vector<std::string> * getSDsPDE(){ return &subdoms_;}

  //! return pointer to vector with solution
  virtual const Vector<Double> & getS() const =0;

  //! return pointer to vector with first derivative of solution
  virtual const Vector<Double> & getS1() const 
  { Error("Not implemented");}

    //! return pointer to vector with first derivative of solution, calculated on previous step
  virtual const Vector<Double> & getS1old() const
  { Error("Not implemented");}

  //! return pointer to vector with second derivative of solution
  virtual const Vector<Double> & getS2() const { Error("Function getS2 is not overloaded in this class");}

 //! return pointer to vector with second derivative of solution, calculated on previous step
  virtual const Vector<Double> & getS2old() const { Error("Function getS2old is not overloaded in this class");}

  //! return parameter beta from Newmark method
  virtual Double getBeta() const { Error("Function getBeta is not overloaded in this class");}

  //! return parameter gamma from Newmark method
  virtual Double getGamma() const { Error("Function getGamma is not overloaded in this class");}

  //! return size of solution
  virtual Integer getSize() const { Error("Function getSize is not overloaded in this class");}  

  //! return pointer to algebraic system
  virtual AbstractAlgebraicSys * getAlgSys() const { return ptalgsys_;}

  //! return id of PDE in algebraic system
  Integer GetSysId() const { return as_sysid_;} 

  //! initialize pointer to algebraic system
  void InitPtAlgSys(AbstractAlgebraicSys * aptalgsys) { ptalgsys_=aptalgsys;}
  
  //! return number of restraints
  Integer GetNumRestraints(BCs * ptBCs, const Integer level=-1);

  //! restore solution from previous step. it is used in time adaptive procedure
  virtual void RestoreSol()
  { Error("Function RestoreSol is not implemented in this class");}  

  //! save received solution as solution on the previous step 
  virtual void SaveSolAsPrevStep()
  { Error("Function RestoreSol is not implemented in this class");}  

  //! calculate energy norm
  virtual  Double CalcEnergyNorm()
   { Error("Function PDE::CalcEnergyNorm is not implemented in this class");}

  //! refine mesh
  virtual void RefineMesh()
  { Error("Function BasePDE::RefineMesh is not implemented in this class");} 

    //! test space error
  virtual Boolean TestError()
    { Error("Function BasePDE::TestError is not implemented in this class");}  

 //! write additional info (marked elements, relative error) to files with mesh
  virtual void PrintMeshesInfo(WriteResults * ptMehes)
  { Error("Function BasePDE::PrintMeshesInfo is not implemented in this class");}

protected:
   //! read from .config-file info about BCs
   void ReadBCs(const std::string eq);

  //! 
  Double StepTime_;

  //! pointer to class Material
  Material * MatFile_;

  //!
  FileType * InFile_;

  //!
  WriteResults * OutFile_;

  //!
  TimeFunc * ptTimeFunc_;

 //!
  Double matrix_factor_[4];

 //!
  TimeErrorEstimator * ptTimeError_;

 //!
  AbstractAlgebraicSys * ptalgsys_;

   //!
  Integer as_sysid_;

  //! boundary conditions
  std::vector<std::string> bcs_hd_;
  std::vector<std::string> bcs_id_;
  std::vector<std::string> bcs_nh_;
  std::vector<std::string> bcs_ni_;
  std::vector<std::string> bcs_rh_;
  std::vector<std::string> bcs_ri_;

  //!
  std::vector<Double> val_id_; 

  //!
  std::vector<std::string> subdoms_;  

};

} // end of namespace

#endif
