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

  //! Base class for partial differential equation
  /*! Class BasePDE is base class, from which different type of PDEs are derived. 
   */

class TimeErrorEstimator;
class SpaceErrorEstimator;

class BasePDE
{
public:

  //! Constructor( read integration parameters, define class Material)
  BasePDE(AbstractAlgebraicSys * aptalgsys, Material * aMatFile, FileType * aInFile, WriteResults * aOutFile, TimeFunc * aptTimeFunc); 

  //! Deconstructor
  virtual ~BasePDE();

  //!
  virtual void SetMatrixFactors()=0;

  //! define algebraic system identifictaion
  virtual void SetAlgSys_id(const Integer AS_sysid);

  //!
  virtual void SetupMatrices(const Integer type)=0;

  //!
  virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, 
                               Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)=0;

  //!
  virtual void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter, Integer &amaxnumit, Integer &numeqcoarse)=0;

  //!
  virtual void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atimestep)=0;

  //!
  virtual void SolveStepStatic(BCs * ptBCs, const Integer level)=0;

  //!
  virtual void SolveStepTrans(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean updatesysmat)
   { Error("Not implemented",__FILE__,__LINE__);}

  //!
  virtual void WriteResultsInFile()=0;  

  //! Calculation parameters in Newmark method
  virtual void CalcParameters(const Double adt)  
  { Error("Not implemented");}

  //! Create pointer to according class of time error estimation
  virtual TimeErrorEstimator * CreatePtTimeError(){ Error("Function CreatePtTimeError is not overloaded in this class");}  

  //! Create pointer to according class of space error estimation
  virtual SpaceErrorEstimator * CreatePtSpaceError(){ Error("Function CreatePtSpaceError is not overloaded in this class");}

  //!
  virtual const Vector<Double> & getS() const =0;

  //!
  virtual const Vector<Double> & getS1() const 
  { Error("Not implemented");}

    //!
  virtual const Vector<Double> & getS1old() const
  { Error("Not implemented");}

  //!
  virtual const Vector<Double> & getS2() const { Error("Function getS2 is not overloaded in this class");}

 //!
  virtual const Vector<Double> & getS2old() const { Error("Function getS2old is not overloaded in this class");}

  //!
  virtual Double getBeta() const { Error("Function getBeta is not overloaded in this class");}

  virtual Double getGamma() const { Error("Function getGamma is not overloaded in this class");}

  //!
  virtual Integer getSize() const { Error("Function getSize is not overloaded in this class");}  

  //!
  virtual AbstractAlgebraicSys * getAlgSys() const { return ptalgsys_;}

  virtual void CalcThirdDerivateFromEquation(Vector<Double> & result)
 { Error("Function getSize is not overloaded in this class");} 

  //!
  Integer GetSysId() const { return as_sysid_;} 

  //!
  Integer GetNumRestraints(BCs * ptBCs, const Integer level);

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
  SpaceErrorEstimator * ptSpaceError_;

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
  std::vector<Integer> val_id_; 

  //!
  std::vector<std::string> subdoms_;  

};

} // end of namespace

#endif
