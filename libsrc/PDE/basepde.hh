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

class BasePDE
{
public:

  //! Constructor( read integration parameters, define class Material)
  BasePDE(Material * aMatFile, FileType * aInFile, WriteResults<Point2D> * aOutFile, TimeFunc * aptTimeFunc); 

  //! Deconstructor
  virtual ~BasePDE(){;}

  //!
  virtual void SetMatrixFactors()=0;

  //! define algebraic system identifictaion
  virtual void SetAlgSys_id(Integer AS_sysid)=0;

  //!
  virtual void SetupMatrices(AbstractAlgebraicSys * algsys, Integer type)=0;

  //!
  virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, 
                               Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints)=0;

  //!
  virtual void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter, Integer &amaxnumit)=0;

  //!
  virtual void SetBCs(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer level, const Integer update, const Double atimestep)=0;

  //!
  virtual void SolveStepStatic(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, Integer level)=0;

  //!
  virtual void SolveStepTrans(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level, const Boolean updatesysmat)=0;

  //!
  virtual void WriteResultsInFile()=0;  

  //! Calculation of integration parameteres 
//  void CalcIntegrationParam(const Double dt);

  //! Calculation parameters in Newmark method
  virtual void CalcParameters(const Double adt)=0;  

  //!
  virtual Vector<Double> & getS()=0;

  //!
  virtual Vector<Double> & getS1()=0;

  //!
  virtual Vector<Double> & getS2()=0;
 
protected:

  //!
//  Integer numnode;

  //!
  Double StepTime_;

  //! pointer to class Material
  Material * MatFile_;

  //!
  FileType * InFile_;

  //!
  WriteResults<Point2D> * OutFile_;

  //!
  TimeFunc * ptTimeFunc_;

 //!
  Double matrix_factor_[4];

};

} // end of namespace

#endif
