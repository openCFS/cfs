#ifndef FILE_ELECTST3DPDE_2001
#define FILE_ELECTST3DPDE_2001

#include "basepde.hh"
 
namespace CoupledField
{

  //! Class for electrostatic equation in 3D
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 3D. 
  */

class Elecst3dPDE: virtual public BasePDE
{
public:

  //!
  Elecst3dPDE(AbstractAlgebraicSys * aptalgsys, Grid * , Material * , TimeFunc * ,FileType * , WriteResults * );

  //!
  virtual ~Elecst3dPDE();

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps,
Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse);

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(const Integer level=0);

    //!
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS();

  //! calculation derivates of solution 
  void CalculationDerivativesSol();

  //! Calculation of energy norm
  Double CalcEnergyNorm();

  //!
  void SolveStepStatic(BCs * ptBCs ,const Integer level);

  //!
   void WriteResultsInFile();

  //!
  virtual const Vector<Double> & getS() const { return sol_;}

  //!
  virtual Integer getSize() const { return size_;}

private:
  //!
  Integer dofspernode_;

  //!
  void CalcCoeff(Vector<Double> & coeff);

  //!
  Grid * ptgrid_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_;  

  //! size of solution and etc.
  Integer size_;

};

class PutElemMatAlgSysElst3d
{
public:

  PutElemMatAlgSysElst3d(AbstractAlgebraicSys * aptalgsys, Grid * aptgrid, const Double acoeffs, const Integer as_sysid, const Integer alevel)
  { sysid_=as_sysid; ptalgsys_=aptalgsys; ptgrid_=aptgrid;
    coeffst_=acoeffs; level_=alevel; matrix_stiff_=2;}

  ~PutElemMatAlgSysElst3d() {}

  // method
  void operator() (Elem t);

private:

     //!
  Grid * ptgrid_;

     //!
  AbstractAlgebraicSys * ptalgsys_;

     //!
  Integer sysid_, level_, matrix_stiff_;

     //!
  Double  coeffst_;

};


} // end of namespace
#endif
