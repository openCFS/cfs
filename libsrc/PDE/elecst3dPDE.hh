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
Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse, Double &coarsealpha);

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

  //! test error of solution
  virtual Boolean TestError();

  //! refine mesh
  virtual void RefineMesh();

private:
  //! initialization of pointer to SpaceErrorEstimator
  void ConstructorError();

  //!  calculation of error for each cell of mesh.
  void CalcErrorMap();

   //! calculation of error for the element of mesh
  void CalcErrorForElem(const Elem* elem, const Vector<Double>* gradSPR, Double & error, Double & normGradSPR);

  //! calculation of electric field
  void CalcElectricField();

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

  //! indicator for WriteResults. TRUE, if we want to print error map.
  Boolean WriteErrorMap_;

  //! indicator that we have to calculate electric field
  Boolean calcElecField_;
  //! arrays we store values of electric field at the center of element
  Vector<Double> elecFieldAtCenterElem_[3];
 
  //! array, in which we store error map
  Vector<Double> errorMap_;
  //! array, in which we store l2 norm of gard SPR
  Vector<Double> gradSPRElemL2norm_;
  //! array, in which we store marked elements
  Vector<Double> markingElems_;
  //! array, in which we store value of relative error
  Vector<Double> relativeErrorMap_;
  //! error tolerance from conf-file
  Double errorTol_;
  //! l2 norm of SPR gradient
  Double normError_;
  
  //! pointer to class of error estimators
  SpaceErrorEstimator<3> * ptError_;
  
};

// class PutElemMatAlgSysElst3d
// {
// public:

//   PutElemMatAlgSysElst3d(AbstractAlgebraicSys * aptalgsys, Grid * aptgrid, const Double acoeffs, const Integer as_sysid, const Integer alevel)
//   { sysid_=as_sysid; ptalgsys_=aptalgsys; ptgrid_=aptgrid;
//     coeffst_=acoeffs; level_=alevel; matrix_stiff_=2;}

//   ~PutElemMatAlgSysElst3d() {}

//   // method
//   void operator() (Elem t);

// private:

//      //!
//   Grid * ptgrid_;

//      //!
//   AbstractAlgebraicSys * ptalgsys_;

//      //!
//   Integer sysid_, level_, matrix_stiff_;

//      //!
//   Double  coeffst_;

// };


} // end of namespace
#endif
