#ifndef FILE_ACOU2DFLOWNOISE_2003
#define FILE_ACOU2DFLOWNOISE_2003

#include "acoustic2dPDE.hh"

namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class Acoustic2dPDE. It is used for solving the acoustic analogy inhomogeneous acoustic equation for flow induced noise on one time step.  We define the right hand side by setting the sound sources using the results from the external fluid computation, set rules for assembling global system matrix according to weak form of this PDE, define right hand side and set first order absorbing boundary conditions. Then we call one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

class Acou2dFlowNoise: public Acoustic2dPDE
{
public:

  //!
  Acou2dFlowNoise(AbstractAlgebraicSys * ptalgsys, Grid * aptgrid, Material * ptMaterial, TimeFunc * aptTimeFunc, FileType * aptFileType, WriteResults * aptOut);

  //!
  virtual ~Acou2dFlowNoise();

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(const Integer level, BCs * ptBCs=NULL);

    //!
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS(const Double atime, BCs * ptBCs=NULL);

  //!
  void preComputeRHS();

  //! Reads at every time the flowdatafile from the Fluids Department LSTM.
  void ReadFlowData(const char * aname, const Integer timestep, Matrix<Double> &nodedata );

  //!
  void SolveStepStatic(BCs * ptBCs ,const Integer level);

  //!
  void SolveStepTrans(BCs * ptBCs ,const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

  //!
   void WriteResultsInFile();

private:

  //!
  void CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff, Vector<Double> & coeffdamp);

  //! list of surfaces which make up the bnd. of the domain
  std::vector<std::string> domain_surfaces_;

  //! list of surfaces, on which we have force
  std::vector<std::string> rhs_surfaces_;

  //! function for RHS
  Integer arg_rhs_;

  //!
  std::vector<Double> directionFnc_;
  
  //! Indicator: is there RHS function
  Boolean SetRHSFnc;

  //! Indicator: is there RHS flow source
  Boolean SetRHSFlowSrc;

};

} // end of namespace
#endif
