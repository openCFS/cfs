#ifndef FILE_ACOU2DFLOWNOISE_2003
#define FILE_ACOU2DFLOWNOISE_2003

#include "basepde.hh"
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
  Acou2dFlowNoise(Grid * aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
		  FileType *aptFileType, WriteResults *aptOut);

  //!
  virtual ~Acou2dFlowNoise();

  //!
  void ComputeRHS(const Double atime);

  //!
  void preComputeRHS();

  //! Reads at every time the flowdatafile from the Fluid's Computation
  void ReadFlowData(const char * aname, const Integer timestep, Matrix<Double> &nodedata );

  //!
  void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

  //!
   void WriteResultsInFile();

  //! Reorganizing grid info for MpCCi and hand over to MpCCI
  void PutExchangeGrid2MpCCI();

private:

  std::vector<std::string> rhs_surfaces_; //!< list of surfaces, on which we have excitation
  Integer arg_rhs_; //!< function for RHS

  Boolean SetRHSFnc; //!< Indicator: is there RHS function
  Boolean SetRHSFlowSrc; //!< Indicator: is there RHS flow source

};

} // end of namespace
#endif
