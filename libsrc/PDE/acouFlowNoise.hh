#ifndef FILE_ACOUFLOWNOISE_2003
#define FILE_ACOUFLOWNOISE_2003

#include "basepde.hh"
#include "acousticPDE.hh"
#include <MpCCIcpl/MpCCIexch.hh>

namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class AcousticPDE. It is used for solving the acoustic analogy inhomogeneous acoustic equation for flow induced noise on one time step.  We define the right hand side by setting the sound sources using the results from the external fluid computation, set rules for assembling global system matrix according to weak form of this PDE, define right hand side and set first order absorbing boundary conditions. Then we call one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

class AcouFlowNoise: public AcousticPDE
{
public:

  //!
  AcouFlowNoise(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, 
		  FileType *aptFileType, WriteResults *aptOut);

  //!
  virtual ~AcouFlowNoise();

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


private:

  std::vector<std::string> rhs_surfaces_; //!< list of surfaces, on which we have excitation
  Integer arg_rhs_; //!< function for RHS

  Boolean SetRHSFnc; //!< Indicator: is there RHS function
  Boolean SetRHSFlowSrc; //!< Indicator: is there RHS flow source
  
  //Flow Data
  Matrix<Double> flowdata_;



  //!MpCCI
#ifdef MpCCI
  MpCCIexch * ptMpCCIexch_;
#endif
  Integer MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
  Integer MpCCI_; //<! if TRUE: coupling via MpCCI to low simulator
  Integer meshId_;
  Integer partId_;
  Integer nNodeIds_;
  Integer *nodeIds_;
  Integer GlobalDim_;
  Integer nElemIds_;
  Integer *elemIds_;
  Integer nElemTypes_;
  Integer *nNodesPerElem_;
  Integer *elemTypes_;
  Integer MpCCIprocess_;

};

} // end of namespace
#endif
