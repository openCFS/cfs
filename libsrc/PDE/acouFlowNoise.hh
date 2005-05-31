#ifndef FILE_ACOUFLOWNOISE_2003
#define FILE_ACOUFLOWNOISE_2003

#include "basePDE.hh"
#include "acousticPDE.hh"
#include <MpCCIcpl/MpCCIexch.hh>

namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class AcousticPDE. It is used for solving 
    the acoustic analogy inhomogeneous acoustic equation for flow induced
    noise on one time step. We define the right hand side by setting the
    sound sources using the results from the external fluid computation, set 
    rules for assembling global system matrix according to weak form of
    this PDE, define right hand side and set first order absorbing boundary
    conditions. Then we call one of methods of LinSystem for solving linear
    system. On the last step we calculate first and second derivatives of
    the solution.
  */

  class AcouFlowNoise: public AcousticPDE
  {
  public:

    //!
    AcouFlowNoise(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut);

    //!
    virtual ~AcouFlowNoise();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //!
    void ComputeRHS(const Double atime);

    //! Reads at every time the flowdatafile from the Fluid's Computation
    void ReadFlowData(const char * aname, const UInt timestep,
                      Matrix<Double> &nodedata );


    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteResultsInFile(const UInt kstep,
                                    const Double asteptime,
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0);
  private:


    UInt arg_rhs_; //!< function for RHS


    //Flow Data
    Matrix<Double> flowdata_;
    //!< name of subdomain to be coupled with MpCCI
    StdVector<RegionIdType> couplSubDomId_; 

    //!type of FlowData
    Boolean nodalSrc_;


    //  BaseNodeStoreSol * rhs_;        //!< For eventual saving of rhs as sol
    
    //!MpCCI
#ifdef MpCCI
    StdVector<UInt> mapSD_;
    MpCCIexch * ptMpCCIexch_;
    UInt MpCCInodes_; //<! number of FE-nodes for MpCCI-domain
    UInt MpCCI_; //<! if TRUE: coupling via MpCCI to low simulator
    UInt meshId_;
    UInt partId_;
    UInt nNodeIds_;
    UInt *nodeIds_;
    UInt GlobalDim_;
    UInt nElemIds_;
    UInt *elemIds_;
    UInt nElemTypes_;
    UInt *nNodesPerElem_;
    UInt *elemTypes_;
    UInt MpCCIprocess_;
#endif

  };

} // end of namespace
#endif
